import numpy as np
import os
from importlib.util import spec_from_file_location, module_from_spec
from astropy.coordinates.matrix_utilities import rotation_matrix
from simfiles import SimFiles
from ._L_align import L_align


def dealias(func):
    def dealias_wrapper(self, key, *args, **kwargs):
        if not key.startswith('_') and hasattr(self._F, '_aliases'):
            key = self._F._aliases.get(key, key)
        return func(self, key, *args, **kwargs)
    return dealias_wrapper


def mask_to_intervals(mask, grouping_ratio=0):
    """
    Convert a boolean mask to a set of 2-tuples defining index intervals.

    Takes a 1D boolean mask and converts it to a set of tuples delimiting the
    intervals where the mask is 'True'. If the grouping ratio is set > 0, gaps
    between 'True' intervals which are a smaller fraction of the mask array
    size than the grouping ratio will be concatenated. A grouping ratio of 1
    groups all the intervals, i.e. returns the indices of the first and last
    'True' values. This is intended for facilitating reading somewhat
    contiguous data from hdf5 files; grouping minimizes the number of reads,
    while avoiding reading entire arrays. A grouping ratio of .2 sets the
    maximum number of reads to 5, which I think will speed up many reads
    without undue slowdown in pathological cases, though I have not tested this
    empirically. If grouping is used, the chunks need to be masked (with the
    corresponding chunks of the mask) after reading.

    Parameters
    ----------
    mask : array-like
        Boolean mask to be converted to intervals.

    grouping_ratio : float
        How much to concatenate intervals with small gaps between them.

    Returns
    -------
    out : list
        List of 2-tuples defining intervals.
    """

    lowers = np.argwhere(np.diff(mask.astype(int)) > 0) + 1
    uppers = np.argwhere(np.diff(mask.astype(int)) < 0) + 1
    if bool(mask[0]) is True:
        lowers = np.concatenate((np.array([[0]]), lowers))
    if bool(mask[-1]) is True:
        uppers = np.concatenate((uppers, np.array([[mask.size]])))
    intervals = np.hstack((lowers, uppers))
    grouped = [tuple(intervals[0])]
    for interval in intervals[1:]:
        if (interval[0] - grouped[-1][1]) / mask.size < grouping_ratio:
            grouped[-1] = (grouped[-1][0], interval[1])
        else:
            grouped.append(tuple(interval))
    return grouped


def usevals(names):
    def usevals_decorator(func):
        def func_wrapper(*args, **kwargs):
            loaded_keys = set()
            loaded_keys.update(
                kwargs['vals'].load(names, verbose=kwargs['verbose'])
            )
            retval = func(*args, **kwargs)
            for k in loaded_keys:
                del kwargs['vals'][k]
            return retval
        return func_wrapper
    return usevals_decorator


def apply_box_wrap(coords, length):
    coords[coords > length / 2.] -= length
    coords[coords < -length / 2.] += length
    return coords


def apply_translate(coords, offset):
    coords += offset
    return coords


def apply_rotmat(coords, rotmat):
    coords = coords.dot(rotmat)
    return coords


def do_recenter(func):
    def rcfunc_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._recenter.keys():
            centre = self[self._recenter[key]]
            self[key] = apply_translate(self[key], -centre)
        return self[key]
    return rcfunc_wrapper


def do_box_wrap(func):
    def wrapfunc_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._box_wrap.keys():
            Lbox = self[self._box_wrap[key]]
            self[key] = apply_box_wrap(self[key], Lbox)
        return self[key]
    return wrapfunc_wrapper


def do_transform_stack(func):
    def tsfunc_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._coord_type.keys():
            for pop in self._transform_stack:
                if pop[0] == 'T' and pop[1] == self._coord_type[key]:
                    self[key] = apply_translate(self[key], pop[2])
                elif pop[0] == 'R':
                    self[key] = apply_rotmat(self[key], pop[1])
        return self[key]
    return tsfunc_wrapper


class MaskDict(dict):

    def __init__(self, SO):
        self.SO = SO
        return

    def __missing__(self, key):
        if key not in self.SO._maskfuncs.keys():
            print(key)
            raise KeyError
        value = self[key] = self.SO._maskfuncs[key](
            *self.SO.init_args['mask_args'],
            **(dict(
                {
                    'vals': self.SO._F,
                    'verbose': self.SO.init_args['verbose'],
                    'SO': self.SO
                },
                **self.SO.init_args['mask_kwargs']
            ))
        )
        return value


class SimObj(dict):
    """
    Code abstraction of an object from a cosmological simulation.

    SimObj is a dict with added features. It provides __getattr__ and
    __setattr__ for ease of access, and automatically loads data from files on
    disk as needed. It can automatically mask the loaded data to isolate
    individual "galaxies". It includes some optimizations for faster loading.

    Parameters
    ----------
    obj_id : index
        An identifier for a specific simulation object. The exact format is
        defined in the configuration for the simulation in question (default:
        None).

    snap_id : index
        An identifier for a specific simulation snapshot. The exact format is
        defined in the configuration for the simulation in question default:
        None).

    mask_type : str
        Key corresponding to one of the masks defined inthe configuration for
        the simulation in question.

    mask_args : tuple
        Any arguments that must be passed to the mask function (default:
        tuple()).

    mask_kwargs : dict
        Any keyword arguments to pass to the mask function (default: dict()).

    configfile : str
        Path to the configuration file to use (default: None).

    simfiles_configfile : str
        Path to the configuration file to use for simfiles. Provide either this
        or a simfiles_instance, not both (default: None).

    simfiles_instance : SimFiles
        Initialized instance of SimFiles class. Provide either this or a
        simfiles_configfile, not both. When initializing SimFiles, setting
        the share_mode flag is recommended (default: None).

    ncpu : int
        Number of processors on which to run (default: 2).
    """

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return

    def __init__(
            self,
            obj_id=None,
            snap_id=None,
            mask_type=None,
            mask_args=None,
            mask_kwargs=None,
            configfile=None,
            simfiles_configfile=None,
            simfiles_instance=None,
            verbose=False,
            ncpu=2,
            grouping_ratio=1,
            autorecenter_off=False
    ):
        if (simfiles_configfile is not None) \
           and (simfiles_instance is not None):
            raise ValueError('Provide either simfiles_configfile or'
                             ' simfiles_instance, not both.')
        if simfiles_configfile is not None:
            self._F = SimFiles(
                snap_id,
                configfile=simfiles_configfile,
                ncpu=ncpu
            )
        elif simfiles_instance is not None:
            self._F = simfiles_instance
        else:
            raise ValueError('One of simfiles_configfile or simfiles_instance'
                             ' is required.')
        self.init_args = dict()
        self.init_args['obj_id'] = obj_id
        self.init_args['snap_id'] = snap_id
        self.init_args['mask_type'] = mask_type
        self.init_args['mask_args'] = \
            tuple() if mask_args is None else mask_args
        self.init_args['mask_kwargs'] = \
            dict() if mask_kwargs is None else mask_kwargs
        self.init_args['configfile'] = configfile
        self.init_args['simfiles_configfile'] = simfiles_configfile
        self.init_args['verbose'] = verbose
        self.init_args['ncpu'] = ncpu
        self.init_args['grouping_ratio'] = grouping_ratio
        self.init_args['autorecenter_off'] = autorecenter_off
        self._transform_stack = list()

        self._read_config()
        self._edit_extractors()

        self._masks = MaskDict(self)

        return

    def _read_config(self):

        try:
            spec = spec_from_file_location('config', os.path.expanduser(
                self.init_args['configfile']))
            config = module_from_spec(spec)
            spec.loader.exec_module(config)
        except FileNotFoundError:
            raise FileNotFoundError(
                "SimObj: configfile '{:s}' not found.".format(
                    self.init_args['configfile']))

        if self.init_args['autorecenter_off']:
            self._recenter = dict()
        else:
            try:
                self._recenter = {self._F._aliases.get(k, k): v
                                  for k, v in config.recenter.items()}
            except AttributeError:
                self._recenter = dict()

        try:
            self._coord_type = {self._F._aliases.get(k, k): v
                                for k, v in config.coord_type.items()}
        except AttributeError:
            self._coord_type = dict()

        try:
            self._box_wrap = {self._F._aliases.get(k, k): v
                              for k, v in config.box_wrap.items()}
        except AttributeError:
            self._box_wrap = dict()

        try:
            self._extractor_edits = config.extractor_edits
        except AttributeError:
            self._extractor_edits = list()

        try:
            config.masks
        except AttributeError:
            raise ValueError("SimObj: configfile missing 'masks' definition.")
        self._maskfuncs = dict()
        for key, maskfunc in config.masks.items():
            self._maskfuncs[key] = maskfunc[self.init_args['mask_type']] \
                if isinstance(maskfunc, dict) \
                else maskfunc

    @dealias
    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    @dealias
    def __getattr__(self, key):
        if '__' in key:
            # avoid requesting reserved keys from SimFiles
            raise AttributeError
        try:
            return self.__dict__[key]
        except KeyError:
            return self[key]

    @dealias
    def __missing__(self, key):
        value = self[key] = self._load_key(key)
        return value

    @dealias
    @do_box_wrap
    @do_transform_stack
    @do_recenter
    def _load_key(self, key):
        if key not in set(self._F.fields(aliases=False)):
            raise KeyError("SimObj: SimFiles member unaware of '"+key+"' key.")

        mask = self._masks[self._F._extractors[key].keytype]
        if (mask is not None) and (not self._F.share_mode):
            if isinstance(mask, slice):
                intervals = ((mask.start, mask.stop), )
            elif isinstance(mask, tuple):
                boolmask = np.zeros(mask[0].max() + 1)
                boolmask[mask] = True
                intervals = mask_to_intervals(
                    boolmask,
                    grouping_ratio=self.init_args['grouping_ratio']
                )
            elif not mask.any():
                intervals = ((0, 0), )
            else:
                intervals = mask_to_intervals(
                    mask, grouping_ratio=self.init_args['grouping_ratio'])
            parts = []
            for interval in intervals:
                loaded_keys = self._F.load((key, ), intervals=(interval, ),
                                           verbose=self.init_args['verbose'])
                if isinstance(mask, slice):
                    parts.append(self._F[key])
                else:
                    parts.append(self._F[key][mask[interval[0]: interval[1]]])
                for k in loaded_keys:
                    del self._F[k]
            try:
                self[key] = np.concatenate([part.value for part in parts]) * \
                    parts[0].unit
            except AttributeError:
                self[key] = np.concatenate(parts)

        elif self._F.share_mode:
            self._F.load((key, ), verbose=self.init_args['verbose'])
            self[key] = self._F[key][mask]
            # del disabled for share_mode

        else:
            loaded_keys = self._F.load(
                (key, ),
                verbose=self.init_args['verbose']
            )
            self[key] = self._F[key]
            for k in loaded_keys:
                del self._F[k]

        return self[key]

    def _edit_extractors(self):
        for condition, field, value in self._extractor_edits:
            for key, extractor in self._F._extractors.items():
                if condition(extractor, self.init_args):
                    self._F._extractors[key] = extractor._replace(
                        **{field: value})
        return

    def rotate(self, axis_angle=None, rotmat=None, L_coords=None):
        """
        Rotate the object.

        If multiple kwargs are given, they are applied in the order of
        (axis_angle, rotmat, L_coords). Note that L_coords will override
        any previous rotations.

        Parameters
        ----------
        axis_angle : tuple
            2-tuple with an axis 'x', 'y' or 'z' and an astropy.units.Quantity
            instance with dimensions of angle, indicating the axis to rotate
            about and the angle to rotate through, respectively.

        rotmat : array-like
            A (3, 3) array specifying a rotation.

        L_coords : tuple
            A 3-tuple containing the keys for the mass, coordinates and
            velocities of the particles. The coordinate system is aligned to
            the angular momentum vector of the inner 1/3 of the particles.
            Empirically, this often results in a good alignment to the galactic
            disc. If finer control is needed, derive the desired rotation
            separately and provide a rotmat instead.

        Returns
        -------
        do_rot : np.ndarray
            The rotation matrix used.
        """

        do_rot = np.eye(3)

        if axis_angle is not None:
            do_rot = rotation_matrix(
                axis_angle[1],
                axis=axis_angle[0]
            ).dot(do_rot)

        if rotmat is not None:
            do_rot = rotmat.dot(do_rot)

        if L_coords is not None:
            mkey, xkey, vkey = L_coords
            do_rot = L_align(
                self[xkey],
                self[vkey],
                self[mkey],
                frac=.3,
                Laxis='z'
            ).dot(do_rot)

        keys = set(self.keys()).intersection(self._coord_type.keys())
        for key in keys:
            self[key] = apply_rotmat(self[key], do_rot)
        self._transform_stack.append(('R', do_rot))
        return do_rot

    def unrotate(self):
        """
        Reverse last coordinate transformation if it was a rotation.

        Returns
        -------
        do_rot : np.ndarray
            The rotation matrix used.
        """

        last_transform = self._transform_stack.pop()
        if last_transform[0] != 'R':
            self._transform_stack.append(last_transform)
            raise RuntimeError('Cannot unrotate if last transformation was not'
                               ' a rotation.')
        do_rot = last_transform[1].T
        keys = set(self.keys()).intersection(self._coord_type.keys())
        for key in keys:
            self[key] = apply_rotmat(self[key], do_rot)
        return do_rot

    def translate(self, translation_type, translation):
        """
        Translate all relevant properties of the object in either position or
        velocity.

        Parameters
        ----------
        translation_type : str
            One of the coord_types defined in the configuration.

        translation : astropy.units.Quantity instance
            Amount by which to translate, with compatible units.
        """
        keys = set(self.keys()).intersection(
            {k: v for k, v in self._coord_type.items()
             if v == translation_type}
        )
        for key in keys:
            self[key] += translation
        self._transform_stack.append(
            ('T', translation_type, translation)
        )
        return

    def transform(self, transform_stack):
        for tf in transform_stack:
            {'T': self.translate, 'R': self.rotate}[tf[0]](*tf[1:])

    def recenter(self, translation_type, new_centre):
        """
        Reposition the coordinate centre in either position or velocity.

        This is implemented as a translation by the negative of the new centre.

        Parameters
        ----------
        translation_type : str
            One of the coord_types defined in the configuration.

        new_centre : astropy.units.Quantity instance
            Location of the new centre, with compatible units.
        """

        self.translate(translation_type, -new_centre)
        return
