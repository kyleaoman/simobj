import numpy as np
import os
from importlib.util import spec_from_file_location, module_from_spec
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy import units as U
from simfiles import SimFiles
from ._L_align import L_align


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

    lowers = np.argwhere(np.diff(mask.astype(np.int)) > 0) + 1
    uppers = np.argwhere(np.diff(mask.astype(np.int)) < 0) + 1
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
    return


def apply_recenter(coords, centre):
    coords -= centre
    return


def apply_rotmat(coords, rotmat):
    coords = coords.dot(rotmat)
    return coords


def do_recenter(func):
    def rcfunc_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._recenter.keys():
            self._F.load(keys=(self._recenter[key], ),
                         verbose=self.init_args['verbose'])
            centres = self._F[self._recenter[key]]
            centre_mask = self._masks[
                self._F._extractors[self._recenter[key]].keytype]
            centre = centres[centre_mask].reshape((1, 3))
            centre -= self.current_translation[self._coord_type[key]]
            apply_recenter(self[key], centre)
            del self._F[self._recenter[key]]
        return self[key]
    return rcfunc_wrapper


def do_box_wrap(func):
    def wrapfunc_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._box_wrap.keys():
            self._F.load(keys=(self._box_wrap[key], ),
                         verbose=self.init_args['verbose'])
            Lbox = self._F[self._box_wrap[key]]
            apply_box_wrap(self[key], Lbox)
            del self._F[self._box_wrap[key]]
        return self[key]
    return wrapfunc_wrapper


def do_rotate(func):
    def rotfunc_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._coord_type.keys():
            self[key] = apply_rotmat(self[key], self.current_rot)
        return self[key]
    return rotfunc_wrapper


class MaskDict(dict):

    def __init__(self, SO):
        self.SO = SO
        return

    def __missing__(self, key):
        if key not in self.SO._maskfuncs.keys():
            raise KeyError
        value = self[key] = self.SO._maskfuncs[key](
            *self.SO.init_args['mask_args'],
            **(dict(
                {'vals': self.SO._F, 'verbose': self.SO.init_args['verbose']},
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
            ncpu=2
    ):
        self.init_args = dict()
        self.init_args['obj_id'] = obj_id
        self.init_args['snap_id'] = snap_id
        self.init_args['mask_type'] = mask_type
        self.init_args['mask_args'] = \
            tuple() if mask_args is None else mask_args
        self.init_args['mask_kwargs'] = \
            dict() if mask_kwargs is None else mask_kwargs
        self.init_args['configfile'] = configfile
        if (simfiles_configfile is not None) \
           and (simfiles_instance is not None):
            raise ValueError('Provide either simfiles_configfile or'
                             ' simfiles_instance, not both.')
        self.init_args['simfiles_configfile'] = simfiles_configfile
        self.init_args['verbose'] = verbose
        self.init_args['ncpu'] = ncpu
        self.current_rot = np.eye(3)
        self.current_translation = {
            'position': np.zeros(3) * U.kpc,
            'velocity': np.zeros(3) * U.km / U.s
        }

        self._read_config()

        if simfiles_configfile is not None:
            self._F = SimFiles(
                self.init_args['snap_id'],
                configfile=self.init_args['simfiles_configfile'],
                ncpu=self.init_args['ncpu']
            )
        elif simfiles_instance is not None:
            self._F = simfiles_instance
        else:
            raise ValueError('One of simfiles_configfile or simfiles_instance'
                             ' is required.')
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

        try:
            self._recenter = config.recenter
        except AttributeError:
            self._recenter = dict()

        try:
            self._coord_type = config.coord_type
        except AttributeError:
            self._coord_type = dict()

        try:
            self._box_wrap = config.box_wrap
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

    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    def __getattr__(self, key):
        if '__' in key:
            # avoid requesting reserved keys from SimFiles
            raise AttributeError
        try:
            return self.__dict__[key]
        except KeyError:
            return self[key]

    def __missing__(self, key):
        value = self[key] = self._load_key(key)
        return value

    @do_box_wrap
    @do_recenter
    @do_rotate
    def _load_key(self, key):
        if key not in self._F.fields():
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
                    grouping_ratio=1
                )
            elif not mask.any():
                intervals = ((0, 0), )
            else:
                intervals = mask_to_intervals(mask, grouping_ratio=1)
            parts = []
            for interval in intervals:
                self._F.load((key, ), intervals=(interval, ),
                             verbose=self.init_args['verbose'])
                if isinstance(mask, slice):
                    parts.append(self._F[key])
                else:
                    parts.append(self._F[key][mask[interval[0]: interval[1]]])
                del self._F[key]
            self[key] = np.concatenate([part.value for part in parts]) * \
                parts[0].unit

        elif self._F.share_mode:
            self._F.load((key, ), verbose=self.init_args['verbose'])
            self[key] = self._F[key][mask]
            # del disabled for share_mode

        else:
            self._F.load((key, ), verbose=self.init_args['verbose'])
            self[key] = self._F[key]
            del self._F[key]

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
        self.current_rot = apply_rotmat(self.current_rot, do_rot)
        return

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
        self.current_translation[translation_type] += translation.reshape((3,))
        return

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
