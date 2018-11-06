from simfiles import SimFiles
import numpy as np
import os
from importlib.util import spec_from_file_location, module_from_spec
from kyleaoman_utilities.L_align import L_align
from astropy.coordinates.matrix_utilities import rotation_matrix


def mask_to_intervals(mask, grouping_ratio=0):
    '''
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
    '''
    lowers = np.argwhere(np.diff(mask.astype(np.int)) > 0) + 1
    uppers = np.argwhere(np.diff(mask.astype(np.int)) < 0) + 1
    if mask[0] is True:
        lowers = np.concatenate((np.array([[0]]), lowers))
    if mask[-1] is True:
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
            loaded_keys.update(kwargs['vals'].load(names))
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
    coords.dot(rotmat)
    return


def do_recenter(func):
    def func_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._recenter.keys():
            self._F.load(keys=(self._recenter[key], ))
            centres = self._F[self._recenter[key]]
            centre_mask = self._masks[
                self._F._extractors[self._recenter[key]].keytype]
            centre = centres[centre_mask]
            centre -= self.current_translation[self._coord_type[key]]
            apply_recenter(self[key], centre)
            del self._F[self._recenter[key]]
        return self[key]
    return func_wrapper


def do_box_wrap(func):
    def func_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._box_wrap.keys():
            self._F.load(keys=(self._box_wrap[key], ))
            Lbox = self._F[self._box_wrap[key]]
            apply_box_wrap(self[key], Lbox)
            del self._F[self._box_wrap[key]]
        return self[key]
    return func_wrapper


def do_rotate(func):
    def func_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._rotate.keys():
            self._F.load(keys=(self._box_wrap[key], ))
            rotmat = self.current_rot
            apply_rotmat(self[key], rotmat)
            del self._F[self._box_wrap[key]]
        return self[key]
    return func_wrapper


class MaskDict(dict):

    def __init__(self, SO):
        self.SO = SO
        return

    def __missing__(self, key):
        if key not in self.SO._maskfuncs.keys():
            raise KeyError
        value = self[key] = self.SO._maskfuncs[key](
            *self.SO.init_args['mask_args'],
            **(dict({'vals': self.SO._F}, **self.SO.init_args['mask_kwargs']))
        )
        return value


class SimObj(dict):

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
        self.init_args['ncpu'] = ncpu
        self.current_rot = np.eye(3)
        self.current_translation = {
            'position': np.zeros(3),
            'velocity': np.zeros(3)
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
                if type(maskfunc) is dict \
                else maskfunc

    def __setattr__(self, key, value):
        return self.__setitem__(key, value)

    def __getattr__(self, key):
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
        if (mask is not None) and (self._F.share_mode is False):
            if isinstance(mask, slice):
                intervals = ((mask.start, mask.stop), )
            elif not mask.any():
                intervals = ((0, 0), )
            else:
                intervals = mask_to_intervals(mask, grouping_ratio=1)
            parts = []
            for interval in intervals:
                self._F.load((key, ), intervals=(interval, ))
                if isinstance(mask, slice):
                    parts.append(self._F[key])
                else:
                    parts.append(self._F[key][mask[interval[0]: interval[1]]])
                del self._F[key]
            self[key] = np.concatenate([part.value for part in parts]) * \
                parts[0].unit

        elif self._F.share_mode is True:
            self._F.load((key, ))
            self[key] = self._F[key][mask]
            del self._F[key]

        else:
            self._F.load((key, ))
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

        do_rot = np.eye(3)

        if axis_angle is not None:
            do_rot = rotation_matrix(
                axis_angle[1],
                axis=axis_angle[0]
            ).dot(do_rot)

        if rotmat is not None:
            do_rot = rotmat.dot(do_rot)

        if L_coords is not None:
            mkey, xkey, vkey, incl, az_rot = L_coords
            do_rot = L_align(
                self[xkey],
                self[vkey],
                self[mkey],
                frac=.3,
                Laxis='z'
            ).dot(do_rot)

        keys = set(self.keys()).intersection(self._coord_type.keys())
        for key in keys:
            self[key].dot(rotmat)
        self.current_rot = self.current_rot.dot(rotmat)

    def translate(self, translation_type, translation):
        keys = set(self.keys()).intersection(
            {k: v for k, v in self._coord_type.items()
             if v == translation_type}
        )
        for key in keys:
            self[key] += translation
        self.current_translation[translation_type] += translation
        return

    def recenter(self, translation_type, new_centre):
        self.translate(translation_type, -new_centre)
        return
