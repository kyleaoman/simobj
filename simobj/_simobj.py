from simfiles import SimFiles
from kyleaoman_utilities.slvars import savevars, loadvars
import numpy as np
import os

#Use (import) SimObj *not* _SimObj.

def usevals(names):
    def usevals_decorator(func):
        def func_wrapper(*args, **kwargs):
            loaded_keys = set()
            loaded_keys.update(args[0].load(names))
            retval = func(*args, **kwargs)
            for k in loaded_keys:
                del args[0][k]
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

def do_recenter(func):
    def func_wrapper(self, key):
        self[key] = func(self, key)
        if key in self._recenter.keys():
            self._F.load(keys=(self._recenter[key], ))
            centres = self._F[self._recenter[key]]
            centre_mask = self._masks[self._F._extractors[self._recenter[key]].keytype]
            centre = self._use_mask(centres, centre_mask)
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

class _SimObj(dict):
    def __init__(self, obj_id=None, snap_id=None, mask_type=None, mask_args=None, mask_kwargs=None, configfile=None, simfiles_configfile=None, cache_prefix='./', disable_cache=False):
        
        self.init_args = dict()
        self.init_args['obj_id'] = obj_id
        self.init_args['snap_id'] = snap_id
        self.init_args['mask_type'] = mask_type
        self.init_args['mask_args'] = tuple() if mask_args is None else mask_args
        self.init_args['mask_kwargs'] = dict() if mask_kwargs is None else mask_kwargs
        self.init_args['configfile'] = configfile
        self.init_args['simfiles_configfile'] = simfiles_configfile
        self.init_args['cache_prefix'] = cache_prefix
        self.init_args['disable_cache'] = disable_cache

        self._locked = False
        
        self._read_config()
        self._path = self.init_args['cache_prefix'] + '/' + 'SimObjCache_' + self._cache_string(**self.init_args)

        if not self.init_args['disable_cache']:
            if os.path.exists(self._path + '.lock'):
                raise RuntimeError(self._path + ".pkl' is locked by another instance.")
            else:
                self._lock()

        if os.path.exists(self._path + '.pkl') and not self.init_args['disable_cache']:
            D, = loadvars(self._path)
            self.update(D)
            self._read_config()
            self._F._read_config()
            self._edit_extractors()

        else:
            self._F = SimFiles(self.init_args['snap_id'], configfile=self.init_args['simfiles_configfile'])
            self._edit_extractors()
            self._init_masks()
            if not self.init_args['disable_cache']:
                self._cache()

        return

    def _read_config(self):
        
        config = dict()
        try:
            execfile(self.init_args['configfile'], config)
        except IOError:
            raise IOError("SimObj: configfile '" + self.init_args['configfile'] + "' not found.")

        try:
            self._cache_string = config['cache_string']
        except KeyError:
            raise valueError("SimObj: configfile missing 'cache_string' definition.")
            
        try:
            self._recenter = config['recenter']
        except KeyError:
            self._recenter = dict()

        try:
            self._box_wrap = config['box_wrap']
        except KeyError:
            self._box_wrap = dict()

        try:
            self._extractor_edits = config['extractor_edits']
        except KeyError:
            self._extractor_edits = list()

        try:
            config['masks']
        except KeyError:
            raise ValueError("SimObj: configfile missing 'masks' definition.")
        self._maskfuncs = dict()
        for key, maskfunc in config['masks'].items():
            self._maskfuncs[key] = maskfunc[self.init_args['mask_type']] if type(maskfunc) is dict else maskfunc

    def _init_masks(self):
        self._masks = dict()
        for key, maskfunc in self._maskfuncs.items():
            self._masks[key] = maskfunc(self._F, *self.init_args['mask_args'], **self.init_args['mask_kwargs'])
        return

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
    
    @do_recenter
    @do_box_wrap
    def _load_key(self, key):
        
        if key not in self._F.fields():
            raise KeyError("SimObj: SimFiles member unaware of '"+key+"' key.")

        self._F.load((key, ))
        mask = self._masks[self._F._extractors[key].keytype]
        self[key] = self._use_mask(self._F[key], mask)
            
        del self._F[key]

        if not self.init_args['disable_cache']:
            self._cache()

        return self[key]

    def _use_mask(self, data, mask):
        return data[mask] if mask is not None else data

    def _cache(self):
        if not self._locked:
            raise RuntimeError("SimObj does not own lock on cache (is it being used outside a 'with ... as ...' block?).")
        del self['_maskfuncs'], self['_cache_string'], self['_extractor_edits'], self._F['_extractors'], self._F['_snapshot']
        savevars([self], self._path + '.pkl')
        self._read_config()
        self._F._read_config()
        self._edit_extractors()
        return

    def _lock(self):
        self._locked = True
        open(self._path + '.lock', 'a').close()
        return

    def _unlock(self):
        self._locked = False
        os.remove(self._path + '.lock')
        return

    def _edit_extractors(self):
        for condition, field, value in self._extractor_edits:
            for key, extractor in self._F._extractors.items():
                if condition(extractor, self.init_args):
                    self._F._extractors[key] = extractor._replace(**{field: value})
        return

class SimObj:

    def __enter__(self):
        self._SO = _SimObj(*self._args, **self._kwargs)
        return self._SO

    def __exit__(self, exc_type, exc_value, traceback):
        if self._SO._locked:
            self._SO._unlock()
        return

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        return
