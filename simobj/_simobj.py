from simfiles import SimFiles
from utilities.slvars import savevars, loadvars
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
        if key in self.recenter.keys():
            self._F.load(keys=(self.recenter[key], ))
            centres = self._F[self.recenter[key]]
            centre_mask = self.masks[self._F._extractors[self.recenter[key]].keytype]
            centre = self._use_mask(centres, centre_mask)
            apply_recenter(self[key], centre)
            del self._F[self.recenter[key]]
        return self[key]
    return func_wrapper

def do_box_wrap(func):
    def func_wrapper(self, key):
        self[key] = func(self, key)
        if key in self.box_wrap.keys():
            self._F.load(keys=(self.box_wrap[key], ))
            Lbox = self._F[self.box_wrap[key]]
            apply_box_wrap(self[key], Lbox)
            del self._F[self.box_wrap[key]]
        return self[key]
    return func_wrapper

class _SimObj(dict):
    def __init__(self, obj_id, snap_id, mask_type=None, mask_args=None, mask_kwargs=None, configfile=None, simfiles_configfile=None, cache_prefix='./', disable_cache=False):
        
        self._locked = False
        self.disable_cache = disable_cache

        self._path = cache_prefix + '/' + 'SimObjCache_' #+ string conversion of snap_id, obj_id, mask_info?

        if not self.disable_cache:
            if os.path.exists(self._path + '.lock'):
                raise RuntimeError("SimObj '" + self._path + ".pkl' is locked by another instance.")
            else:
                self._lock()

        if os.path.exists(self._path + '.pkl') and not self.disable_cache:
            D, = loadvars(self._path)
            self.update(D)
            self._F._read_config()

        else:
            self.snap_id = snap_id
            self.obj_id = obj_id
            self.mask_type = mask_type
            self.mask_args = tuple() if mask_args is None else mask_args
            self.mask_kwargs = dict() if mask_kwargs is None else mask_kwargs
            self.configfile = configfile

            self._F = SimFiles(snap_id, configfile=simfiles_configfile)

            self._read_config()

            if not self.disable_cache:
                self._cache()

        return

    def _read_config(self):
        
        config = dict()
        try:
            execfile(self.configfile, config)
        except IOError:
            raise IOError("SimObj: configfile '" + self.configfile + "' not found.")

        try:
            self.recenter = config['recenter']
        except KeyError:
            self.recenter = dict()

        try:
            self.box_wrap = config['box_wrap']
        except KeyError:
            self.box_wrap = dict()

        try:
            config['masks']
        except KeyError:
            raise ValueError("SimObj: configfile missing 'masks' definition.")
        self.masks = dict()
        for key, value in config['masks'].items():
            if type(value) is dict:
                value = value[self.mask_type]
            self.masks[key] = value(self._F, *self.mask_args, **self.mask_kwargs)

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
        mask = self.masks[self._F._extractors[key].keytype]
        self[key] = self._use_mask(self._F[key], mask)
            
        del self._F[key]

        if not self.disable_cache:
            self._cache()

        return self[key]

    def _use_mask(self, data, mask):
        return data[mask] if mask is not None else data

    def _cache(self):
        if not self._locked:
            raise RuntimeError("SimObj does not own lock on cache (is it being used outside a 'with ... as ...' block?).")
        del self._F['_extractors'], self._F['_snapshot']
        savevars([self], self._path + '.pkl')
        self._F._read_config()
        return

    def _lock(self):
        self._locked = True
        open(self._path + '.lock', 'a').close()
        return

    def _unlock(self):
        self._locked = False
        os.remove(self._path + '.lock')
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
