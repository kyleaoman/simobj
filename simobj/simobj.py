from simfiles import SimFiles
from utilities.slvars import savevars, loadvars
import numpy as np
import os

#Use (import) SimObj *not* _SimObj.

class _SimObj(dict):
    def __init__(self, snap_id, obj_id, mask_type=None, mask_args={}, cache_prefix='./', disable_cache=False, configfile=None):
        
        self._path = prefix + '/' + 'SimObjCache_' + #string conversion of snap_id, obj_id, mask_info??
        if not disable_cache:
            if os.path.exists(self._path + '.lock'):
                raise RuntimeError("SimObj '" + self._path + ".pkl' is locked by another instance.")
            else:
                self._lock()

        if os.path.exists(self._path + '.pkl') and not disable_cache:
            D, = loadvars(self._path)
            self.update(D)
            self._F._init_extractors()

        else:
            self.snap_id = snap_id
            self.obj_id = obj_id
            self.mask_type, self.mask_args = mask_type, mask_args

            self._F = SimFiles(snap_id, configfile=configfile)
            self._define_masks()
            if not diable_cache:
                self._cache()

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

    def _define_masks(self):
        loaded_keys = set()
        loaded_keys.update(self._F.load('group', ('gns', 'sgns', 'nfof', 'nID', 'offID')))
        self.gmask = np.logical_and(self._F.gns == self.fof, self._F.sgns == self.sub)
        self.fofmask = np.arange(1, self._F.nfof + 1) == self.fof
        self.idmask = np.s_[self._F.offID[np.logical_and(self._F.gns == self.fof, self._F.sgns == self.sub)][0] : self._F.offID[np.logical_and(self._F.gns == self.fof, self._F.sgns == self.sub)][0] + self._F.nID[np.logical_and(self._F.gns == self.fof, self._F.sgns == self.sub)][0]]
        self.pmasks = {}
        if self.mask_type == 'fof_sub':
            types = ['g', 'dm', 's', 'bh']
            loaded_keys.update(self._F.load('particle', tuple([field + typesuffix for field in ['ng_', 'nsg_'] for typesuffix in types])))
            for typesuffix in types:
                self.pmasks[typesuffix] = np.logical_and(self._F['ng_' + typesuffix] == self.fof, self._F['nsg_' + typesuffix] == self.sub)
        elif self.mask_type == 'fof':
            types = ['g', 'dm', 's', 'bh']
            loaded_keys.update(self._F.load('particle', tuple(['ng_' + typesuffix for typesuffix in types])))
            for typesuffix in types:
                self.pmasks[typesuffix] = self._F['ng_' + typesuffix] == self.fof
        elif self.mask_type == 'aperture':
            loaded_keys.update(self._F.load('group', ('cops', 'vcents')))
            loaded_keys.update(self._F.load('snapshot', ('xyz_g', 'xyz_dm', 'xyz_b2', 'xyz_b3', 'xyz_s', 'xyz_bh', 'Lbox')))
            for typesuffix in T.keys():
                self._F['xyz_' + typesuffix] = self._F['xyz_' + typesuffix] - self._F.cops[self.gmask]
                self._F['xyz_' + typesuffix][self._F['xyz_' + typesuffix] > self._F.Lbox / 2.] -= self._F.Lbox
                self._F['xyz_' + typesuffix][self._F['xyz_' + typesuffix] < self._F.Lbox / 2.] += self._F.Lbox
                cube = (np.abs(self._F['xyz_' + typesuffix]) < self.aperture).all(axis=1)
                self.pmasks[typesuffix] = np.zeros(self._F['xyz_' + typesuffix].shape[0], dtype=np.bool)
                self.pmasks[typesuffix][cube] = np.sum(np.power(self._F['xyz_' + typesuffix][cube], 2), axis=1) < np.power(self.aperture, 2)
        for k in loaded_keys:
            del self._F[k]
        return

    def _load_key(self, key):

        if key not in self._F.fields():
            raise KeyError

        loaded_keys = set()

        loaded_keys.update(self._F.load((key, )))
        self[key] = self._F[key][self.masks[self._F._extractors[key].keytype]]

        if key.split('_')[0] == 'xyz':
            loaded_keys.update(self._F.load('Lbox', 'cops'))
            self[key] -= self._F.cops[self.masks[self._F._extractors['cops'].keytype]]
            self[key][self[key] > self._F.Lbox / 2.] -= self._F.Lbox
            self[key][self[key] < -self._F.Lbox / 2.] += self._F.Lbox
        elif key.split('_')[0] == 'vxyz':
            loaded_keys.update(self._F.load(('vcents',)))
            self[key] -= self._F.vcents[self.masks[self._F._extractors['vcents'].keytype]]

        for k in loaded_keys:
            del self._F[k]

        self._cache()

        return self[key]

    def _cache(self):
        if not self._locked:
            raise RuntimeError("SimObj does not own lock on cache (is it being used outside a 'with ... as ...' block?).")
        del self._F['_extractors'], self._F['_snapshots'], self._F['_snapshot']
        savevars([self], self._path + '.pkl')
        self._F._read_config(self._F.snap_id, self._F.configfile)
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
        self._SO._unlock()
        return

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        return
