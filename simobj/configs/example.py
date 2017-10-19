from collections import namedtuple

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
        

#define object unique id tuple format
obj_id = namedtuple('obj_id', ['fof', 'sub'])

T = ['g', 'dm', 'b2', 'b3', 's', 'bh']

recenter = {pos_vel[0] + '_' + t: pos_vel[1] for pos_vel in [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T}

box_wrap = {'xyz_' + t: 'Lbox' for t in T}

import numpy as np
from simobj._simobj import apply_recenter, apply_box_wrap

masks = {}

@usevals(tuple())
def header_mask(vals, obj_id, aperture=None):
    return None

@usevals(('gns', 'sgns'))
def group_mask(vals, obj_id, aperture=None):
    return np.logical_and(vals.gns == obj_id.fof, vals.sgns == obj_id.sub)

@usevals(('nfof', ))
def fof_mask(vals, obj_id, aperture=None):
    return np.arange(1, vals.nfof + 1) == obj_id.fof

@usevals(('offID', 'nID'))
def id_mask(vals, obj_id, aperture=None):
    gmask = group_mask(vals, obj_id)
    start = vals.offID[gmask]
    end = start + vals.nID[gmask]
    return np.s_[start:end]

def particle_mask_fofsub(ptype):
    if ptype in ['b2', 'b3']:
        return lambda vals, obj_id: None
    @usevals(('ng_'+ptype, 'nsg_'+ptype))
    def mask(vals, obj_id, aperture=None):
        return np.logical_and(vals['ng_'+ptype] == obj_id.fof, vals['nsg_'+ptype] == obj_id.sub)
    return mask

def particle_mask_fof(ptype):
    if ptype in ['b2', 'b3']:
        return lambda vals, obj_id: None
    @usevals(('ng_'+ptype, ))
    def mask(vals, obj_id, aperture=None):
        return vals['ng_'+ptype] == obj_id.fof
    return mask

def particle_mask_aperture(ptype):
    def mask(vals, obj_id, aperture=None):
        key = 'xyz_'+ptype
        loaded_keys = set()
        loaded_keys.update(vals.load(keys=('cops', 'Lbox')))
        loaded_keys.update(vals.load(keys=(key, ), filetype='snapshot'))
        apply_recenter(vals[key], vals.cops[group_mask(vals, obj_id)])
        apply_box_wrap(vals[key], vals.Lbox)
        cube = (np.abs(vals[key]) < aperture).all(axis=1)
        retval = np.zeros(vals[key].shape[0], dtype=np.bool)
        retval[cube] = np.sum(np.power(vals[key][cube], 2), axis=1) < np.power(aperture, 2)
        #need to force reading from snapshot files, can't use extractors because cache will re-read config files and overwrite changes!
        for k in loaded_keys:
            del vals[k]
        return retval
    return mask

masks['header'] = header_mask
masks['group'] = group_mask
masks['fofgroup'] = fof_mask
masks['idgroup'] = id_mask
for ptype in T:
    masks['particle_'+ptype] = {
        'fofsub': particle_mask_fofsub(ptype),
        'fof': particle_mask_fof(ptype),
        'aperture': particle_mask_aperture(ptype)
    }
