from collections import namedtuple

def cache_string(**kwargs):
    snap_id = kwargs['snap_id']
    obj_id = kwargs['obj_id']
    return '_'.join((
        str(snap_id.res),
        str(snap_id.phys),
        str(snap_id.vol),
        str(snap_id.snap),
        str(obj_id.fof),
        str(obj_id.sub)
    ))

#define suffix mnemonics for EAGLE/APOSTLE particle types
T = ['g', 'dm', 'b2', 'b3', 's', 'bh']

recenter = {pos_vel[0] + '_' + t: pos_vel[1] for pos_vel in [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T}

box_wrap = {'xyz_' + t: 'Lbox' for t in T}

# The functions below are written with an obj_id such as the following in mind:
# obj_id = namedtuple('obj_id', ['fof', 'sub'])

extractor_edits = [
    (
        lambda E, A: ('particle' in E.keytype) and (A['mask_type'] == 'aperture'),
        'filetype',
        'snapshot'
    )
]

import numpy as np
from simobj import apply_recenter, apply_box_wrap, usevals

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

# For an aperture mask use snapshot files to ensure inclusion of *all* particles within aperture, this
# requires manual handling of loading/unloading keys instead of using '@usevals'.
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
