import numpy as np
from simobj import usevals, apply_translate, apply_box_wrap

# define suffix mnemonics for EAGLE particle types
T = ['g', 'dm', 's', 'bh']

# keys for recentering
recenter = {
    pos_vel[0] + '_' + t: pos_vel[1] for pos_vel in
    [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T
}


# keys for coordinate differentiation
coord_type = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'position'), ('vxyz', 'velocity')] for t in T
}

# keys for box wrapping
box_wrap = {'xyz_' + t: 'Lbox' for t in T}

# adjustments for aperture masking
extractor_edits = [
    (
        lambda E, A:
        ('particle' in E.keytype) and (A['mask_type'] == 'aperture'),
        'filetype',
        'snapshot'
    )
]


# mask functions
def particle_mask_fof(ptype):
    @usevals(('ng_'+ptype, ))
    def mask(obj_id, vals=None, **kwargs):
        return np.abs(vals['ng_'+ptype]) == obj_id.fof
    return mask


def particle_mask_fofsub(ptype):
    @usevals(('ng_'+ptype, 'nsg_'+ptype))
    def mask(obj_id, vals=None, **kwargs):
        return np.logical_and(
            np.abs(vals['ng_'+ptype]) == obj_id.fof,
            vals['nsg_'+ptype] == obj_id.sub
        )
    return mask


def particle_mask_aperture(ptype):
    @usevals(('xyz_'+ptype, ))
    def mask(obj_id, vals=None, aperture=None, **kwargs):
        key = 'xyz_'+ptype
        gmask = group_mask(obj_id, vals=vals, **kwargs)
        apply_translate(vals[key], -vals['cops'][gmask])
        apply_box_wrap(vals[key], vals['Lbox'])
        retval = np.zeros(vals[key].shape[0], dtype=np.bool)
        outer_cube = (np.abs(vals[key]) < aperture).all(axis=1)
        inner_cube = np.zeros(vals[key].shape[0], dtype=np.bool)
        inner_cube[outer_cube] = (
            np.abs(vals[key][outer_cube]) < aperture / np.sqrt(3)
        ).all(axis=1)
        need_D = np.logical_and(outer_cube, np.logical_not(inner_cube))
        retval[inner_cube] = np.ones(np.sum(inner_cube))
        retval[need_D] = np.sum(np.power(vals[key][need_D], 2), axis=1) < \
            np.power(aperture, 2)
        return retval
    return mask


@usevals(('gns', 'sgns'))
def group_mask(obj_id, vals=None, **kwargs):
    return np.logical_and(vals.gns == obj_id.fof, vals.sgns == obj_id.sub)


@usevals(('nfof', ))
def fof_mask(obj_id, vals=None, **kwargs):
    return np.arange(1, vals.nfof + 1) == obj_id.fof


@usevals(('offID', 'nID'))
def id_mask(obj_id, vals=None, **kwargs):
    gmask = group_mask(obj_id, vals=vals)
    start = vals.offID[gmask][0]
    end = start + vals.nID[gmask][0]
    return np.s_[start:end]


@usevals(tuple())
def header_mask(*args, **kwargs):
    return None


# define masks
masks = {}

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
