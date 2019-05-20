import numpy as np
from simobj import usevals, apply_recenter, apply_box_wrap

# define suffix mnemonics for Auriga particle types
T = ['g', 'dm', 'b2', 'b3', 's', 'bh']

recenter = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1]
    for pos_vel in [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T
}

coord_type = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1]
    for pos_vel in [('xyz', 'position'), ('vxyz', 'velocity')] for t in T
}

box_wrap = {'xyz_' + t: 'Lbox' for t in T}


def particle_mask_fof(ptype):
    if ptype in ['b2', 'b3']:
        return lambda obj_id, vals=None: None

    @usevals(('l_'+ptype, ))
    def mask(obj_id, vals=None, **kwargs):
        fmask = fof_mask(obj_id, vals=vals, **kwargs)
        lasts = np.cumsum(vals['l_'+ptype], dtype=np.int).value
        firsts = np.r_[0, lasts[:-1]]
        return np.s_[firsts[fmask][0]: lasts[fmask][0]]

    return mask


def particle_mask_fofsub(ptype):
    if ptype in ['b2', 'b3']:
        return lambda obj_id, vals=None: None

    @usevals(('l_'+ptype, 'sl_'+ptype))
    def mask(obj_id, vals=None, **kwargs):
        fmask = fof_mask(obj_id, vals=vals, **kwargs)
        gmask = group_mask(obj_id, vals=vals, **kwargs)
        fof_offs = np.cumsum(vals['l_'+ptype], dtype=np.int).value
        fof_offs = np.r_[0, fof_offs[:-1]]
        fof_off = fof_offs[fmask][0]
        isub = np.argmax(gmask)
        sub_off = np.sum(
            vals['sl_'+ptype][isub - obj_id.sub:isub], dtype=np.int
        )
        sub_len = np.array(vals['sl_'+ptype], dtype=np.int)[gmask][0, 0]
        return np.s_[fof_off + sub_off:fof_off + sub_off + sub_len]

    return mask


def particle_mask_aperture(ptype):

    @usevals(('xyz_'+ptype, 'cops', 'Lbox'))
    def mask(obj_id, vals=None, aperture=None, **kwargs):
        key = 'xyz_'+ptype
        gmask = group_mask(obj_id, vals=vals, **kwargs)
        apply_recenter(vals[key], vals['cops'][gmask])
        apply_box_wrap(vals[key], vals['Lbox'])
        retval = np.zeros(vals[key].shape[0], dtype=np.bool)
        outer_cube = (np.abs(vals[key]) < aperture).all(axis=1)
        inner_cube = np.zeros(vals[key].shape[0], dtype=np.bool)
        inner_cube[outer_cube] = (np.abs(vals[key][outer_cube]) < aperture
                                  / np.sqrt(3)).all(axis=1)
        need_D = np.logical_and(outer_cube, np.logical_not(inner_cube))
        retval[inner_cube] = np.ones(np.sum(inner_cube))
        retval[need_D] = np.sum(np.power(vals[key][need_D], 2), axis=1) \
            < np.power(aperture, 2)
        return retval

    return mask


@usevals(('firstsub', 'nsubhalos'))
def group_mask(obj_id, vals=None, **kwargs):
    fmask = fof_mask(obj_id, vals=vals, **kwargs)
    if obj_id.sub >= vals.nsubhalos[fmask]:
        raise ValueError('sub exceeds number of subhalos.')
    return (
        np.array([(vals.firstsub[fmask] + obj_id.sub).value], dtype=np.int),
    )


@usevals(('nfof', ))
def fof_mask(obj_id, vals=None, **kwargs):
    return (np.array([obj_id.fof], dtype=np.int), )


@usevals(tuple())
def header_mask(obj_id, vals=None, **kwargs):
    return None


masks = {}

masks['header'] = header_mask
masks['group'] = group_mask
masks['fofgroup'] = fof_mask
for ptype in T:
    masks['particle_'+ptype] = {
        'fofsub': particle_mask_fofsub(ptype),
        'fof': particle_mask_fof(ptype),
        'aperture': particle_mask_aperture(ptype)
    }
