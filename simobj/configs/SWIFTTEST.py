import numpy as np
from simobj import usevals

# define suffix mnemonics for EAGLE particle types
T = ['g', 'dm', 's']

# keys for coordinate differentiation
coord_type = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'position'), ('vxyz', 'velocity')] for t in T
}

# keys for box wrapping
box_wrap = {'xyz_' + t: 'Lbox' for t in T}


def particle_mask_cube(ptype):
    @usevals(('xyz_'+ptype, 'Lbox'))
    def mask(obj_id, vals=None, Xcube=None, Lcube=None, **kwargs):
        from simobj import apply_translate, apply_box_wrap
        key = 'xyz_'+ptype
        apply_translate(vals[key], -Xcube)
        apply_box_wrap(vals[key], vals['Lbox'])
        return np.logical_and(
            (vals[key] < Lcube / 2).all(axis=1),
            (vals[key] > -Lcube / 2).all(axis=1)
        )
    return mask


@usevals(tuple())
def header_mask(*args, **kwargs):
    return None


# define masks
masks = {}

masks['header'] = header_mask
for ptype in T:
    masks['particle_'+ptype] = {
        'cube': particle_mask_cube(ptype),
    }
