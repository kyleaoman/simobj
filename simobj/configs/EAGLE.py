from simobj.configs._EAGLE_toolkit import header_mask, group_mask, fof_mask, \
    id_mask, particle_mask_fofsub, particle_mask_fof, particle_mask_aperture, \
    particle_mask_pyread_eagle


# define suffix mnemonics for EAGLE particle types
T = ['g', 'dm', 's', 'bh']

# keys for recentering
recenter = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
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
        ('particle' in E.keytype)
        and (A['mask_type'] in ('aperture', 'pyread_eagle')),
        'filetype',
        'snapshot'
    )
]

# define masks
masks = dict()

masks['header'] = header_mask
masks['group'] = group_mask
masks['fofgroup'] = fof_mask
masks['idgroup'] = id_mask
for ptype in T:
    masks['particle_'+ptype] = {
        'fofsub': particle_mask_fofsub(ptype),
        'fof': particle_mask_fof(ptype),
        'aperture': particle_mask_aperture(ptype),
        'pyread_eagle': particle_mask_pyread_eagle(ptype)
    }
