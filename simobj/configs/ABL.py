from simobj.configs._EAGLE_toolkit import header_mask, group_mask, fof_mask, \
    id_mask, particle_mask_fofsub, particle_mask_fof, particle_mask_aperture, \
    particle_mask_pyread_eagle


# define suffix mnemonics for EAGLE particle types
T = dict(g=0, dm=1, s=4, bh=5)

# keys for recentering
recenter = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T.keys()
}

# keys for coordinate differentiation
coord_type = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'position'), ('vxyz', 'velocity')] for t in T.keys()
}

# keys for box wrapping
box_wrap = {'xyz_' + t: 'Lbox' for t in T.keys()}

# adjustments for aperture masking
extractor_edits = [
    (
        lambda E, A:
        (E.filetype == 'particle'),
        'filetype',
        'snapshot'
    )
]

# define masks
masks = dict()

masks['meta'] = header_mask
masks['group'] = group_mask
masks['fofgroup'] = fof_mask
masks['idgroup'] = id_mask
for ptype, pnum in T.items():
    masks['particle{:d}'.format(pnum)] = {
        'fofsub': particle_mask_fofsub(ptype),
        'fof': particle_mask_fof(ptype),
        'aperture': particle_mask_aperture(ptype),
        'pyread_eagle': particle_mask_pyread_eagle(ptype)
    }
