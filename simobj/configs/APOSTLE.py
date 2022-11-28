from simobj.configs._EAGLE_toolkit import header_mask, group_mask, fof_mask, \
    id_mask, particle_mask_fofsub, particle_mask_fof, particle_mask_aperture

# define suffix mnemonics for EAGLE particle types
T = dict(g=0, dm=1, b2=2, b3=3, s=4, bh=5)

recenter = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T.keys()
}

coord_type = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'position'), ('vxyz', 'velocity')] for t in T.keys()
}

box_wrap = {'xyz_' + t: 'Lbox' for t in T.keys()}

extractor_edits = [
    (
        lambda E, A:
        ('particle' in E.keytype) and (A['mask_type'] == 'aperture'),
        'filetype',
        'snapshot'
    ),
    (
        lambda E, A:
        ('particle' in E.keytype) and (A['mask_type'] != 'aperture'),
        'filetype',
        'particle'
    )
]

masks = dict()

masks['meta'] = header_mask
masks['group'] = group_mask
masks['fofgroup'] = fof_mask
masks['idgroup'] = id_mask
for ptype, pnum in T.items():
    masks['particle{:d}'.format(pnum)] = {
        'fofsub': particle_mask_fofsub(ptype),
        'fof': particle_mask_fof(ptype),
        'aperture': particle_mask_aperture(ptype)
    }
