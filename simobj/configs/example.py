T = ['g', 'dm', 'b2', 'b3', 's', 'bh']

recenter = {pos_vel[0] + '_' + t: pos_vel[1] for pos_vel in [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T}

box_wrap = {'xyz_' + t: 'Lbox' for t in T}
