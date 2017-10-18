T = ['g', 'dm', 'b2', 'b3', 's', 'bh']

recenter = {pos_vel[0] + '_' + t: pos_vel[1] for pos_vel in [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T}

box_wrap = {'xyz_' + t: 'Lbox' for t in T}

masks['header'] = 
masks['group'] = 
masks['fofgroup'] = 
masks['idgroup'] = 
for ptype in T:
    mask['particle_'+ptype] =

"""
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
"""
