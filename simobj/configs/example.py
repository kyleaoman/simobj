import numpy as np
from simobj import usevals

# ------------------------------ OBJ_ID FORMAT --------------------------------
#
# The format of obj_id's does not need to be defined explicitly in this
# configfile, but many functions to be defined below must be written with the
# format of the obj_id in mind. In this example, the assumed format is:
#
#     from collections import namedtuple
#     obj_id = namedtuple('obj_id', ['fof', 'sub'])
#     my_obj = obj_id(fof = 1, sub = 0) #for instance
#
# The functions below may also make use of the snap_id. This ID *is* explicitly
# (user) defined, in the simfiles configfile. This example assumed the format
# from the simfiles example configuration:
#
#     from collections import namedtuple
#     snap_id = namedtuple('snap_id', ['res', 'phys', 'vol', 'snap'])
#     my_snap = snap_id(res = 3, phys = 'hydro', vol = 1, snap = 127)  # e.g.
#
# The examples below are intended for the APOSTLE simulations, used in
# conjunction with the configfile included as an example with the simfiles
# package.
#
# -----------------------------------------------------------------------------

# ------------------------ RECENTER AND BOX_WRAP ------------------------------
#
# In this portion two dicts may be defined, called 'recenter' and 'box_wrap'.
# If these are omitted the behavior described here simply won't occur.
#
# The keys of 'recenter' should be coordinate (position or velocity) keys for
# particle tables. The corresponding values should be keys pointing to the
# corresponding centroid (position or velocity). When loaded, these coordinates
# will have their centroid subtracted to recenter the coordinate system.
#
# The keys of 'recenter_type' should be coordinates (position or velocity) keys
# for particle tables. The corresponding values must be either 'position' or
# 'velocity'; this is used to determine which coordinate transformations should
# be applied to a given coordinate array.
#
# Similarly, the keys of 'box_wrap' should be coordinate (position) keys for
# particle tables in the case of a periodic box simulation. The corresponding
# values should be a key pointing to the box side length. Coordinate values
# which exceed half a box length (or are less than minus half a box length) in
# any direction will be 'wrapped'. This is crucial after changing the position
# centroid in a periodic box.
#
# -----------------------------------------------------------------------------

# define suffix mnemonics for EAGLE/APOSTLE particle types
T = ['g', 'dm', 'b2', 'b3', 's', 'bh']

recenter = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'cops'), ('vxyz', 'vcents')] for t in T
}

coord_type = {
    '{:s}_{:s}'.format(pos_vel[0], t): pos_vel[1] for pos_vel in
    [('xyz', 'position'), ('vxyz', 'velocity')] for t in T
}

box_wrap = {'xyz_' + t: 'Lbox' for t in T}

# -----------------------------------------------------------------------------

# ------------------------- EXTRACTOR EDITS -----------------------------------
#
# FOR ADVANCED USE ONLY. If a particular scenario requires modifying the values
# of simfiles extractor keys for certain particle types, a list of changes may
# be supplied as 'extractor_edits'. Each list entry should be a change
# described as 3-tuple. The first is a function which takes an extractor and
# the list of kwargs used to initialize the SimObj object. If this function
# returns true (all extractors will be tested), the next 2 entries will be used
# to modify the extractor. The second entry is the name of the key to modify,
# and the third is the new value for that key.
#
# In the example here, extractor_edits is used to force particle information to
# be read from the 'raw' snapshot files (i.e. pre-group finding) when the
# 'aperture' mask type is selected.
#
# -----------------------------------------------------------------------------


extractor_edits = [
    (
        lambda E, A:
        ('particle' in E.keytype) and (A['mask_type'] == 'aperture'),
        'filetype',
        'snapshot'
    )
]

# -----------------------------------------------------------------------------
# -------------------------- MASK FUNCTIONS -----------------------------------
#
# In this portion functions should be defined which allow the particles (and
# other properties) corresponding to a simulated object/galaxy to be selected
# from the raw tables as read by a SimFiles object. This example assumes that
# three mask_type's will be used:
#   1. A mask to select particles belonging to a particular fof group
#      (particle_mask_fof).
#   2. A mask to select particles belonging to a particular subfind subhalo
#      (particle_mask_fofsub).
#   3. A mask to select particles within a fixed distance of a point
#      (particle_mask_aperture).
# In addition, more functions are defined to enable selecting values
# corresponding to the group of interest from tables of:
#   4. Subhalo properties tables, e.g. centroids, etc. (group_mask).
#   5. FOF group properties tables, e.g. M200, R200, etc. (fof_mask).
#   6. Particle ID tables sorted by group (id_mask).
# Finally, a 'null mask' is defined for values that require no masking, e.g.
# constants:
#   7. Null mask (header_mask).
#
# The details of the function of course depend sensitively on the details of
# the simulation and/or group finding algorithm used. All masks, however, share
# some common properties:
#   - The function will receive as arguments and keyword arguments the values
#     of the parameters mask_args and mask_kwargs, respectively, used to
#     initialize the SimObj instance. The functions should be prepared to
#     receive these. In this example, one mask_arg is exepcted: the obj_id. The
#     aperture mask also expects a keyword argument: aperture.
#   - Functions may make use of the @usevals(names) function decorator. This
#     should contain a list of any keys to extractors which the function may
#     make use of. These will be explicitly loaded using simfiles, then
#     unloaded at the end of the mask function call. If the decorator is used,
#     the function should accept one additional keyword argument 'vals=None'.
#     The data can then be accessed within the function as attributes or dict
#     entries of 'vals'. Explicitly no values may be specified as
#     @usevals(tuple()), though this is not required.
#   - In addition to the vals, any previously computed masks are also available
#     via a kwarg called 'masks'. For instance, put masks=None in the signature
#     and then gmask = masks['group'] in the body to access the group mask. The
#     mask will be computed if not yet available, otherwise just accessed.
#     Recursion is to be avoided, however!
#   - Functions may either return a boolean array which will be used to mask
#     tables loaded via simfiles, or None, in which case no mask will be
#     applied.
#
# In the example below, in some cases a function is used to generate several
# similar mask functions, e.g. for different particle types.
#
# -----------------------------------------------------------------------------


def particle_mask_fof(ptype):
    if ptype in ['b2', 'b3']:
        return lambda obj_id, vals=None, **kwargs: None

    @usevals(('ng_'+ptype, ))
    def mask(obj_id, vals=None, **kwargs):
        return np.abs(vals['ng_'+ptype]) == obj_id.fof
    return mask


def particle_mask_fofsub(ptype):
    if ptype in ['b2', 'b3']:
        return lambda obj_id, vals=None, **kwargs: None

    @usevals(('ng_'+ptype, 'nsg_'+ptype))
    def mask(obj_id, vals=None, **kwargs):
        return np.logical_and(
            np.abs(vals['ng_'+ptype]) == obj_id.fof,
            vals['nsg_'+ptype] == obj_id.sub
        )
    return mask


def particle_mask_aperture(ptype):

    @usevals(('xyz_'+ptype, 'cops', 'Lbox'))
    def mask(obj_id, vals=None, aperture=None, masks=None, **kwargs):
        from simobj import apply_translate, apply_box_wrap
        key = 'xyz_'+ptype
        gmask = masks['group']
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
def id_mask(obj_id, vals=None, masks=None, **kwargs):
    gmask = masks['group']
    start = vals.offID[gmask][0]
    end = start + vals.nID[gmask][0]
    start, end = int(start.value), int(end.value)
    return np.s_[start:end]


@usevals(tuple())
def header_mask(obj_id, vals=None, **kwargs):
    return None

# -----------------------------------------------------------------------------

# ----------------------------- DEFINE MASKS ----------------------------------
#
# In this portion, the functions defined in the previous portion must be
# assigned to a dict called 'masks'. There keys of the dict should be the
# different keytype values known to simfiles (defined in the extractors). The
# values may be either of:
# - A mask function. This mask function will be applied to fields matching the
#   keytype regardless of the mask_type.
# - A dict, its keys should be the mask_type values intended to be used with
#   SimObj objects. Its values are then the corresponding mask functions.
#
# -----------------------------------------------------------------------------


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

# -----------------------------------------------------------------------------
