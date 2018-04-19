# simobj
Software abstraction of an object from a cosmological simulation (a "galaxy"), setup via a configuration file. Examine [`example.py`](https://github.com/kyleaoman/simobj/blob/master/simobj/configs/example.py) for simple and advanced configuration examples. Note that you must first set up a [configuration]() for [`simfiles`](https://github.com/kyleaoman/simfiles) before using [simobj](https://github.com/kyleaoman/simobj).

## Installation:
 - Download via web UI, or `git clone https://github.com/kyleaoman/simobj.git`
 - Install dependencies if necessary (see [`setup.py`](https://github.com/kyleaoman/simobj/blob/master/setup.py)), some may be found in [other repositories by kyleaoman](https://github.com/kyleaoman?tab=repositories)
 - Global install (Linux): 
   - `cd` to directory with [`setup.py`](https://github.com/kyleaoman/simobj/blob/master/setup.py)
   - run `sudo pip install -e .` (`-e` installs via symlink, so pulling repository will do a 'live' update of the installation)
 - User install (Linux):
   - `cd` to directory with [`setup.py`](https://github.com/kyleaoman/simobj/blob/master/setup.py)
   - ensure `~/lib/python3.6/site-packages` or similar is on your `PYTHONPATH` (e.g. `echo $PYTHONPATH`), if not, add it (perhaps in `.bash_profile` or similar)
   - run `pip install --prefix ~ -e .` (`-e` installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch `python`; you should be able to do `from simobj import SimObj`
 
 ## Usage:

Ensure you have [`simfiles`](https://github.com/kyleaoman/simfiles) properly installed and configured for your simulation before beginning.

```python
from simobj import SimObj
SO = SimObj(obj_id=None, snap_id=None, mask_type=None, mask_args=None, mask_kwargs=None, configfile=None, simfiles_configfile=None, simfiles_instance=None, ncpu=2):
```

Initializes a `SimObj` object, representing an object/galaxy as defined in the simulation configured using the configfiles.
 - `obj_id`: (default: `None`) unique (within a snapshot) object identifier. This will be passed to user-defined functions in the `configfile`; these functions should be written with the desired `obj_id` format in mind.
 - `snap_id`: (default: `None`) unique snapshot identifier, which must be recognizable to a `SimFiles` object initialized with the `simfiles_configfile` supplied, or the `simfiles_instance` supplied.
 - `mask_type`: (default: `None`) configuration may specify multiple types of masks, specify which to use.
 - `mask_args`: (default: `None`) if mask selected takes arguments in addition to those provided by default, provide them here. Prefer keyword arguments where possible.
 - `mask_kwargs`: (default: `None`) if mask selected takes additional keyword arguments, provide them here as a `dict`.
 - `configfile`: (default: `None`) path to a configfile; see [example configuration](https://github.com/kyleaoman/simobj/blob/master/simobj/configs/example.py) for a description of what should be defined in this file.
 - `simfiles_configfile`: (default: `None`) path to a configfile for [`simfiles`](https://github.com/kyleaoman/simfiles); see the [example configuration](https://github.com/kyleaoman/simfiles/blob/master/simfiles/configs/example.py) for that package for a description of what should be defined in this file. Specify this parameter, or `simfiles_instance`, not both.
 - `simfiles_instance`: (default: `None`) a [`SimFiles`](https://github.com/kyleaoman/simfiles) object. Specify this parameter, or `simfiles_configfile`, not both. Providing an object instead of a configuration is useful if use of the `share_mode` option of the SimFiles object is desired - when this mode is set, data persists within the SimFiles object, which will prevent redundant reading from disk by multiple SimObj instances corresponding to multiple simulation objects which are part of an underlying simulation fileset. In other words, this should be used when looping over many objects in one simulation volume. `share_mode` should be set `False` for use with a single SimObj, or a small number (containing a small fraction of the simulation particles), to enable a different set of speed optimizations for this scenario. Note that the SimFiles object can be "shared" only in a serial sense, by one SimObj instance after another; parallel access is not supported.
 - `ncpu`: (default: `2`) number of processes to use for parallel reading of data files. A value of 0 defaults to the number of CPUs less one. A value of 2 seemed optimal for one system I tested on, but this depends on hardware & configuration. A value of 1 (serial mode) is required if the SimFiles object will be used inside a parallel (portion of an) application.
 
 ```python
 key1 = SO.key1
 key2 = SO.key2
 key3 = SO['key3']
 ...
 ```
Access to data is trivially simple via attribute syntax (or getitem, i.e. []). All loading, dependencies, masking, etc. is handled internally based on information in the configfiles. SimObj is also compatible with the use of a context manager, i.e. `with SimObj(...) as SO:`.
