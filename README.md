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
   - ensure `~/lib/python2.7/site-packages` or similar is on your `PYTHONPATH` (e.g. `echo $PYTHONPATH`), if not, add it (perhaps in `.bash_profile` or similar)
   - run `pip install --prefix ~ -e .` (`-e` installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch `python`; you should be able to do `from simobj import SimObj`
 
 ## Usage:

Ensure you have [`simfiles`](https://github.com/kyleaoman/simfiles) properly installed and configured for your simulation before beginning.

```python
from simobj import SimObj
with SimObj(obj_id=None, snap_id=None, mask_type=None, mask_args=None, mask_kwargs=None, configfile=None, simfiles_configfile=None, cache_prefix='./', disable_cache=False)) as SO:
```

Initializes a `SimObj` object, representing an object/galaxy as defined in the simulation configured using the configfiles. A context manager (`with ... as ...:`) should always be used; failing to do so will usually result in an error and a reminder. By default, a cache will be written to speed up repeated reads of the same fields. To guard against concurrent writes to the cache file, a `.lock` file is created in the same directory as the cache file. Multiple scripts attempting to write the same cache simultaneously will result in an error (for all but the first script, which will successfully acquire the lock). Usually the lock is released on code exit, but abnormal code termination (e.g. `SIGINT`) may prevent it from being cleaned up. It can be safely deleted if no running processes are interacting with the cache. The cache can be safely deleted as well. Note that the cache is not in general portable across systems.
 - `obj_id`: (default: `None`) unique (within a snapshot) object identifier. This will be passed to user-defined functions in the `configfile`; these functions should be written with the desired `obj_id` format in mind.
 - `snap_id`: (default: `None`) unique snapshot identifier, which must be recognizable to a `SimFiles` object initialized with the `simfiles_configfile` supplied.
 - `mask_type`: (default: `None`) configuration may specify multiple types of masks, specify which to use.
 - `mask_args`: (default: `None`) if mask selected takes arguments in addition to those provided by default, provide them here. Prefer keyword arguments where possible.
 - `mask_kwargs`: (default: `None`) if mask selected takes additional keyword arguments, provide them here as a `dict`.
 - `configfile`: (default: `None`) path to a configfile; see [example configuration](https://github.com/kyleaoman/simobj/blob/master/simobj/configs/example.py) for a description of what should be defined in this file.
 - `simfiles_configfile`: (default: `None`) path to a configfile for [`simfiles`](https://github.com/kyleaoman/simfiles); see the [example configuration](https://github.com/kyleaoman/simfiles/blob/master/simfiles/configs/example.py) for that package for a description of what should be defined in this file.
 - `cache_prefix`: (default: `./`) location to store a cache file ([`.pkl` format](https://docs.python.org/2/library/pickle.html)) to (greatly) speed up loading values which have previously been loaded for this object, including across multiple scripts. See [example configuration](https://github.com/kyleaoman/simobj/blob/master/simobj/configs/example.py) for more information.
 - `disable_cache`: (default: `False`) set to `True` to disable caching behavior.
 
 ```python
     key1 = SO.key1
     key2 = SO.key2
     ...
 ```
Within the `with ... as ...` block, access to data is trivially simple via attribute syntax. All loading, dependencies, masking, etc. is handled internally based on information in the configfiles.
