# simobj
Software abstraction of an object from a cosmological simulation (a "galaxy"), setup via a configuration file. Examine 'simobj/simobj/configs/example.py' for simple and advanced configuration examples. Note that you must first set up a configuration for simfiles (repository by kyleaoman) before using simobj.

**Installation:**
 - Download via web UI, or 'git clone https://github.com/kyleaoman/simobj.git'
 - Install dependencies if necessary (see 'setup.py'), some may be found in other repositories by kyleaoman.
 - Global install (Linux): 
   - cd to directory with 'setup.py'
   - run 'sudo pip install -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - User install (Linux):
   - cd to directory with 'setup.py'
   - ensure '~/lib/python2.7/site-packages' or similar is on your PYTHONPATH (e.g. 'echo $PYTHONPATH'), if not, add it (perhaps in .bash_profile or similar)
   - run 'pip install --prefix ~ -e .' (-e installs via symlink, so pulling repository will do a 'live' update of the installation)
 - cd to a directory outside the module and launch python; you should be able to do 'from simobj import SimObj'
