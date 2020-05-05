# Run all test scripts and doctests.

import subprocess

files_to_test = [
  'projection_tools.py',
  'pj_healpix.py',
  'pj_rhealpix.py',
  'ellipsoids.py',
  'projection_wrapper.py',
  'dggs_rhealpix.py',
  'distortion.py',
]

# Run test scripts.
for f in files_to_test:
    print 'Testing %s...' % f
    command = 'python ' + 'tests/test_' + f
    subprocess.call(command, shell=True)

# Run doctests. 
for f in files_to_test:
    print 'Doctesting %s...' % f
    command = 'python -m doctest ' + f
    subprocess.call(command, shell=True)
