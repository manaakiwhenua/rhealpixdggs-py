# Run all test scripts and doctests.
import unittest, doctest, importlib

package_name = 'rhealpix_dggs'
module_names = [
  'utils',
  'pj_healpix',
  'pj_rhealpix',
  'ellipsoids',
  'projection_wrapper',
  'distortion',
  'dggs',
]

# Run tests
suite = unittest.TestSuite()
for name in module_names:
    m = package_name + '.tests.test_' + name
    try:
        suite.addTest(unittest.defaultTestLoader.loadTestsFromName(m))
    except (ImportError, AttributeError):
        print("Couldn't find tests for module %s" % m)
        continue        
unittest.TextTestRunner().run(suite)

# Run doctests
for name in module_names:
    m = importlib.import_module(package_name + '.' + name)
    doctest.testmod(m)

# Doctest the introduction of the manual
doctest.testfile('docs/source/introduction.rst')