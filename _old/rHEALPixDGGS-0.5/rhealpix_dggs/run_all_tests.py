# Run all test scripts and doctests.
import unittest, doctest, importlib

module_names = [
  'projection_tools',
  'pj_healpix',
  'pj_rhealpix',
  'ellipsoids',
  'projection_wrapper',
  'rhealpix_dggs',
  'distortion',
]

# Run tests
suite = unittest.TestSuite()
for name in module_names:
    t = 'tests.test_' + name
    try:
        # If the module defines a suite() function, call it to get the suite.
        mod = __import__(t, globals(), locals(), ['suite'])
        suitefn = getattr(mod, 'suite')
        suite.addTest(suitefn())
    except (ImportError, AttributeError):
        # else, just load all the test cases from the module.
        suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))
unittest.TextTestRunner().run(suite)

# Run doctests
for name in module_names:
    m = importlib.import_module(name)
    doctest.testmod(m)