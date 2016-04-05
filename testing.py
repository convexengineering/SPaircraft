"""Script for running all gpkit-models test() methods"""
import unittest
from gpkit.tests.helpers import run_tests
MODULES = []

import beam
MODULES.append(beam)

import troposphere
MODULES.append(troposphere)
import troposphere_sp
MODULES.append(troposphere_sp)

from aircraft import breguet_range, cost
MODULES.append(breguet_range)
MODULES.append(cost)


def test_generator(model_class):
    """Returns method that tests model_class, to be attached to a TestCase"""
    def test(self):
        m = model_class()
        m.test()
    return test


class TestModels(unittest.TestCase):
    """Single TestCase for all Models in gpkit-models"""
    pass


def attach_tests():
    """Gather Models that have test() methods; attach to TestModels"""
    for module in MODULES:
        for attr_name in dir(module):
            candidate = getattr(module, attr_name)
            try:
                if candidate.__module__ != module.__name__:
                    # skip things imported by this module, like np
                    continue
                if hasattr(candidate, 'test'):
                    setattr(TestModels,
                            "test_%s" % candidate.__name__,
                            test_generator(candidate))
            except AttributeError:
                # dict attributes like __builtins__ don't have __module__
                pass


def run(xmloutput=False):
    """Run all gpkit-models unit tests.

    Arguments
    ---------
    xmloutput: bool
        If true, generate xml output files for continuous integration
    """
    if xmloutput:
        raise NotImplementedError("gpkit-models CI has not yet been set up")
    attach_tests()
    run_tests([TestModels])


if __name__ == '__main__':
    run()
