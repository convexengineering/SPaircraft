"""Script for running all gpkit-models test() methods"""
import unittest
from gpkit.tests.helpers import run_tests, StdoutCaptured
MODELS = []

import beam
MODELS.append(beam.Beam)

import atmosphere
MODELS.append(atmosphere.Troposphere)
MODELS.append(atmosphere.Tropopause)
MODELS.append(atmosphere.Sutherland)

# from aircraft import breguet_range, cost
# MODELS.append(breguet_range)
# MODELS.append(cost)


class TestModels(unittest.TestCase):
    """Single TestCase for all Models in gpkit-models"""
    pass


def generate_test(modelclass):
    def test(self):
        with StdoutCaptured():
            modelclass.test()
    return test


def attach_tests():
    """Gather Models that have test() methods; attach to TestModels"""
    for modelclass in MODELS:
        setattr(TestModels, "test_%s" % modelclass.__name__,
                generate_test(modelclass))


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
