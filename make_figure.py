from __future__ import absolute_import
from math import sin, cos, pi
from itertools import combinations

"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS

# only needed for the local bounded debugging tool
from collections import defaultdict
from gpkit.small_scripts import mag

# importing from D8_integration
from .stand_alone_simple_profile import FlightState, Altitude, Atmosphere
from D8_VT_yaw_rate_and_EO_simple_profile import VerticalTail, VerticalTailPerformance
from Wing_simple_performance import Wing, WingPerformance
from D8_integration import Engine, EnginePerformance
from CFP_Fuselage_Performance_int_HT import Mission,Fuselage, FuselagePerformance, Aircraft, AircraftP, HorizontalTail, HorizontalTailPerformance


def make_figure():
    m = Mission()
    # sol = m.localsolve(verbosity=2)

    ac = Aircraft()
    acP = AircraftP(ac, FlightState())
    subsystems = [Engine(),Wing(),HorizontalTail(),VerticalTail(),Fuselage()]

    # keychain = {}
    # for vk in ac.varkeys:
    #     keychain[vk.str_without(['models'])] = vk.models
    #     model = vk.descr.get("models")
    #     if 'Aircraft' in model or 'AircraftP' in model:  # and vk not in ac.substitutions:
    #         v = vk.str_without(['models'])
    #         keychain[v] = []
    #         for ss in subsystems:
    #             for i in ss.varkeys:
    #                 if v == i.str_without(['models']):
    #                     keychain[v].append(type(ss).__name__)
    #         if len(keychain[v]) <= 1:
    #             del keychain[v]
    # for vk in acP.varkeys:
    #     keychain[vk.str_without(['models'])] = vk.models
    # for vk in m.varkeys:
    #      keychain[vk.str_without(['models'])] = vk.models
    #
    # for key in keychain.keys():
    #     if 'VerticalTailNoStruct' in keychain[key] :
    #         keychain[key].remove('VerticalTailNoStruct')
    #     if 'HorizontalTailNoStruct' in keychain[key]:
    #         keychain[key].remove('HorizontalTailNoStruct')
    #     if 'WingNoStruct' in keychain[key]:
    #         keychain[key].remove('WingNoStruct')
    #     if 'WingBox' in keychain[key]:
    #         keychain[key].remove('WingBox')

    keychain = {}

    constraints = ac.flat(constraintsets=False)
    constraintsP = acP.flat(constraintsets=False)
    counter = 0


    linkedModels = []
    linkedVariables = []
    for i in constraints:
        counter += 1
        # Getting the variables in each constraint
        vars = i.varkeys#
        # print vars
        for i in vars:
            linkedModels[counter].append(i.descr.get('models'))
            linkedVariables[counter].append(i.str_without('models'))
        # varList = []
        # for i in vars:
        #     model.add(i.str)
        # keychain[counter] = i
    #
    # for i in constraintsP:
    #     counter += 1
    #     vars = i.varkeys
    #     keychain[counter] = i


    # # counter = 0
    # # modellinks = {}
    # # while counter < len(keychain):
    # #     counter += 1
    # #     modellinks[counter]
    #
    #
    # # Get all combinations of sub-modelsa
    # modellist = [type(ss).__name__ for ss in subsystems]
    # modellistcombo = []
    # for n in range(2, len(modellist)):
    #     modellistcombo += combinations(modellist, n)
    #
    # # Create a dictionary that has each combination (tuple) as a key
    # zidane = {}
    # for modelgroup in modellistcombo:
    #     zidane[modelgroup] = []
    #
    # # Fill this dictionary in with all the varkeys for which each combo applies
    # for key in keychain:
    #     zidane[tuple(keychain[key])] += [key]
    #
    # # Get rid of any empty entries
    # viera = {k: v for k, v in zidane.items() if v}
    #
    # with open('figs/linked_vars.tex', 'w') as outfile:
    #
    #     outfile.write('')
    #
    #     outfile.write('\\begin{center}\n' +
    #                   '\\begin{tikzpicture}[thick]\n')
    #
    #     # Create a circular ring of nodes, one for each model
    #     nodepos = {}
    #     i = 0
    #     I = len(modellist)
    #     R = 6
    #     for model in modellist:
    #         nodepos[model] = (R * sin(2 * pi * i / I), R * cos(2 * pi * i / I))
    #         outfile.write('\\node at {0}'.format(str(nodepos[model])) +
    #                       '[circle, minimum size=3cm, fill=blue!20]' +
    #                       '({0}){{{0}}};\n'.format(model))
    #         i = i + 1
    #
    #     j = 0
    #     colours = ['red', 'blue', 'cyan ', 'magenta', 'brown', 'gray',
    #                'olive', 'orange']
    #     for key in viera:
    #         # Create a node for every group of variables
    #         vargroup = viera[key]
    #         varnodeposx = 0
    #         varnodeposy = 0
    #         for k in key:
    #             varnodeposx += nodepos[k][0]
    #             varnodeposy += nodepos[k][1]
    #         varnodeposx = varnodeposx / len(key)
    #         varnodeposy = varnodeposy / len(key)
    #         outfile.write('\\node at (%.2f,%.2f)' % (varnodeposx, varnodeposy) +
    #                       '[draw,rectangle, color=%s, align=left, fill=white]({'
    #                       % colours[j] +
    #                       ' '.join(vargroup).replace('\\', '') + '}){$' +
    #                       '$\\\\$'.join(vargroup) + '$};\n')
    #         # Create edges from vargroups to models they appear in
    #         for k in key:
    #             outfile.write('\\begin{pgfonlayer}{background}\n' +
    #                           '\\draw[color=%s] ({' % colours[j] +
    #                           ' '.join(vargroup).replace('\\', '') +
    #                           '})--(' + k + ');\n' +
    #                           '\\end{pgfonlayer}\n')
    #
    #         j += 1
    #
    #     outfile.write('\\end{tikzpicture}\n' +
    #                   '\\end{center}')

if __name__ == "__main__":
    make_figure()
