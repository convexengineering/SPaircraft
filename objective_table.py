"""
script to generate the values in the SP Aircraft different objective table
"""

from SP_aircraft import run_optimal_737

#solve all the cases
base_sol = run_optimal_737('W_{f_{total}}', False, True)
emptyW_sol = run_optimal_737('W_{dry}', False, True)
span_sol = run_optimal_737('b', False, True)
AR_sol = run_optimal_737('AR', False, True)
engineW_sol = run_optimal_737('W_{engine}', False, True)
time_sol = run_optimal_737('Total_Time', False, True)
LD_sol = run_optimal_737('L/D', False, True)
LGw_sol = run_optimal_737('W_{lg}', False, True)

#compute times
base_time = sum(base_sol('thr')['thr_Mission/ClimbSegment/ClimbP/AircraftP'])+sum(base_sol('thr')['thr_Mission/CruiseSegment/CruiseP/AircraftP.1'])
emptyW_time = sum(emptyW_sol('thr')['thr_Mission.1/ClimbSegment.1/ClimbP.1/AircraftP.2'])+sum(emptyW_sol('thr')['thr_Mission.1/CruiseSegment.1/CruiseP.1/AircraftP.3'])
span_time = sum(span_sol('thr')['thr_Mission.2/ClimbSegment.2/ClimbP.2/AircraftP.4'])+sum(span_sol('thr')['thr_Mission.2/CruiseSegment.2/CruiseP.2/AircraftP.5'])
AR_time = sum(AR_sol('thr')['thr_Mission.3/ClimbSegment.3/ClimbP.3/AircraftP.6'])+sum(AR_sol('thr')['thr_Mission.3/CruiseSegment.3/CruiseP.3/AircraftP.7'])
engineW_time = sum(engineW_sol('thr')['thr_Mission.4/ClimbSegment.4/ClimbP.4/AircraftP.8'])+sum(engineW_sol('thr')['thr_Mission.4/CruiseSegment.4/CruiseP.4/AircraftP.9'])
LD_time = sum(LD_sol('thr')['thr_Mission.6/ClimbSegment.6/ClimbP.6/AircraftP.12'])+sum(LD_sol('thr')['thr_Mission.6/CruiseSegment.6/CruiseP.6/AircraftP.13'])
LGw_time = sum(LGw_sol('thr')['thr_Mission.7/ClimbSegment.7/ClimbP.7/AircraftP.14'])+sum(LGw_sol('thr')['thr_Mission.7/CruiseSegment.7/CruiseP.7/AircraftP.15'])

#output the columns of the table
print "column 2"
print "\n"
print [emptyW_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), emptyW_sol('W_{dry}')/base_sol('W_{dry}'), emptyW_sol('b')/base_sol('b'), emptyW_sol('AR')/base_sol('AR'), emptyW_sol('W_{engine}')/base_sol('W_{engine}'), emptyW_time/base_time, emptyW_sol('L/D')['L/D_Mission.1/CruiseSegment.1/CruiseP.1/AircraftP.3'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0],  emptyW_sol('W_{lg}')/base_sol('W_{lg}')]
print "\n"
print "\n"
print "column 3"
print "\n"
print [span_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), span_sol('W_{dry}')/base_sol('W_{dry}'), span_sol('b')/base_sol('b'), span_sol('AR')/base_sol('AR'), span_sol('W_{engine}')/base_sol('W_{engine}'), span_time/base_time, span_sol('L/D')['L/D_Mission.2/CruiseSegment.2/CruiseP.2/AircraftP.5'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0], span_sol('W_{lg}')/base_sol('W_{lg}')]
print "\n"
print "\n"
print "column 4"
print [AR_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), AR_sol('W_{dry}')/base_sol('W_{dry}'), AR_sol('b')/base_sol('b'), AR_sol('AR')/base_sol('AR'), AR_sol('W_{engine}')/base_sol('W_{engine}'), AR_time/base_time, AR_sol('L/D')['L/D_Mission.3/CruiseSegment.3/CruiseP.3/AircraftP.7'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0], AR_sol('W_{lg}')/base_sol('W_{lg}')]
print "\n"
print "\n"
print "column 5"
print [engineW_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), engineW_sol('W_{dry}')/base_sol('W_{dry}'), engineW_sol('b')/base_sol('b'), engineW_sol('AR')/base_sol('AR'), engineW_sol('W_{engine}')/base_sol('W_{engine}'), engineW_time/base_time, engineW_sol('L/D')['L/D_Mission.4/CruiseSegment.4/CruiseP.4/AircraftP.9'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0], engineW_sol('W_{lg}')/base_sol('W_{lg}')]
print "\n"
print "\n"
print "column 6"
print [time_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), time_sol('W_{dry}')/base_sol('W_{dry}'), time_sol('b')/base_sol('b'), time_sol('AR')/base_sol('AR'), time_sol('W_{engine}')/base_sol('W_{engine}'), time_sol('TotalTime')/base_time, time_sol('L/D')['L/D_Mission.5/CruiseSegment.5/CruiseP.5/AircraftP.11'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0], time_sol('W_{lg}')/base_sol('W_{lg}')]
print "\n"
print "\n"
print "column 7"
print [LD_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), LD_sol('W_{dry}')/base_sol('W_{dry}'), LD_sol('b')/base_sol('b'), LD_sol('AR')/base_sol('AR'), LD_sol('W_{engine}')/base_sol('W_{engine}'), LD_time/base_time, LD_sol('L/D')['L/D_Mission.6/CruiseSegment.6/CruiseP.6/AircraftP.13'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0], LD_sol('W_{lg}')/base_sol('W_{lg}')]
print "\n"
print "\n"
print "column 8"
print [LGw_sol('W_{f_{total}}')/base_sol('W_{f_{total}}'), LGw_sol('W_{dry}')/base_sol('W_{dry}'), LGw_sol('b')/base_sol('b'), LGw_sol('AR')/base_sol('AR'), LGw_sol('W_{engine}')/base_sol('W_{engine}'), LGw_time/base_time, LGw_sol('L/D')['L/D_Mission.7/CruiseSegment.7/CruiseP.7/AircraftP.15'][0]/base_sol('L/D')['L/D_Mission/CruiseSegment/CruiseP/AircraftP.1'][0], LGw_sol('W_{lg}')/base_sol('W_{lg}')]

