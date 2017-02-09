# RUN IN CONSOLE FOR BEST RESULTS

import csv

Wempty = sol('W_{dry}')
Wpay = sol('W_{payload}')
Wfuel = sol('W_{f_{total}}')
max = sol('W_{f_{total}}')
WMTO = sol('W_{total}')

Whpesys = 0. * sol('W_{dry}')
Wlgnose = 0. * sol('W_{dry}')
Wlgmain = 0. * sol('W_{dry}')
Wtotadd = 0. * sol('W_{dry}')

Wfix = sol('W_{fix}')
Wapu = sol('W_{apu}')
Wpadd = sol('W_{padd}')
Wshell = sol('W_{shell}')
Wcone = sol('W_{cone}')
Whbend = sol('W_{hbend}')
Wvbend = sol('W_{vbend}')
Wwindow = sol('W_{window}')
Winsul = sol('W_{insul}')
Wfloor = sol('W_{floor}')
Wseat = sol('W_{seat}')
Wfuse = sol('W_{fuse}')

Wcap = 0. * sol('W_{dry}')
Wweb = sol('W_{db}')
Wflap = 0. * sol('W_{dry}')
Wslat = 0. * sol('W_{dry}')
Waile = 0. * sol('W_{dry}')
Wlete = 00. * sol('W_{dry}')
Wribs = 0. * sol('W_{dry}')
Wspoi = 0. * sol('W_{dry}')
Wwatt = 0. * sol('W_{dry}')
Wwing = sol('W_{wing}')

Wstrut = 0. * sol('W_{dry}')
Whtail = sol('W_{HT}')
Wvtail = sol('W_{VT}')

Webare = 0. * sol('W_{dry}')
Weadd = 0. * sol('W_{dry}')
Wnace = 0. * sol('W_{dry}')
Wpylon = 0. * sol('W_{dry}')
Weng = sol('W_{engine}').to('lbf')

arr = [Wempty, Wpay, Wfuel, max, WMTO,
       Whpesys, Wlgnose, Wlgmain, Wtotadd,
       Wfix, Wapu, Wpadd, Wshell, Wcone, Whbend, Wvbend, Wwindow, Winsul, Wfloor, Wseat, Wfuse,
       Wcap, Wweb, Wflap, Wslat, Waile, Wlete, Wribs, Wspoi, Wwatt, Wwing,
       Wstrut, Whtail, Wvtail,
       Webare, Weadd, Wnace, Wpylon, Weng]

with open('TASOPT_weight_fractions.csv', 'rb') as csvfile:
    data = list(csv.reader(csvfile))
ct = 0
with open('TASOPT_weight_fractions.csv','wb') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    for row in data:
        ct += 1
        if ct == 2:
            row.append(row[1])
            row.append(row[2])
        if ct >= 3:
            row.append(str(arr[ct - 3].magnitude))
            row.append(str(arr[ct - 3].magnitude/WMTO.magnitude))
        print row
        writer.writerow(row)

