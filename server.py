from SimpleWebSocketServer import SimpleWebSocketServer, WebSocket
import json
from SPaircraft import Mission
from saveSol import gendes, gencsm
from shutil import copyfile

from gpkit import units

from subs.optimalD8 import get_optimalD8_subs
from aircraft import Mission
from SPaircraft import optimize_aircraft

EXIT = [False]
ID = 0
LASTSOL = None

def genfiles(m, sol):
    global ID
    gensoltxt(m, sol, ID)
    gencsm(m, sol, 'optimalD8', ID)
    copyfile("d82.csm", "d82_%03i.csm" % ID)
    ID += 1


def gensoltxt(m, sol, ID):
    with open("sol_%03i.txt" % ID, "w") as f:
        for var in sorted(m.varkeys, key=str):
            f.write("%s [%s]\t\t%f\n" % (var, var.unitstr(options=":~",
                                                          dimless="-"),
                                         sol["variables"][var]))


class SPaircraftServer(WebSocket):

    def handleMessage(self):
        print "< received", repr(self.data)
        try:
            self.data = json.loads(self.data)
            print self.data

            config = 'optimalD8'
            substitutions = get_optimalD8_subs()
            fixedBPR = False
            pRatOpt = True
            mutategparg = False
            m = Mission(3, 2, config, 1)
            m.cost = m['W_{f_{total}}']

            for name, value in self.data.items():
                try:
                    key = m.design_parameters[name]
                    m.substitutions[key] = value
                except KeyError as e:
                    print repr(e)
            substitutions = get_optimalD8_subs()
            substitutions.update({'R_{req}': 3000.*units('nmi'),
                         'n_{pass}': 180.})
            sol = optimize_aircraft(m, substitutions, fixedBPR, pRatOpt, mutategparg)
            LASTSOL = sol
            genfiles(m, sol)

            self.send({"status": "optimal",
                       "msg": ("Successfully optimized."
                               " Optimal objective: %.1f "
                               % sol(m.cost))})
        except Exception as e:
            self.send({"status": "unknown", "msg": "The last solution"
                      " raised an exception; tweak it and send again."})
            print type(e), e

    def send(self, msg):
        print "> sent", repr(msg)
        self.sendMessage(unicode(json.dumps(msg)))

    def handleConnected(self):
        print self.address, "connected"

    def handleClose(self):
        print self.address, "closed"
        EXIT[0] = True


if __name__ == "__main__":
    Nclimb = 3 # number of climb segments
    Ncruise = 2 # number of cruise segments
    Nmission = 1 # number of missions
    config = 'optimalD8' # String describing configuration:
    # currently one of: 'D8_eng_wing', 'optimal737', 'optimal777', 'optimalD8', 'D8_no_BLI', 'M072_737'
    m = Mission(Nclimb, Ncruise, config, Nmission)
    # Objective
    m.cost = m['W_{f_{total}}'].sum()
    # Inputs to the model
    substitutions = get_optimalD8_subs()
    substitutions.update({'R_{req}': 3000.*units('nmi'), #6000*units('nmi'),
                         'n_{pass}': 180.})              #450.,)

    # Additional options
    fixedBPR = False
    pRatOpt = True
    mutategparg = False
    sol = optimize_aircraft(m, substitutions, fixedBPR, pRatOpt, mutategparg)
    LASTSOL = sol
    genfiles(m, sol)
    server = SimpleWebSocketServer('', 8000, SPaircraftServer)
    while not EXIT[0]:
        server.serveonce()
    print "Python server has exited."
