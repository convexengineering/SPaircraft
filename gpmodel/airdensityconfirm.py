from gp_d8 import Mission, Wing, Aircraft,  init_subs

wing = Wing()
ac = Aircraft(wing)
M = Mission(ac)
init_subs(M)
M.cost = M["W_{fuel-tot}"]
sol = M.solve("mosek")

print "air density INPUT: %.3f" % sol("\\rho")[0].magnitude

M.substitutions.update({"S": 3})
for vk in M.varkeys["\\rho"]:
    del M.substitutions[vk]
sol = M.solve("mosek")

print "air density optimum: %.3f" % sol("\\rho")[0].magnitude
for i, sv in enumerate(sol("\\rho")):
    print "cruise %d air density: %.3f" % (i, sv.magnitude)
