import numpy as np
import matplotlib.pyplot as plt
import sys
# Change this
#sys.path.append("/Users/philippekirschen/Documents/MIT/Research/SPaircraft")
from gpkit import Model, Variable
from SP_aircraft import run_optimal_737
from SP_aircraft import run_optimal_777
from SP_aircraft import run_optimal_D8

def convergence_plots(aircraft):
    if aircraft == "737":
        sol = run_optimal_737('W_{f_{total}}', True, False, False)
    if aircraft == "777":
        sol = run_optimal_777('W_{f_{total}}', True, False, False)    
    if aircraft == "D8":
        sol = run_optimal_D8('W_{f_{total}}', True, False, False)    

    fuel = []
    fuel_only = []
    relax = []
    cost = []
    i = 0
    xaxis_fuel = []
    xaxis_relax = []
    xaxis = []
    print sol.program.gps
    for gp in sol.program.gps:
        cost += [gp.result['cost']]
        hold = []
        i = i+1
        xaxis.append(i)
        varkeys = [k for k in gp.varlocs if "Relax" in k.models and
                                            gp.result(k) >= 1.00001]
        if varkeys:
            for v in varkeys:
                hold.append(gp.result['freevariables'][v])
            relax.append(sum(hold)**20)
            xaxis_relax.append(i)
            fuel += [gp.result['cost']/sum(hold)**20]
        if not varkeys:
            fuel += [gp.result['cost']]
            fuel_only += [gp.result['cost']]
            xaxis_fuel.append(i)
            relax.append(1)
        print varkeys

    print relax
    print fuel
    #plot cost contribution of both fuel burn and relaxed variables
    #on top of one another
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].semilogy(xaxis, relax, '-o')
    axarr[1].semilogy(xaxis, fuel, '-o')
    if aircraft == "737":
        axarr[0].set_title('737 Relaxed Variables Cost Contribution', fontsize=18)
        axarr[1].set_title('737 Fuel Burn Cost Contribution', fontsize=18)
    if aircraft == "777":
        axarr[0].set_title('777 Relaxed Variables Cost Contribution', fontsize=18)
        axarr[1].set_title('777 Fuel Burn Cost Contribution', fontsize=18)   
    if aircraft == "D8":
        axarr[0].set_title('D8 Relaxed Variables Cost Contribution', fontsize=18)
        axarr[1].set_title('D8 Fuel Burn Cost Contribution', fontsize=18)
    axarr[0].set_ylabel('Cost Contribution', fontsize=18)
    axarr[1].set_ylabel('Cost Contribution', fontsize=18)
    f.text(0.5, 0.04, 'GP Iteration', ha='center', va='center', fontsize=18)
    plt.savefig('fuel_and_relax_cost.pdf')
    if aircraft == "737":
        plt.savefig('737_fuel_and_relax_cost.pdf', bbox_inches="tight")
    if aircraft == "777":
        plt.savefig('777_fuel_and_relax_cost.pdf', bbox_inches="tight")  
    if aircraft == "D8":
        plt.savefig('D8_fuel_and_relax_cost.pdf', bbox_inches="tight")
    plt.show()

    #plot the fuel burn cost only when no relaxed variables
    plt.plot(xaxis_fuel, fuel_only, '-o')
    if aircraft == "737":
        plt.title('737 Fuel Burn Cost Contribution', fontsize=18)
    if aircraft == "777":
        plt.title('777 Fuel Burn Cost Contribution', fontsize=18)  
    if aircraft == "D8":
        plt.title('D8 Fuel Burn Cost Contribution', fontsize=18)
    plt.xticks(xaxis_fuel)
    plt.ylabel('Fuel Burn Cost Contribution', fontsize=18)
    plt.xlabel('GP Iteration', fontsize=18)
    if aircraft == "737":
        plt.savefig('737_fuel_cost.pdf', bbox_inches="tight")
    if aircraft == "777":
        plt.savefig('777_fuel_cost.pdf', bbox_inches="tight")
    if aircraft == "D8":
        plt.savefig('D8_fuel_cost.pdf', bbox_inches="tight")
    plt.show()

    #plot the total cost
    plt.semilogy(xaxis, cost, '-o')
    if aircraft == "737":
        plt.title('737 Total Cost', fontsize=18)
    if aircraft == "777":
        plt.title('777 Total Cost', fontsize=18) 
    if aircraft == "D8":
        plt.title('D8 Total Cost', fontsize=18)
    plt.xticks(xaxis)
    plt.ylabel('Total Cost', fontsize=18)
    plt.xlabel('GP Iteration', fontsize=18)
    if aircraft == "737":
        plt.savefig('737_total_cost.pdf', bbox_inches="tight")
    if aircraft == "777":
        plt.savefig('777_total_cost.pdf', bbox_inches="tight")
    if aircraft == "D8":
        plt.savefig('D8_total_cost.pdf', bbox_inches="tight")
    plt.show()
