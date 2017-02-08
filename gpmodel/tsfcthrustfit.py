"tasopt output fits"
import numpy as np
from atm import atm
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':15})

def text_to_df(filename):
    "parse TASOPT out file to dataFrame"
    lines = list(open(filename))
    ind2 = 1e6
    ind1 = 1e6
    for i, l in enumerate(lines):
        lines[i] = l.split("\n")[0].replace("\r", "")
        if "Mission profile summary" in l:
            ind1 = i
        if i > ind1 and i <= ind2 and "-------" in l:
            ind2 = i
            break

    lines = lines[ind1+2:ind2]
    del lines[1]
    for i, l in enumerate(lines):
        for j in 10-np.arange(9):
            if " "*j in lines[i]:
                lines[i] = lines[i].replace(" "*j, " ")

    data = {}
    titles = ["FS"] + lines[0].split(" ")[1:]
    for t in titles:
        data[t] = []

    for l in lines[1:]:
        for i, v in enumerate(l.split(" ")[1:]):
            if v[0].isdigit():
                v = float(v)
            data[titles[i]].append(v)

    df = pd.DataFrame(data)
    return df

if __name__ == "__main__":
    df = text_to_df("737.out")
    fig, ax = plt.subplots()
    rho, Tatm, _ = atm(df["h"]*0.3048)
    V = df["Mach"]*(1.4*287*Tatm)**0.5
    drag = 0.5*df["CL"]/df["L/D'"]*rho*V**2*125.67
    ax.plot(drag[df["h"]==0.0], df["TSFC'"][df["h"]==0.0], "o")
    ax.plot(drag[(df["h"]>0.0) & (df["h"]<30000)], df["TSFC'"][(df["h"]>0.0) & (df["h"]<30000)], "o")
    ax.plot(drag[df["h"]>30000], df["TSFC'"][df["h"]>30000], "o")
    fig.savefig("thusttotsfc.pdf")
