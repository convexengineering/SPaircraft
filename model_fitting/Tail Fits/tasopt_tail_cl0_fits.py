"tasopt_tail_polarfits.py"
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from gpfit.fit import fit
plt.rcParams.update({'font.size':15})

def text_to_df(filename):
    "parse XFOIL polars and concatente data in DataFrame"
    lines = list(open(filename))
    for i, l in enumerate(lines):
        lines[i] = l.split("\n")[0]
        for j in 10-np.arange(9):
            if " "*j in lines[i]:
                lines[i] = lines[i].replace(" "*j, " ")
            if "---" in lines[i]:
                start = i
    data = {}
    titles = lines[start-1].split(" ")[1:]
    for t in titles:
        data[t] = []

    for l in lines[start+1:]:
        for i, v in enumerate(l.split(" ")[1:]):
            data[titles[i]].append(v)

    df = pd.DataFrame(data)
    df = df.astype(float)
    return df

def fit_setup(thick_range, re_range, M_range):
    "set up x and y parameters for gp fitting"
    cd = []
    tau = []
    mach = []
    re = []
    for m in M_range:
        for n in thick_range:
            for r in re_range:
                dataf = text_to_df("blade.t%s.cl0.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CD"]) != 0:
                    cd.append(dataf["CD"])
                    re.append(r)
                    tau.append(float(n))
                    mach.append(m)

    u1 = np.hstack(re)
    u2 = np.hstack(tau)
    print(u2)
    u3 = np.hstack(mach)
    w = np.hstack(cd)
    u1 = u1.astype(np.float)
    u2 = u2.astype(np.float)
    u3 = u3.astype(np.float)
    w = w.astype(np.float)
    u = [u1, u2, u3]
    x = np.log(u)
    y = np.log(w)
    return x, y

def return_fit(u_1, u_2, u_3):
    "naca tau, M, and reynolds fit"
    """returned fit is:
        w**6.48983 = 5.28751e-20 * (u_1)**0.900672 * (u_2)**0.912222 * (u_3)**8.64547
        + 1.67605e-28 * (u_1)**0.350958 * (u_2)**6.29187 * (u_3)**10.2559
        + 7.09757e-25 * (u_1)**1.39489 * (u_2)**1.96239 * (u_3)**0.567066
        + 3.73076e-14 * (u_1)**-2.57406 * (u_2)**3.12793 * (u_3)**0.448159
        + 1.44343e-12 * (u_1)**-3.91046 * (u_2)**4.66279 * (u_3)**7.68852

        u1 = Re
        u2 = tau
        u3 = M

        only covers 0008-0020
        RMS is 0.01140494297, 5 term SMA
    """
    w = (5.28751e-20 * (u_1)**0.900672 * (u_2)**0.912222 * (u_3)**8.64547
    + 1.67605e-28 * (u_1)**0.350958 * (u_2)**6.29187 * (u_3)**10.2559
    + 7.09757e-25 * (u_1)**1.39489 * (u_2)**1.96239 * (u_3)**0.567066
    + 3.73076e-14 * (u_1)**-2.57406 * (u_2)**3.12793 * (u_3)**0.448159
    + 1.44343e-12 * (u_1)**-3.91046 * (u_2)**4.66279 * (u_3)**7.68852)**(1/6.48983)
    return w

def make_fit(thick_range, re_range, M_range):
    #call the fit setup function
    x, y = fit_setup(thick_range, re_range, M_range)

    cstrt, rms = fit(x, y, 5, 'SMA')
    print("RMS")
    print(rms)

def plot_fits(thick_range, re_range, M_range):
    "plot fit compared to data"

    colors = ["k", "m", "b"]
    assert len(colors) == len(thick_range)
    for m in M_range:
        fig, ax = plt.subplots()
        for n, col in zip(thick_range, colors):
            res = np.linspace(re_range[0], re_range[-1], 50)
            cd = []
            re_array = np.array(re_range)
            i = 0
            delcount = 0
            for i in range(len(re_range)):
                r = re_range[i]
                dataf = text_to_df("blade.t%s.cl0.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CD"]) != 0:
                    cd.append(dataf["CD"])
                else:
                    re_array = np.delete(re_array,i - delcount, 0)
                    delcount = delcount + 1
                i = i+1
            re_range_plot = np.ndarray.tolist(re_array)
            re_range_plot = [i * 1000 for i in re_range_plot]
            ax.plot(re_range_plot, cd, "o", mec=col, c="None", mew=1.5)
            w = return_fit(res, float(n), float(m))
            res = [i * 1000 for i in res]
            ax.plot(res, w, c=col, label="NACA %s" % n, lw=2)
        ax.legend(thick_range, loc=1, fontsize=15)
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.set_ylim([0,0.018])
        ax.set_xlabel("$Re$")
        ax.set_ylabel("$c_{dp}$")
        ax.grid()
        ax.set_title('Profile Drag Coefficient vs Re for M %s' % m)
        fig.savefig("tail_fits/tasopt_taildragpolar_fit_M%s.pdf" % m, bbox_inches="tight")

    
    colors = ["k", "m", "b"]
    res = np.linspace(re_range[0], re_range[-1], 50)
    for m in M_range:
        fig, ax = plt.subplots()
        for n, col in zip(thick_range, colors):
            res = np.linspace(re_range[0], re_range[-1], 50)
            cd = []
            re_array = np.array(re_range)
            i = 0
            delcount = 0
            for i in range(len(re_range)):
                r = re_range[i]
                dataf = text_to_df("blade.t%s.cl0.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CD"]) != 0:
                    cd.append(dataf["CD"])
                else:
                    re_array = np.delete(re_array,i - delcount, 0)
                    delcount = delcount + 1
                i = i+1
            re_range_plot = np.ndarray.tolist(re_array)
            re_range_plot = [i * 1000 for i in re_range_plot]
            ax.plot(np.log(re_range_plot), np.log(cd), "o", mec=col, c="None", mew=1.5)
            w = return_fit(res, float(n), float(m))
            res = [i * 1000 for i in res]
            ax.plot(np.log(res), np.log(w), c=col, label="NACA %s" % n, lw=2)
        ax.legend(thick_range, loc=1, fontsize=15)
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.set_xlabel("Log of $Re$")
        ax.set_ylabel("Log of $c_{dp}$")
        ax.grid()
        ax.set_title('Profile Drag Coefficient vs Re for M %s' % m)
        fig.savefig("tail_fits/tasopt_log_log_taildragpolarfit_M%s.pdf" % m, bbox_inches="tight")

def plot_data(thick_range, re_range, M_range):
    "plot x foil data"

    colors = ["k", "m", "b"]
    assert len(colors) == len(thick_range)
    res = np.linspace(re_range[0], re_range[-1], 50)
    for m in M_range:
        fig, ax = plt.subplots()
        for n, col in zip(thick_range, colors):
            cd = []
            re_array = np.array(re_range)
            i = 0
            delcount = 0
            for i in range(len(re_range)):
                r = re_range[i]
                dataf = text_to_df("blade.t%s.cl0.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CD"]) != 0:
                    cd.append(dataf["CD"])
                else:
                    re_array = np.delete(re_array,i - delcount, 0)
                    delcount = delcount + 1
                i = i+1
            re_range_plot = np.ndarray.tolist(re_array)
            re_range_plot = [i * 1000 for i in re_range_plot]
            ax.plot(re_range_plot, cd, "o", mec=col, c="None", mew=1.5)
        ax.legend(thick_range, loc=1, fontsize=15)
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.set_ylim([0,0.011])
        ax.set_xlabel("$Re$")
        ax.set_ylabel("$c_{dp}$")
        ax.grid()
        ax.set_title('Profile Drag Coefficient vs Re for M %s' % m)
        fig.savefig("tail_fits/tasopt_taildragpolar_data_M%s.pdf" % m, bbox_inches="tight")

    
    colors = ["k", "m", "b"]
    res = np.linspace(re_range[0], re_range[-1], 50)
    for m in M_range:
        fig, ax = plt.subplots()
        for n, col in zip(thick_range, colors):
            cd = []
            re_array = np.array(re_range)
            i = 0
            delcount = 0
            for i in range(len(re_range)):
                r = re_range[i]
                dataf = text_to_df("blade.t%s.cl0.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CD"]) != 0:
                    cd.append(dataf["CD"])
                else:
                    re_array = np.delete(re_array,i - delcount, 0)
                    delcount = delcount + 1
                i = i+1
            re_range_plot = np.ndarray.tolist(re_array)
            re_range_plot = [i * 1000 for i in re_range_plot]
            ax.plot(np.log(re_range_plot), np.log(cd), "o", mec=col, c="None", mew=1.5)
        ax.legend(thick_range, loc=1, fontsize=15)
        ax.set_xlabel("Log of $Re$")
        ax.set_ylabel("Log of $c_{dp}$")
        ax.grid()
        ax.set_title('Log of Profile Drag Coefficient vs log of Re for M %s' % m)
        fig.savefig("tail_fits/tasopt_log_log_taildragpolar_data_M%s.pdf" % m, bbox_inches="tight")


if __name__ == "__main__":
    Re = range(500, 9500, 500)
    thick = ["100", "120", "140"]
    M = [0.4, 0.6, 0.8]
##    X, Y = fit_setup(thick, Re) # call fit(X, Y, 4, "SMA") to get fit
##    make_fit(thick, Re, M)
##    plot_data(thick, Re, M)
    plot_fits(thick, Re, M)
    
