"naca_polarfits.py"
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

def fit_setup(naca_range, re_range):
    "set up x and y parameters for gp fitting"
    tau = [[float(n)]*len(re_range) for n in naca_range]
    re = [re_range]*len(naca_range)
    cd = []
    for n in naca_range:
        for r in re_range:
            dataf = text_to_df("naca%s.cl0.Re%dk.pol" % (n, r))
            cd.append(dataf["CD"])


    u1 = np.hstack(re)
    u2 = np.hstack(tau)
    w = np.hstack(cd)
    u1 = u1.astype(np.float)
    u2 = u2.astype(np.float)
    w = w.astype(np.float)
    u = [u1, u2]
    x = np.log(u)
    y = np.log(w)
    return x, y

def return_fit(u_1, u_2):
    "naca tau and reynolds fit"
    w = (7.42688e-90 * (u_1)**-33.0637 * (u_2)**18.0419
         + 5.02826e-163 * (u_1)**-18.7959 * (u_2)**53.1879
         + 4.22901e-77 * (u_1)**-41.1704 * (u_2)**28.4609)**(1/70.5599)
    # SMA function, K=3, max RMS error = 0.0173
    return w

def make_fit(naca_range, re_range):
    #call the fit setup function
    x, y = fit_setup(naca_range, re_range)
##    print np.size(x)
##    print np.size(x)
##    print np.size(y)
    fit(x, y, 3, 'SMA')

def plot_fits(naca_range, re_range):
    "plot fit compared to data"

    fig, ax = plt.subplots()
    colors = ["k", "m", "b", "g", "y", "r"]
    assert len(colors) == len(naca_range)
    res = np.linspace(re_range[0], re_range[-1], 50)
    for n, col in zip(naca_range, colors):
        cd = []
        re_array = np.array(re_range)
        i = 0
        delcount = 0
        for i in range(len(re_range)):
            r = re_range[i]
            dataf = text_to_df("naca%s.cl0.Re%dk.pol" % (n, r))
            if len(dataf["CD"]) != 0:
                cd.append(dataf["CD"])
            else:
                re_array = np.delete(re_array,i - delcount, 0)
                delcount = delcount + 1
            i = i+1
        re_range_plot = np.ndarray.tolist(re_array)
        ax.plot(re_range_plot, cd, "o", mec=col, c="None", mew=1.5)
##        w = return_fit(res, float(n))
##        ax.plot(res, w, c=col, label="NACA %s" % n, lw=2)
##    ax.legend(naca_range, loc=4, fontsize=15)
##    labels = ["k" + item.get_text() for item in ax.get_xticklabels()]
##    labels = ["%dk" % l for l in np.linspace(3500, 9000, len(labels))]
##    ax.set_xticklabels(labels)
    ax.set_xlabel("$Re$")
    ax.set_ylabel("$c_{dp}$")
    ax.grid()
##    ax.set_title('Log of Profile Drag Coefficient vs log of Re')
    fig.savefig("tail_fits/taildragpolar_fits.pdf", bbox_inches="tight")

    fig, ax = plt.subplots()
    colors = ["k", "m", "b", "g", "y", "r"]
    res = np.linspace(re_range[0], re_range[-1], 50)
    for n, col in zip(naca_range, colors):
        cd = []
        re_array = np.array(re_range)
        i = 0
        delcount = 0
        for i in range(len(re_range)):
            r = re_range[i]
            dataf = text_to_df("naca%s.cl0.Re%dk.pol" % (n, r))
            if len(dataf["CD"]) != 0:
                cd.append(dataf["CD"])
            else:
                re_array = np.delete(re_array,i - delcount, 0)
                delcount = delcount + 1
            i = i+1
        re_range_plot = np.ndarray.tolist(re_array)
        ax.plot(np.log(re_range_plot), np.log(cd), "o", mec=col, c="None", mew=1.5)
##        w = return_fit(res, float(n))
##        ax.plot(res, w, c=col, label="NACA %s" % n, lw=2)
##    ax.legend(naca_range, loc=4, fontsize=15)
##    labels = ["k" + item.get_text() for item in ax.get_xticklabels()]
##    labels = ["%dk" % l for l in np.linspace(3500, 9000, len(labels))]
##    ax.set_xticklabels(labels)
    ax.set_xlabel("$Re$")
    ax.set_ylabel("$c_{dp}$")
    ax.grid()
    ax.set_title('Log of Profile Drag Coefficient vs log of Re')
    fig.savefig("tail_fits/log_log_taildragpolar_fits.pdf", bbox_inches="tight")

def plot_data(naca_range, re_range):
    "plot fit compared to data"

    fig, ax = plt.subplots()
    colors = ["k", "m", "b", "g", "y", "r"]
    assert len(colors) == len(naca_range)
    res = np.linspace(re_range[0], re_range[-1], 50)
    for n, col in zip(naca_range, colors):
        cd = []
        re_array = np.array(re_range)
        i = 0
        delcount = 0
        for i in range(len(re_range)):
            r = re_range[i]
            dataf = text_to_df("naca%s.cl0.Re%dk.pol" % (n, r))
            if len(dataf["CD"]) != 0:
                cd.append(dataf["CD"])
            else:
                re_array = np.delete(re_array,i - delcount, 0)
                delcount = delcount + 1
            i = i+1
        re_range_plot = np.ndarray.tolist(re_array)
        ax.plot(re_range_plot, cd, "o", mec=col, c="None", mew=1.5)
    ax.set_xlabel("$Re$")
    ax.set_ylabel("$c_{dp}$")
    ax.grid()
    ax.set_title('Profile Drag Coefficient vs Re')
    fig.savefig("tail_fits/taildragpolar_data.pdf", bbox_inches="tight")

    fig, ax = plt.subplots()
    colors = ["k", "m", "b", "g", "y", "r"]
    res = np.linspace(re_range[0], re_range[-1], 50)
    for n, col in zip(naca_range, colors):
        cd = []
        re_array = np.array(re_range)
        i = 0
        delcount = 0
        for i in range(len(re_range)):
            r = re_range[i]
            dataf = text_to_df("naca%s.cl0.Re%dk.pol" % (n, r))
            if len(dataf["CD"]) != 0:
                cd.append(dataf["CD"])
            else:
                re_array = np.delete(re_array,i - delcount, 0)
                delcount = delcount + 1
            i = i+1
        re_range_plot = np.ndarray.tolist(re_array)
        ax.plot(np.log(re_range_plot), np.log(cd), "o", mec=col, c="None", mew=1.5)
    ax.set_xlabel("$Re$")
    ax.set_ylabel("$c_{dp}$")
    ax.grid()
    ax.set_title('Log of Profile Drag Coefficient vs log of Re')
    fig.savefig("tail_fits/log_log_taildragpolar_data.pdf", bbox_inches="tight")


if __name__ == "__main__":
    Re = range(500, 9500, 500)
    NACA = ["0005", "0008", "0009", "0010", "0015", "0020"]
##    X, Y = fit_setup(NACA, Re) # call fit(X, Y, 4, "SMA") to get fit
##    make_fit(NACA, Re)
    plot_data(NACA, Re)
    
