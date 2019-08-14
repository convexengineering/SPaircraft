"TASOPT c series airfoil fits"
from __future__ import print_function
from builtins import zip
from builtins import range
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
    cl = []
    for m in M_range:
        for n in thick_range:
            for r in re_range:
                dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                for i in range(len(dataf["CD"])):
                    if dataf["CD"][i] and dataf["CL"][i] != 0:
                        cd.append(dataf["CD"][i])
                        cl.append(dataf["CL"][i])
                        re.append(r)
                        tau.append(float(n)/1000)
                        mach.append(m)

    u1 = np.hstack(re)
    print(u1)
    u2 = np.hstack(tau)
    print(u2)
    u3 = np.hstack(mach)
    u4 = np.hstack(cl)
##    print len(u1)
##    print len(u2)
##    print len(u3)
##    print len(u4)
    w = np.hstack(cd)
##    print len(w)
    u1 = u1.astype(np.float)
    u2 = u2.astype(np.float)
    u3 = u3.astype(np.float)
    w = w.astype(np.float)
    u = [u1, u2, u3, u4]
    x = np.log(u)
    y = np.log(w)
    return x, y

def return_fit(u_1, u_2, u_3, u_4):
    "c series airfoil tau, M, cl, and reynolds fit"
    """returned 4 term SMA fit w/RMS of 0.0517691407091:
        w**1.65625 = 191.003 * (u_1)**-0.219123 * (u_2)**3.95226 * (u_3)**19.2722 * (u_4)**1.15491
            + 0.046236 * (u_1)**-0.390044 * (u_2)**0.78635 * (u_3)**-0.340387 * (u_4)**0.953522
            + 2.72291e-12 * (u_1)**1.18241 * (u_2)**-1.75793 * (u_3)**0.10586 * (u_4)**-1.44159
            + 1.60657 * (u_1)**-0.551398 * (u_2)**1.29389 * (u_3)**3.0428 * (u_4)**1.78078

        returned 3 term SMA fir w/RMS of 0.0524738754518 is
        w**0.238549 = 1.05729 * (u_1)**-0.0897274 * (u_2)**0.165956 * (u_3)**0.0331981 * (u_4)**0.208971
            + 1.11829e-06 * (u_1)**0.701152 * (u_2)**-1.09254 * (u_3)**-0.192159 * (u_4)**-1.08356
            + 42.2841 * (u_1)**-0.0792097 * (u_2)**2.09529 * (u_3)**11.7479 * (u_4)**0.480427

        u1 = Re
        u2 = tau
        u3 = M
        u4 = cl

        fitted ranges are defined by:
            Re = range(10000, 35000, 5000)
            thick = ["100", "110", "120", "130", "140", "145"]
            M = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
            cl = np.linspace(0.35, 0.70, 8)
    """
    w = (1.05729 * (u_1)**-0.0897274 * (u_2)**0.165956 * (u_3)**0.0331981 * (u_4)**0.208971
            + 1.11829e-06 * (u_1)**0.701152 * (u_2)**-1.09254 * (u_3)**-0.192159 * (u_4)**-1.08356
            + 42.2841 * (u_1)**-0.0792097 * (u_2)**2.09529 * (u_3)**11.7479 * (u_4)**0.480427)**(1/0.238549)
    return w

def make_fit(thick_range, re_range, M_range):
    #call the fit setup function
    x, y = fit_setup(thick_range, re_range, M_range)

    cstrt, rms = fit(x, y, 4, 'SMA')
    print("RMS")
    print(rms)

def plot_fits(thick_range, re_range, M_range, cl_range):
    "plot fit compared to data"
##    colors = ["k", "m", "b", "g", "y"]#, "r"]#, "c", "m", "k"]
##    res = np.linspace(re_range[0], re_range[-1], 50)
##    for m in M_range:
##        for n in thick_range:
##            i = 0
##            fig, ax = plt.subplots()
##            for i in range(len(re_range)):
##                cd = []
##                cl = []
##                
##                r = re_range[i]
##                dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
##                if len(dataf["CL"]) != 0:
##                    cd.append(dataf["CD"])
##                    cl.append(dataf["CL"])
##                    ax.plot(cl, cd, "o", mec=colors[i], c="None", mew=1.5)
##                i = i+1
##            ax.set_xlabel("$C_{l}$")
##            ax.set_ylabel("$c_{dp}$")
##            ax.grid()
##            ax.set_title('NC%s Drag Polar for M=%s and Re of %sk' % (n, m, r))
##            fig.savefig("wing__data_fits/tasopt_c_series_data_M%s_Re%s.pdf" % (m, r), bbox_inches="tight")
##
##    for m in M_range:
##        for n in thick_range:
##            i = 0
##            fig, ax = plt.subplots()
##            for i in range(len(re_range)):
##                cd = []
##                cl = []
##                
##                r = re_range[i]
##                dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
##                if len(dataf["CL"]) != 0:
##                    cd.append(dataf["CD"])
##                    cl.append(dataf["CL"])
##                    ax.plot(np.log(cl), np.log(cd), "o", mec=colors[i], c="None", mew=1.5)
##                i = i+1                
##            ax.set_xlabel("Log of $C_{l}$")
##            ax.set_ylabel("Log of $c_{dp}$")
##            ax.grid()
##            ax.set_title('NC%s Drag Polar for M=%s and Re of %sk' % (n, m, r))
##            fig.savefig("wing_data_fits/log_tasopt_c_series_data_M%s_Re%s.pdf" % (m, r), bbox_inches="tight")
##                    
##    colors = ["k", "m", "b", "g", "y", "r", "c", "m", "k"]
##    for r in re_range:
##        fig, ax = plt.subplots()
##        for n, col in zip(thick_range, colors):
##            i = 0
##            for i in range(len(cl_range)):
##                cd = []
##                m_vec = []
##                cl = []
##                for m in M_range:
##                    dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
##                    for j in range(len(dataf["CL"])):
##                        if dataf["CL"][j] <= cl_range[i]+0.01 and dataf["CL"][j] >= cl_range[i]-0.01:
##                            cd.append(dataf["CD"][j])
##                            m_vec.append(m)
##                            cl.append(dataf["CL"][j])
##                ax.plot(m_vec, cd, "o", mec=colors[i], c="None", mew=1.5)
##                i = i+1
##            ax.legend(cl_range, loc=2, fontsize=15)
##            ax.set_xlabel("M")
##            ax.set_ylabel("$c_{dp}$")
##            ax.grid()
##            ax.set_title('NC%s Drag Rise for Re of %sk' % (n, r))
##            fig.savefig("m_fits/tasopt_NC%s_Re%s_drag_rise.pdf" % (n, r), bbox_inches="tight")
    
    colors = ["k", "m", "b", "g", "y", "r", "c", "m", "k"]
    res = np.linspace(re_range[0], re_range[-1], 50)
    for r in re_range:
        for n, col in zip(thick_range, colors):
            i = 0
            fig, ax = plt.subplots()
            for i in range(len(cl_range)):
                cd = []
                m_vec = []
                cl = []
                refit = []
                w = []
                for m in M_range:
                    ms = res = np.linspace(M_range[0], M_range[-1], len(res))
                    dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                    for j in range(len(dataf["CL"])):
                        if dataf["CL"][j] <= cl_range[i]+0.01 and dataf["CL"][j] >= cl_range[i]-0.01:
                            cd.append(dataf["CD"][j])
                            m_vec.append(m)
                            cl.append(dataf["CL"][j])
                            refit.append(r)
                ax.plot(m_vec, cd, "o", mec=colors[i], c="None", mew=1.5)
                for h in range(len(refit)):
                    w.append(return_fit(refit[h], float(n)/1000., m_vec[h], cl[h]))
                    h = h+1
                ax.plot(m_vec, w, c=colors[i], label="NC%s" % n, lw=2)
                i = i+1
##                ax.legend(cl_range, loc=2, fontsize=15)
            ax.set_xlabel("Log of M")
            ax.set_ylabel("Log of $c_{dp}$")
            ax.grid()
            ax.set_title('NC%s Drag Rise for Re of %sk' % (n, r))
            fig.savefig("m_fits/log_tasopt_NC%s_Re%s_drag_rise.pdf" % (n, r), bbox_inches="tight")

    #plot fixed cl diff curves based off of mach number cd vs re
##    colors = ["k", "m", "b", "g", "y"]#, "r"]#, "c", "m", "k"]
##    for m in M_range:
##        fig, ax = plt.subplots()
##        for n, col in zip(thick_range, colors):
##            i = 0
##            
##            for i in range(len(cl_range)):
##                cd = []
##                cl = []
##                re_plot = []
##                for r in re_range:
##                    dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
##                    for j in range(len(dataf["CL"])):
##                        if dataf["CL"][j] <= cl_range[i]+0.01 and dataf["CL"][j] >= cl_range[i]-0.01:
##                            cd.append(dataf["CD"][j])
##                            re_plot.append(r)
##                ax.plot(re_plot, cd, "o", mec='k', c="None", mew=1.5)
##                i = i+1
####                ax.legend(cl_range, loc=2, fontsize=15)
##            ax.set_xlabel("Log of M")
##            ax.set_ylabel("Log of $c_{dp}$")
##            ax.grid()
##            ax.set_title('NC%s Drag Rise for Re of %sk' % (n, r))
##            fig.savefig("re_fits/log_tasopt_NC%s_Re%s_drag_rise.pdf" % (n, r), bbox_inches="tight")

def plot_data(thick_range, re_range, M_range, cl_range):
    "plot x foil data"

    colors = ["k", "m", "b", "g", "y"]#, "r"]#, "c", "m", "k"]
    res = np.linspace(re_range[0], re_range[-1], 50)
    for m in M_range:
        for n in thick_range:
            i = 0
            fig, ax = plt.subplots()
            for i in range(len(re_range)):
                cd = []
                cl = []
                
                r = re_range[i]
                dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CL"]) != 0:
                    cd.append(dataf["CD"])
                    cl.append(dataf["CL"])
                    ax.plot(cl, cd, "o", mec=colors[i], c="None", mew=1.5)
                i = i+1
            ax.set_xlabel("$C_{l}$")
            ax.set_ylabel("$c_{dp}$")
            ax.grid()
            ax.set_title('NC%s Drag Polar for M=%s and Re of %sk' % (n, m, r))
            fig.savefig("wing_data/tasopt_c_series_data_M%s_Re%s.pdf" % (m, r), bbox_inches="tight")

    for m in M_range:
        for n in thick_range:
            i = 0
            fig, ax = plt.subplots()
            for i in range(len(re_range)):
                cd = []
                cl = []
                
                r = re_range[i]
                dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                if len(dataf["CL"]) != 0:
                    cd.append(dataf["CD"])
                    cl.append(dataf["CL"])
                    ax.plot(np.log(cl), np.log(cd), "o", mec=colors[i], c="None", mew=1.5)
                i = i+1                
            ax.set_xlabel("Log of $C_{l}$")
            ax.set_ylabel("Log of $c_{dp}$")
            ax.grid()
            ax.set_title('NC%s Drag Polar for M=%s and Re of %sk' % (n, m, r))
            fig.savefig("wing_data/log_tasopt_c_series_data_M%s_Re%s.pdf" % (m, r), bbox_inches="tight")
                    
    colors = ["k", "m", "b", "g", "y", "r", "c", "m", "k"]
    for r in re_range:
        fig, ax = plt.subplots()
        for n, col in zip(thick_range, colors):
            i = 0
            for i in range(len(cl_range)):
                cd = []
                m_vec = []
                cl = []
                for m in M_range:
                    dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                    for j in range(len(dataf["CL"])):
                        if dataf["CL"][j] <= cl_range[i]+0.01 and dataf["CL"][j] >= cl_range[i]-0.01:
                            cd.append(dataf["CD"][j])
                            m_vec.append(m)
                            cl.append(dataf["CL"][j])
                ax.plot(m_vec, cd, "o", mec=colors[i], c="None", mew=1.5)
                i = i+1
            ax.legend(cl_range, loc=2, fontsize=15)
            ax.set_xlabel("M")
            ax.set_ylabel("$c_{dp}$")
            ax.grid()
            ax.set_title('NC%s Drag Rise for Re of %sk' % (n, r))
            fig.savefig("m_data/tasopt_NC%s_Re%s_drag_rise.pdf" % (n, r), bbox_inches="tight")

        for r in re_range:
            fig, ax = plt.subplots()
            for n, col in zip(thick_range, colors):
                i = 0
                for i in range(len(cl_range)):
                    cd = []
                    m_vec = []
                    cl = []
                    for m in M_range:
                        dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                        for j in range(len(dataf["CL"])):
                            if dataf["CL"][j] <= cl_range[i]+0.01 and dataf["CL"][j] >= cl_range[i]-0.01:
                                cd.append(dataf["CD"][j])
                                m_vec.append(m)
                                cl.append(dataf["CL"][j])
                    ax.plot(np.log(m_vec), np.log(cd), "o", mec=colors[i], c="None", mew=1.5)
                    i = i+1
##                ax.legend(cl_range, loc=2, fontsize=15)
                ax.set_xlabel("Log of M")
                ax.set_ylabel("Log of $c_{dp}$")
                ax.grid()
                ax.set_title('NC%s Drag Rise for Re of %sk' % (n, r))
                fig.savefig("m_data/log_tasopt_NC%s_Re%s_drag_rise.pdf" % (n, r), bbox_inches="tight")

    #plot fixed cl diff curves based off of mach number cd vs re
    colors = ["k", "m", "b", "g", "y"]#, "r"]#, "c", "m", "k"]
    for m in M_range:
        fig, ax = plt.subplots()
        for n, col in zip(thick_range, colors):
            i = 0
            
            for i in range(len(cl_range)):
                cd = []
                cl = []
                re_plot = []
                for r in re_range:
                    dataf = text_to_df("blade.c%s.Re%dk.M%s.pol" % (n, r, m))
                    for j in range(len(dataf["CL"])):
                        if dataf["CL"][j] <= cl_range[i]+0.01 and dataf["CL"][j] >= cl_range[i]-0.01:
                            cd.append(dataf["CD"][j])
                            re_plot.append(r)
                ax.plot(re_plot, cd, "o", mec='k', c="None", mew=1.5)
                i = i+1
##                ax.legend(cl_range, loc=2, fontsize=15)
            ax.set_xlabel("Log of M")
            ax.set_ylabel("Log of $c_{dp}$")
            ax.grid()
            ax.set_title('NC%s Drag Rise for Re of %sk' % (n, r))
            fig.savefig("re_data/log_tasopt_NC%s_Re%s_drag_rise.pdf" % (n, r), bbox_inches="tight")

if __name__ == "__main__":
    Re = list(range(10000, 35000, 5000))
    thick = ["100", "110", "120", "130", "140", "145"]
    M = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    cl = np.linspace(0.35, 0.70, 8)
##    X, Y = fit_setup(thick, Re, M) # call fit(X, Y, 4, "SMA") to get fit
##    make_fit(thick, Re, M)
##    plot_data(thick, Re, M, cl)
    plot_fits(thick, Re, M, cl)
    
