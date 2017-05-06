import numpy as np
import pandas as pd
from gasmaleperf import Mission
from gpkit.small_scripts import unitstr
from gpkit import Variable
from gen_tex import find_submodels
import xlsxwriter

def mission_vars(M, sol):
    """
    This ouputs variables relevant accross a mission
    """
    mission = ["Climb", "Cruise", "Loiter"]
    mns = [10, 1, 5, 1]
    data = {}

    vks = []
    for m in M.varkeys:
        for fs in mission:
            if fs in m.models and "shape" in m.descr:
                data[m.name + "_" + ", ".join([mname for mname in m.models if mname != fs])] = [unitstr(m.descr["units"])] + [""]*17 \
                               + [m.descr["label"]]
                vks.append(m)

    for vk in vks:
        fs = [m for m in vk.models if m in mission][0]
        mn = vk.modelnums[vk.models.index(fs)]
        if mn == 1:
            ind = -2
        else:
            if fs == "Climb":
                ind = 1 + vk.idx[0]
            elif fs == "Cruise":
                ind = 11
            elif fs == "Loiter":
                ind = 12 + vk.idx[0]
        data[vk.name + "_" + ", ".join([mname for mname in vk.models if mname != fs])][ind] = sol(vk).magnitude


    colnames = ["Units"] + ["Climb%d" % i for i in range(10)] + ["Cruise"] + ["Loiter%d" % i for i in range(5)] + ["Cruise"] + ["Label"]

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    return df

def bd_vars(M, sol, varname, morevars):

    colnames = ["Value", "Units", "Margin", "Margin Sens", "Label"]

    data = {}
    for sv in sol(varname):
        name = max(list(sv.keys), key=len)
        data[name] = [sol(sv).magnitude]
        data[name].append(unitstr(sv.units))
        for mfac in sol("m_{fac}"):
            if not sv.models == mfac.models:
                continue
            data[name].append(sol(mfac).magnitude)
            data[name].append(sol["sensitivities"]["constants"][mfac])
        if len(data[name]) != 4:
            data[name] += [""]*2
        data[name].append(sv.label)

    for name in morevars:
        sv = sol(name)
        data[name] = [sv.magnitude]
        data[name].append(unitstr(M[name].descr["units"]))
        for mfac in sol("m_{fac}"):
            if not M[name].descr["models"] == mfac.models:
                continue
            data[name].append(sol(mfac).magnitude)
            data[name].append(sol["sensitivities"]["constants"][mfac])
        if len(data[name]) != 4:
            data[name] += [""]*2
        data[name].append(M[name].descr["label"])

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    return df

def sketch_params(M, sol, varnames, othervars=None, pointmasses=None):

    data = {}
    for vname in varnames:
        data[vname] = [sol(vname).magnitude, unitstr(M[vname].descr["units"]),
                       M[vname].descr["label"]]

    if othervars:
        data.update(othervars)

    if hasattr(M, "get_cgs"):
        xnp, xcg, SM = M.get_cgs()
        data["x_{np}"] = [xnp.magnitude, xnp.units, "neutral point"]
        data["x_{cg}"] = [xcg.magnitude, xcg.units, "center of gravity"]
        data["SM"] = [SM.magnitude, "-", "static margin"]

    if pointmasses:
        for pm in pointmasses:
            data[pm] = []

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = ["Value", "Units", "Label"]
    return df

def write_to_excel(path, filename, df, sens_formatting):

    coldepth = df.count()[0]
    rowdepth = len(df.columns)

    colind = []
    for colname in df.columns:
        if "Sens" in colname:
            colind.append(list(df.columns).index(colname))

    rowind = []
    for rowname in df.index:
        if "Sens" in rowname:
            rowind.append(list(df.index).index(rowname))

    alp = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]

    writer = pd.ExcelWriter("%s%s" % (path, filename), engine="xlsxwriter")
    df.to_excel(writer, sheet_name="Sheet1")

    workbook = writer.book
    ws = writer.sheets["Sheet1"]
    # Light red fill
    format1 = workbook.add_format({'bg_color': '#FFC7CE'})

    # Green fill
    format2 = workbook.add_format({'bg_color': '#FFCC99'})

    # Green fill
    format3 = workbook.add_format({'bg_color': '#C6EFCE'})

    for i in colind:
        ws.conditional_format("%s2:%s%d" % (alp[i+1], alp[i+1], coldepth+1),
                              {"type": "cell",
                               "criteria": ">=",
                               "value": sens_formatting["bad"],
                               "format": format1})

        ws.conditional_format("%s2:%s%d" % (alp[i+1], alp[i+1], coldepth+1),
                              {"type": "cell",
                               "criteria": "between",
                               "minimum": sens_formatting["bad"],
                               "maximum": sens_formatting["good"],
                               "format": format2})

        ws.conditional_format("%s2:%s%d" % (alp[i+1], alp[i+1], coldepth+1),
                              {"type": "cell",
                               "criteria": "<",
                               "value": sens_formatting["good"],
                               "format": format3})

    for i in rowind:
        ws.conditional_format("%s%d:%s%d" % (alp[2], i+2, alp[rowdepth-1], i+2),
                              {"type": "cell",
                               "criteria": ">=",
                               "value": sens_formatting["bad"],
                               "format": format1})

        ws.conditional_format("%s%d:%s%d" % (alp[2], i+2, alp[rowdepth-1], i+2),
                              {"type": "cell",
                               "criteria": "between",
                               "minimum": sens_formatting["bad"],
                               "maximum": sens_formatting["good"],
                               "format": format2})

        ws.conditional_format("%s%d:%s%d" % (alp[2], i+2, alp[rowdepth-1], i+2),
                              {"type": "cell",
                               "criteria": "<",
                               "value": sens_formatting["good"],
                               "format": format3})

    writer.save()

def model_params(subM, sol):

    data = {}
    for v in subM.varkeys:
        if "idx" not in v.descr or v.idx == (0,):
            if "Cruise" not in v.descr["models"]:
                if "TailAero" not in v.descr["models"]:
                    data[v] = [sol(v.name + "_" + ", ".join(v.models)).magnitude]
                    data[v].append(unitstr(M[v].units))
                    data[v].append(v.descr["label"])

    if data:
        df = pd.DataFrame(data)
        df = df.transpose()
        df.columns = ["Value", "Units", "Label"]
    else:
        df = None
    return df

if __name__ == "__main__":
    M = Mission(DF70=True)
    M.cost = 1/M["t_Mission, Loiter"]
    subs = {"b_Mission, Aircraft, Wing": 24,
            "l_Mission, Aircraft, Empennage, TailBoom": 7.0,
            "AR_v": 1.5, "AR": 24, "SM_{corr}": 0.5, "AR_h": 4, "k": 0.0,
            "(1-k/2)": 1, "d_0": 1}
    M.substitutions.update(subs)
    for p in M.varkeys["P_{avn}"]:
        M.substitutions.update({p: 65})
    for t in M.varkeys["\\theta_{max}"]:
        M.substitutions.update({t: 0.2})
    # JHO.debug(solver="mosek")
    sol = M.solve("mosek")
    PATH = "/Users/mjburton11/Dropbox (MIT)/16.82GasMALE/Management/GpkitReports/"

    # Mission_vars = ["RPM", "BSFC", "V", "P_{shaft}",
    #                 "P_{shaft-tot}", "h_{dot}", "h", "T_{atm}", "\\mu",
    #                 "\\rho", "W_{fuel}", "W_{N}", "W_{N+1}", "C_D", "C_L",
    #                 "\\eta_{prop}", "T", "h_{loss}", "P_{shaft-max}", "t",
    #                 "C_{f-fuse}", "C_{D-fuse}", "c_{dp}", "V_{wind}"]
    # Margins = ["BSFC", "c_{dp}"]
    Sens_boundaries = {"bad": 0.8, "good": 0.2}
    DF = mission_vars(M, sol)
    DF.to_csv(PATH + "mission_params.csv")
    # write_to_excel(PATH, "Mission_params.xlsx", DF, Sens_boundaries)
    DF = bd_vars(M, sol, "W", ["MTOW", "W_{fuel-tot}", "W_{zfw}"])
    write_to_excel(PATH, "W_breakdown.xlsx", DF, Sens_boundaries)

    m, mn = find_submodels([M], ["Mission"])
    for subm, name in zip(m, mn):
        df = model_params(subm, sol)
        if df is not None:
            df.to_csv(PATH + "%s.csv" % name)
