from gpkit.small_scripts import unitstr

def gen_model_tex(model, modelname, texname=None):
    if texname:
        filename = texname
    else:
        filename = modelname
    with open('tex/%s.vars.generated.tex' % filename, 'w') as f:
        f.write("\\begin{longtable}{llll}\n \\toprule\n")
        f.write("\\toprule\n")
        f.write("Variables & Value & Units & Description \\\\ \n")
        f.write("\\midrule\n")
        #f.write("\\multicolumn{3}{l}\n")
        varnames = ["firstname"]
        for var in model.varkeys:
            name = var.name
            if name in varnames:
                pass
            else:
                if var.models[0] == modelname:
                    varnames.append(name)
                    unitstring = var.unitstr()[1:]
                    unitstring = "$[%s]$" % unitstring if unitstring else ""
                    val = "%0.3f" % var.value if var.value else ""
                    f.write("$%s$ & %s & %s & %s \\\\\n" %
                            (var.name, val, unitstring, var.label))
                else:
                    pass
        f.write("\\bottomrule\n")
        f.write("\\end{longtable}\n")

    with open('tex/%s.cnstrs.generated.tex' % filename, 'w') as f:
        l1 = model.latex(excluded=["models"]).replace("[ll]", "{ll}")
        models, modelnames = find_submodels([model], [])
        for m in modelnames:
            if m in l1:
                l1 = l1.replace("_{%s}" % m, "")
        lines = l1.split("\n")
        modeltex = "\n".join(lines[:1] + lines[3:])
        f.write("$$ %s $$" % modeltex)

def find_submodels(models, modelnames):
    runagain = 0
    for m in models:
        if "submodels" in m.__dict__.keys():
            for sub in m.submodels:
                if sub.__class__.__name__ not in modelnames:
                    models.append(sub)
                    modelnames.append(sub.__class__.__name__)
                    runagain += 1
                else:
                    pass
        else:
            pass
    if runagain > 0:
        return find_submodels(models, modelnames)
    else:
        return models, modelnames

def gen_tex_fig(fig, filename, caption=None):
    fig.savefig("figs/%s.pdf" % filename)
    with open("tex/%s.fig.generated.tex" % filename, "w") as f:
        f.write("\\begin{figure}[H]\n")
        f.write("\\label{f:%s}\n" % filename)
        f.write("\\begin{center}\n")
        f.write("\\includegraphics[scale=0.5]{figs/%s}\n" % filename)
        if caption:
            f.write("\\caption{%s}\n" % caption)
        f.write("\\end{center}\n")
        f.write("\\end{figure}\n")

def gen_fixvars_tex(model, solution, fixvars, filename=None):
    if filename:
        texname = "%s.table.generated.tex" % filename
    else:
        texname = "tex/fixvars.table.generated.tex"
    with open(texname, 'w') as f:
        f.write("\\begin{longtable}{lllll}\n \\toprule\n")
        f.write("\\toprule\n")
        f.write("Variables & Value & Units & Description \\\\ \n")
        f.write("\\midrule\n")
        for varname in fixvars:
            name = model[varname].descr["name"]
            val = "%0.3f" % solution(varname).magnitude
            unitstring = unitstr(model[varname].units)
            label = model[varname].descr["label"]
            if varname in solution["sensitivities"]["constants"]:
                sens = "%0.3f" % solution["sensitivities"]["constants"][varname]
            else:
                sens = ""
            f.write("$%s$ & %s & %s & %s & %s \\\\\n" %
                    (name, val, unitstring, sens, label))
        f.write("\\bottomrule\n")
        f.write("\\end{longtable}\n")
