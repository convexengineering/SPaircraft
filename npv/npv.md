# NPV Model
```python
#inPDF: skip

from numpy import pi
from gpkit import VectorVariable, Variable, Model, units
import gpkit
from gpkit.tools import te_exp_minus1

class NPV(Model):
    def __init__(self, **kwargs):


```
# Variables

```python
#inPDF: replace with vartable.generated.tex
        N = 3 # number of payments
        npv = Variable('NPV', 'USD', 'net present value')
        PV = Variable('PV', 'USD', 'present value')
        C = VectorVariable(N, 'C', 'USD', 'cash flow')
        C_r = Variable('C_r', 5e6, 'USD/years', 'cash flow rate')
        r = Variable('r', 0.1, '1/years', 'interest rate')
        t = VectorVariable(N, 't', 'years', 'time')

```

# Problem Set Up

Let's assume that we want the Net Present Value (NPV) to be \$10 Million.  NPV can be expressed by:

$$ \text{NPV} = \displaystyle\sum\limits_{i=0}^N \text{PV}_i $$ .

Let's assume that the time periods of the payments are not equal but that the PV$_{i}$ are equal and that there are 3 payments.  PV can be expressed by

$$\text{PV} = \text{C} e^{-rt} $$.

With these assumptions we can get rid of the $i$ subscript and claim that

$$\text{NPV} \leq 3\text{PV} $$

Now the problem becomes solving for the length of each time period such that each PV is equal.  Let's assume that the same payment, C, is made at each payment period and that C is given.  This allows us to write

$$ \text{C} \leq \text{PV}e^{r t_i} $$

Let's assume that $t_0$ is the time of evaulation of the NPV and that $t_i$ is when every payment is made and every PV evaluated. This means that

$$ \begin{bmatrix} \Delta t_1 = t_1 - t_0 \\
\Delta t_2 = t_2 - t_1 \\
\Delta t_3 = t_3 - t_2
\end{bmatrix} $$

If there is a cash flow rate, C$_r$ then,

$$ \text{C} \leq\text{C}_r \Delta t $$


# Constraints and Objective

```python
#inPDF: replace with model.generated.tex

        constraints = [N*PV >= npv,
                       C >= PV*(te_exp_minus1(r*t, 4) + 1),
                       t >= t.left + C/C_r]

        cost = 1/npv
```
```python
#inPDF: replace with sol.generated.tex

        Model.__init__(self, cost, constraints, **kwargs)

    @classmethod
    def test(cls):
        "test the model by solving it"
        return cls().solve()

if __name__ == "__main__":
    SOL = NPV.test()
    with open("vartable.generated.tex", "w") as f:
        f.write("\\begin{tabbing}\n  XXXXXXXXXX \\= \\kill\n")
        for var in M.varkeys:
            f.write("$%s$ : [%s] %s \\\\\n" % (var.name, var.units, var.label))
	f.write("\\end{tabbing}")
    with open("model.generated.tex", "w") as f:
        f.write("$$ %s $$" % M.latex(excluded=["models"]).replace("[ll]", "{ll}"))
    with open("sol.generated.tex", "w") as f:
        f.write(SOL.table(latex=True))
```
