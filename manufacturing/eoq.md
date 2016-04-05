```python
#inPDF: skip
from gpkit import Variable
```
# Economic Order Quantity
Reference: Hopp, Wallace J., and Mark L. Spearman. Factory physics. Waveland Press, 2011.

Econonic order quantity (EOQ) expresses the fundamental tradeoff between
setup costs and holding costs. Small lot sizes have large setup costs per unit,
but large lot sizes incur holding costs.

The variables are:
```python
D = Variable("D", "count/year", "Demand rate")
c = Variable("c", "USD/count", "per unit production cost")
A = Variable("A", "USD", "setup cost")
h = Variable("h", "USD/count/year", "holding cost")
Q = Variable("Q", "count", "lot size")
Y = Variable("Y", "USD/year", "cost per year")
```

The cost per year is simply the sum of holding, setup, and production costs,
$$Y(Q) = \frac{hQ}{2} + \frac{AD}{Q} + cD$$
Or in gpkit,
```python
eoq = [Y >= h*Q/2 + A*D/Q + c*D]
```

Now create a model that minimizes cost
and substitute in the values given in Factory Physics.
```python
from gpkit import Model, units
m = Model(Y, eoq)
m.substitutions.update({D: 1000,  # /units("year"),
                        c: 250,   # *units("USD"),
                        A: 500,   # *units("USD"),
                        h: 35})   # *units("USD/year")})
sol = m.solve()
# would like to assert here that
# sol.subinto((2*A*D/h)**0.5 =~= sol.subinto(Q)
```
```python
#inPDF: replace with sol.generated.tex
with open("sol.generated.tex", "w") as f:
    f.write(m.solution.table(latex=True))
```
Now sweep $Q$ to see the tradoff curve from Factory Physics
```python
import numpy as np
m.substitutions.update({
    Q: ("sweep", np.logspace(np.log10(10), np.log10(500), 25))})
sol = m.solve()
```
```python
#inPDF: skip
import matplotlib.pyplot as plt
plt.figure()
plt.plot(sol("Q"), sol("Y")/sol("D") - sol("c"))
plt.ylim((0, 30))
plt.xlabel(Q)
plt.ylabel("Holding plus Setup Cost [$/unit]")
plt.savefig("qsweep.pdf")
plt.close()
```
![Cost vs Order Quantity](qsweep.pdf)
