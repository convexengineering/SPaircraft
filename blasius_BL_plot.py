import numpy as np
import matplotlib.pyplot as plt

#values from http://www.calpoly.edu/~kshollen/ME347/Handouts/Blasius.pdf
x = [0, 0.1659, 0.3298, 0.4868, 0.6298, 0.7513, 0.8461, 0.9131, 0.9555, 0.9795, 0.9916]
y = np.linspace(0,5,11)

plt.plot(x,y)
plt.ylabel('$\\frac{u}{U}$', fontsize=24)
plt.xlabel('$\\eta = y \\sqrt{\\frac{U}{\\nu x}}$', fontsize=24)
plt.savefig('blasius_BL.pdf', bbox_inches="tight")
plt.show()

print sum(x)/11
