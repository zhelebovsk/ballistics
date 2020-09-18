import numpy as np
from matplotlib import pyplot as plt


def Fdf(u, v, nu, d):
    Rep = np.abs(u - v) * d / nu
    if Rep > 1000.0:
        C = 0.11 * Rep / 6.0
    else:
        C = 1 + np.power(Rep, 2.0 / 3.0) / 6.0
    Cd = 24.0 / Rep * C
    return Cd * rol * np.pi * np.power(d, 2.0) / 4.0 * (u - v) * np.abs(u - v) / 2.0


def Fgf(rop, d, g):
    return - rop * np.pi * np.power(d, 3.0) * g / 6.0

# dydx = (q(la, k) - q(la + num, k)) / num

def ueq(u, v, nu, d, rop, g):
    num = 0.00000000000001
    la = u
    while np.abs(Fgf(rop, d, g) + Fdf(la, v, nu, d)) > num:
        dydx = (Fdf(la + num, v, nu, d) - (Fdf(la, v, nu, d))) / num
        la = la - (Fdf(la, v, nu, d) + Fgf(rop, d, g)) / dydx
    return la


g = 9.81
d = np.arange(100.0/1000000.0, 1000.0/1000000.0, 10.0/1000000.0)       # m
u = 1.0                    # m/s
mu = 1.78 / 100000.0         # dynamical viscosity
rol = 1.24                   # density
nu = mu / rol
rop = 1500.0
v = 0.0
u = np.zeros([np.size(d)])

k = 0

for i in d:
    u[k] = ueq(1.0, v, nu, i, rop, g)
    print(ueq(1.0, v, nu, i, rop, g))
    k = k + 1
plt.plot(u, d)
plt.show()

