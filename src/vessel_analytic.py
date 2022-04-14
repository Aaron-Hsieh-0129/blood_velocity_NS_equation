import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def FD2(ur, delta_x, ax):
    du_dx = (np.roll(ur, -1, axis=ax) - ur) / delta_x
    return du_dx


def SD3(ur, delta_x, ax):
    du_dx_2 = (np.roll(ur, -1, axis=ax) - 2 * ur + np.roll(ur, 1, axis=ax)) / delta_x ** 2
    return du_dx_2


# initial value
mu = 3 * 10 ** (-3)
rho = 1.06 * 10 ** 3
nu = mu / rho

n = 101                 # 切幾等分
R = 0.0125              # 血管半徑
y = np.linspace(-R, R, n)
dy = y[1] - y[0]
z = np.linspace(-R, R, n)
dz = z[1] - z[0]
r = np.zeros([101, 101])

for i in range(101):
    for j in range(101):
        r[i, j] = np.sqrt(y[i] ** 2 + z[j] ** 2)

r = np.where(r <= R, r, 0)
dt = 0.001

x = 0.5                 # 血管長度
P1 = 15998.6842         # 初始壓力
P = np.zeros((n, n, n))
dx = x / n
for i in range(n):
    P[i, :, :] = P1 - 50 * dx * i

dP_dx = FD2(P, dx, 0)

ur = np.zeros([n, n, n])
for i in range(100):
    ur[i, :, :] = 1/(4*mu) * (r**2-R**2) * dP_dx[i, :, :]


plt.xlabel("The length of the blood vessel (split)")
plt.ylabel("velocity")
plt.title("The velocity of the blood [analytic solution]")
plt.plot(ur[50, 50, :])
plt.show()