import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm


def FD2(ur, delta_x, ax):
    du_dx = (np.roll(ur, -1, axis=ax) - ur) / delta_x
    return du_dx


def SD3(ur, delta_x, ax):
    du_dx_2 = (np.roll(ur, -1, axis=ax) - 2 * ur + np.roll(ur, 1, axis=ax)) / delta_x ** 2
    return du_dx_2


def BC(u):
    r = np.zeros(np.shape(u[0, :, :]))
    ny = np.squeeze(np.shape(u[0, :, 0]))
    nz = np.squeeze(np.shape(u[0, 0, :]))

    yc = dy * (ny - 1) / 2  # Center of y (R)
    zc = dz * (nz - 1) / 2  # Center of z (R)
    for i in range(ny):
        for j in range(nz):
            y = i * dy
            z = j * dz
            r[i, j] = ((y - yc) ** 2 + (z - zc) ** 2) ** 0.5  # Calculate the distance from center
            if r[i, j] >= dy * (ny - 1) / 2:  # If distance >R, set u to zero
                u[:, i, j] = 0
    return u

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
t = np.arange(0, 1 + dt, dt)

x = 1               # 血管長度
P1 = 15998.6842         # 初始壓力
P = np.zeros((n, n, n))
dx = x / n
for i in range(n):
    P[i, :, :] = P1 - 666.611842 * dx * i

dP_dx = FD2(P, dx, 0)

ur = np.zeros([n, n, n])
for i in range(100):
    ur = BC(ur)
    ur[i, :, :] = 1/(4*mu) * (r**2-R**2) * dP_dx[i, :, :]


plt.title("The velocity of the blood [analytic solution]")
a = plt.contourf(ur[50, :, :], cmap=cm.jet, levels=np.arange(0, 9, 0.1))
plt.colorbar(a)
plt.show()