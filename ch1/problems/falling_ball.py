import matplotlib.pyplot as plt
import numpy as np

r = 0.25
rho = 2700
h = 10_000
g = -9.81
dt = 1e-5
cd = 0.47
rho_air = 1.225
A = 3.14159 * r**2
m = (4/3) * 3.14159 * r**3 * rho
k = 0.5 * rho_air * cd * A
v_terminal = (m * abs(g) / k) ** 0.5

y = h
v = 0
v -= g*dt/2
t = 0

v_coeff = rho_air * cd * A * dt / (2 * m)

ts = []
ys_sim = []
ys_analytical = []
vs_sim = []
vs_analytical = []

print(f"Terminal velocity: {v_terminal:.2f} m/s")

for step in range(3_0000000):
    ts.append(t)
    ys_sim.append(y)
    vs_sim.append(v)

    # v_analytical = -v_terminal * np.tanh(abs(g) * t / v_terminal)
    # y_analytical = h - (v_terminal**2 / abs(g)) * np.log(np.cosh(abs(g) * t / v_terminal))

    # ys_analytical.append(y_analytical)
    # vs_analytical.append(v_analytical)

    v_old = v
    v += g*dt
    v +=  v_old * v_old * v_coeff * (-1 if v_old > 0 else 1)
    y += v*dt
    t += dt

    if y < 0:
        y = -y
        v = -v

fig, axs = plt.subplots(2, 1, sharex=True)

axs[0].plot(ts, ys_sim)
# axs[0].plot(ts, ys_analytical, "--", label="y (analytical)")
axs[0].set_ylabel("Position y (m)")
axs[0].set_title("Falling Ball Simulation")

axs[1].plot(ts, vs_sim)
# axs[1].plot(ts, vs_analytical, "--", label="v (analytical)")
axs[1].set_ylabel("Velocity v (m/s)")
axs[1].set_xlabel("Time (s)")

plt.tight_layout()
plt.show()
