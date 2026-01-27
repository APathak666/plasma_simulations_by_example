import pylab as pl
import numpy as np

pl.rcdefaults()
pl.rc('text', usetex=False)

# read ef.csv produced by logElectricField (columns: x, ef, ef_a)
data = np.loadtxt('ef.csv', delimiter=',', skiprows=1)

x = data[:, 0]
ef = data[:, 1]
ef_a = data[:, 2]

pl.close('all')
fig, ax = pl.subplots(figsize=(7, 4))

ax.plot(x, ef, '-', linewidth=2, color='k', label="electric field (ef)")
ax.plot(x, ef_a, '--', linewidth=2, color='0.5', label="analytic ef")

ax.set_xlabel('position (m)')
ax.set_ylabel('electric field (V/m)')
ax.legend(loc='best')
ax.grid(axis='both')

pl.tight_layout()
pl.show()
