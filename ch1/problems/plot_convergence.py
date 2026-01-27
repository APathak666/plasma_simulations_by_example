import pylab as pl
import numpy as np

pl.rcdefaults()
pl.rc('text', usetex=False)

data = np.loadtxt('convergence.csv', skiprows=1, delimiter=',')
iters = data[:,0]
residual = data[:,1]

pl.close('all')
fig, ax = pl.subplots(figsize=(7,4))

ax.semilogy(iters, residual, '-', linewidth=1.5, color='k')
ax.set_xlabel('iteration')
ax.set_ylabel('residual (log scale)')
ax.set_title('Solver Convergence')
ax.grid(True, which='both')
pl.tight_layout()
pl.show()
