import pylab as pl
import numpy as np

pl.rcdefaults()
pl.rc('text', usetex=False)

#read results.txt (tab-separated: x, phi, rho, ef)
data = np.loadtxt('results.txt', skiprows=1)

x = data[:,0]
phi = data[:,1]  # electric potential (V)
ef = data[:,3]   # electric field

#create a new figure
pl.close('all')
fig, ax1 = pl.subplots(figsize=(7,4))

#plot potential and electric field vs x
ax1.plot(x,phi,'-',linewidth=2,color='k',label="potential (phi)")
ax1.plot(x,ef,'--',linewidth=2,color='0.5',label="electric field")

ax1.set_xlabel('position (m)')
ax1.set_ylabel('potential (V)')
ax1.legend(loc='best')
pl.tight_layout()

ax1.grid(axis='both')
pl.show()
