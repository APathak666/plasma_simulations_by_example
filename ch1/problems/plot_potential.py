import pylab as pl
import numpy as np

pl.rcdefaults()
pl.rc('text', usetex=False)

#read results.txt (tab-separated: x, phi, rho, ef)
data = np.loadtxt('potential.csv', skiprows=1, delimiter=',')

x = data[:,0]
phi = data[:,1]  # electric potential (V)
phi_a = data[:,2]   # electric field

#create a new figure
pl.close('all')
fig, ax1 = pl.subplots(figsize=(7,4))

#plot potential and electric field vs x
ax1.plot(x,phi,'-',linewidth=2,color='k',label="potential (phi)")
ax1.plot(x,phi_a,'--',linewidth=2,color='0.5',label="potential (analytical)")

ax1.set_xlabel('position (m)')
ax1.set_ylabel('potential (V)')
ax1.legend(loc='best')
pl.tight_layout()

ax1.grid(axis='both')
pl.show()
