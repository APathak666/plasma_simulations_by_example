import pylab as pl
import numpy as np

pl.rcdefaults()
pl.rc('text', usetex=False)

data = np.loadtxt('potential.csv', skiprows=1, delimiter=',')

x = data[:,0]
y = data[:,1]
z = data[:,2]
phi = data[:,3]
rho = data[:,4]

# get unique coordinates
xu = np.unique(x)
yu = np.unique(y)
zu = np.unique(z)

ni, nj, nk = len(xu), len(yu), len(zu)

# reshape into 3D arrays (data is stored in i,j,k order with k varying fastest)
phi_3d = phi.reshape(ni, nj, nk)
rho_3d = rho.reshape(ni, nj, nk)

# take midplane slices
jmid = nj // 2
kmid = nk // 2

pl.close('all')
fig, axes = pl.subplots(1, 3, figsize=(16, 4))

# slice along x at y=mid, z=mid
ax = axes[0]
ax.plot(xu, phi_3d[:, jmid, kmid], '-', linewidth=2, color='k', label='phi')
ax.set_xlabel('x (m)')
ax.set_ylabel('potential (V)')
ax.set_title('Potential along x (y=mid, z=mid)')
ax.legend()
ax.grid(True)

# 2D contour: phi in x-y plane at z=mid
ax = axes[1]
X, Y = np.meshgrid(xu, yu, indexing='ij')
c = ax.contourf(X, Y, phi_3d[:, :, kmid], levels=20, cmap='RdBu_r')
pl.colorbar(c, ax=ax, label='phi (V)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Potential in x-y plane (z=mid)')

# 2D contour: rho in x-y plane at z=mid
ax = axes[2]
c = ax.contourf(X, Y, rho_3d[:, :, kmid], levels=20, cmap='RdBu_r')
pl.colorbar(c, ax=ax, label='rho (C/m^3)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('Charge density in x-y plane (z=mid)')

pl.tight_layout()
pl.show()
