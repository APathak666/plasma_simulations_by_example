import pylab as pl
import numpy as np

pl.rcdefaults()
pl.rc('text', usetex=False)

data = np.loadtxt('ef.csv', skiprows=1, delimiter=',')

x = data[:,0]
y = data[:,1]
z = data[:,2]
ef_x = data[:,3]
ef_y = data[:,4]
ef_z = data[:,5]

# get unique coordinates
xu = np.unique(x)
yu = np.unique(y)
zu = np.unique(z)

ni, nj, nk = len(xu), len(yu), len(zu)

# reshape into 3D arrays (data is stored in i,j,k order with k varying fastest)
efx_3d = ef_x.reshape(ni, nj, nk)
efy_3d = ef_y.reshape(ni, nj, nk)
efz_3d = ef_z.reshape(ni, nj, nk)
ef_mag = np.sqrt(efx_3d**2 + efy_3d**2 + efz_3d**2)

jmid = nj // 2
kmid = nk // 2

pl.close('all')
fig, axes = pl.subplots(1, 3, figsize=(16, 4))

# slice along x at y=mid, z=mid
ax = axes[0]
ax.plot(xu, efx_3d[:, jmid, kmid], '-', linewidth=2, color='r', label='ef_x')
ax.plot(xu, efy_3d[:, jmid, kmid], '--', linewidth=2, color='g', label='ef_y')
ax.plot(xu, efz_3d[:, jmid, kmid], ':', linewidth=2, color='b', label='ef_z')
ax.set_xlabel('x (m)')
ax.set_ylabel('electric field (V/m)')
ax.set_title('E-field along x (y=mid, z=mid)')
ax.legend()
ax.grid(True)

# 2D contour: |E| in x-y plane at z=mid
ax = axes[1]
X, Y = np.meshgrid(xu, yu, indexing='ij')
c = ax.contourf(X, Y, ef_mag[:, :, kmid], levels=20, cmap='hot')
pl.colorbar(c, ax=ax, label='|E| (V/m)')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('|E| in x-y plane (z=mid)')

# 2D quiver: E-field vectors in x-y plane at z=mid
ax = axes[2]
skip = max(1, ni // 10)  # subsample for readability
ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
          efx_3d[::skip, ::skip, kmid], efy_3d[::skip, ::skip, kmid],
          ef_mag[::skip, ::skip, kmid], cmap='hot')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('E-field vectors in x-y plane (z=mid)')
ax.set_aspect('equal')

pl.tight_layout()
pl.show()
