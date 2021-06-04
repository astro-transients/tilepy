import matplotlib.pyplot as plt
import os
import numpy as np
import astropy.units as u
from gammapy.spectrum.models import Absorption
import sys
sys.path.append('/Users/mseglar/Documents/CurrentPhD/CTA/grb_paper-master/scripts')
from GWCTAPointingTools import GRB
import sys
sys.path.append('/Users/mseglar/Documents/GitHub/GWfollowup/include')
from GWHESSPointingTools import *
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
filepath = '/Users/mseglar/Documents/CurrentPhD/CTA/GWCosmos/GammaCatalogV1.0/run0017_ID000132.fits'

absorption = Absorption.read_builtin('dominguez')
grb = GRB.from_fitsfile(filepath=filepath, absorption=absorption)

for i in range(0, len(grb.spectral_model)):
    grb.spectral_model[i].plot(energy_range=(0.001, 10) * u.TeV)
    #grb.spectral_model[i].plot(time_range=(0.001, 10) * u.s)
#plt.show()

data = fits.open(filepath)
energy_interval = data[1].data.field(0)
time_interval = data[2].data.field(0)
flux = data[3].data
Flux=[]
for i in range(0,len(energy_interval)):
        flux = data[3].data.field(i)
        Flux.append(flux)
Flux=np.array(Flux)
print(time_interval.ndim)
print(energy_interval.ndim)
print(Flux.ndim)

#flux1 = data[3].data.field(1)

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot the surface.
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
surf = ax.plot_surface(X,Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
