import matplotlib.pyplot as plt
import astropy.units as u
from gammapy.spectrum.models import Absorption
from gammagwfollowup.include.GWCTAPointingTools import GRB
from gammagwfollowup.include.GWHESSPointingTools import *
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import interpolate
from scipy import integrate

def spectrum(x):
    return (x/1)**(-2.2)


filepath = '/Users/mseglar/Documents/GitLab/dataset/GammaCatalogV2.0/run0378_ID000234.fits'

tpointingFile = '/Users/mseglar/Documents/GitLab/output/PGWonFoV/ScheduledObs/run0378_MergerID000234.txt'
ObsDelay = np.genfromtxt(tpointingFile, usecols=(10), skip_header=1,delimiter=' ',unpack=True)  # ra, dec in degrees


hdu_list = fits.open(filepath)
#print(hdu_list[3].data.field(0))

datalc = hdu_list[3].data
datatime = hdu_list[2].data
dataenergy = hdu_list[1].data

lc = datalc.field(0)
time = datatime.field(0)
energy = dataenergy.field(0)

E1 = 30
E2 = 10000
# Integral of the spectrum
intl, errl = integrate.quad(lambda x: spectrum(x), E1, E2)  # GeV

flux = intl * lc
# Interpolation of the flux with time
#flux = interpolate.interp1d(time, lc)
# fluencen, errorn = integrate.quad(lambda x: flux(x), tstart, t)  # ph/cm2/GeV

#energy_interval = hdu_list[1].data.field(0)
#time_interval = hdu_list[2].data.field(0)
#flux = hdu_list[3].data
#Flux=[]


plt.plot(time,flux,color = 'k', label = 'Simulated GRB emission (30 GeV-10 TeV) ')
plt.axvline(x=210, color='red', linestyle='dashed') #Alert is recevied
plt.axvline(x=238, color='red', linestyle='dashed')  # Telescopes moving
for i in range(0,len(ObsDelay)-2):
    plt.axvline(x=ObsDelay[i], color = 'orange')
    plt.axvline(x=ObsDelay[i]+3, color = 'orange')
    plt.axvline(x=ObsDelay[i]+5, color = 'orange')
    plt.axvline(x=ObsDelay[i]+8, color = 'orange')
    plt.axvline(x=ObsDelay[i]+10, color = 'orange')
plt.yscale('log')
plt.xscale('log')
plt.ylim([0.00000000001,0.00001])
plt.xlim([10,1000])
plt.ylabel('Flux [ph cm$^{-2}$ s$^{-1}$]')
plt.xlabel('Time [s]')
plt.grid()
plt.legend(loc = 'upper left')
plt.savefig("run0378_ID000234_LC.png")

#absorption = Absorption.read_builtin('dominguez')
#grb = GRB.from_fitsfile(filepath=filepath, absorption=absorption)

#for i in range(0, len(grb.spectral_model)):
#    grb.spectral_model[i].plot(energy_range=(0.001, 10) * u.TeV)
    #grb.spectral_model[i].plot(time_range=(0.001, 10) * u.s)
#plt.show()

'''
intl, errl = integrate.quad(lambda x: spectrum(x), x1, x2)  # GeV

# Interpolation of the flux with time
flux = interpolate.interp1d(time, lc)

for i in range(0,len(energy_interval)):
        flux = data[3].data.field(i)
        Flux.append(flux)
Flux=np.array(Flux)
print(time_interval.ndim)
print(energy_interval.ndim)
print(Flux.ndim)
'''

'''
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
'''