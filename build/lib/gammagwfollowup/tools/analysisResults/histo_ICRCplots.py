import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import seaborn as sns
import seaborn
from scipy import stats
import numpy as np
from scipy import stats
import statistics
from matplotlib.mlab import griddata
def func(x, a, b, c):
    return a * np.sin(b * x)+c

filename=('/Users/mseglar/Documents/GitLab/gw-follow-up-simulations/tools/analysisResults/SimuGW_OnAxis.txt')
Theta,Area,luminosity, durationPointing,scheduledPointings, possiblePointings,foundFirst,foundTimes,pGW= np.genfromtxt(filename, usecols=(3,4,5,6,7,8,9,10,11),unpack=True,skip_header=1)

print(Theta,Area,luminosity, durationPointing,scheduledPointings, possiblePointings,foundFirst,foundTimes,pGW)


print("--------------------- Statistics -------------")

V2=[]

V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity<10**47.5)])/len(foundFirst[luminosity<10**47.5]))
V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity>10**47.5)&(luminosity<10**48.5)])/len(foundFirst[(luminosity>10**47.5)&(luminosity<10**48.5)]))
V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity>10**48.5)&(luminosity<10**49.5)])/len(foundFirst[(luminosity>10**48.5)&(luminosity<10**49.5)]))
V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity>10**49.5)&(luminosity<10**50.5)])/len(foundFirst[(luminosity>10**49.5)&(luminosity<10**50.5)]))
V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity>10**50.5)&(luminosity<10**51.5)])/len(foundFirst[(luminosity>10**50.5)&(luminosity<10**51.5)]))
V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity>10**51.5)&(luminosity<10**52.5)])/len(foundFirst[(luminosity>10**51.5)&(luminosity<10**52.5)]))
V2.append(100*len(foundFirst[(foundFirst>-1)&(luminosity>10**52.5)])/len(foundFirst[luminosity>10**52.5]))
'''
#Cumulative
V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**47)])/len(foundFirst2[luminosity2<10**47]))
V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**48)])/len(foundFirst2[luminosity2<10**48]))
V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**49)])/len(foundFirst2[luminosity2<10**49]))
V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**50)])/len(foundFirst2[(luminosity2<10**50)]))
V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**51)])/len(foundFirst2[(luminosity2<10**51)]))
V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**52)])/len(foundFirst2[(luminosity2<10**52)]))
#V2.append(100*len(foundFirst2[(foundFirst2>-1)&(luminosity2<10**53)])/len(foundFirst2[luminosity2>10**53]))


fig, ax = plt.subplots(figsize =(9, 7))
sns.violinplot(ax = ax, data = [foundFirst[(foundFirst>-1)&(Area<31)],foundFirst[(foundFirst>-1)&(Area>31)&(Area<316)],foundFirst[(foundFirst>-1)&(Area>316)]])
plt.grid()
plt.savefig('./Plots/test.png')

'''
bin=2

fig, ax = plt.subplots()
plt.hist(pGW,100,cumulative=True,density=True,histtype='step',color='limegreen',label='% found')
ax.set_ylabel('%')
ax.legend(loc='upper left')
plt.grid()
plt.savefig('./Plots/PGW_cum.png')

fig, ax = plt.subplots()
plt.hist(luminosity,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('#')
ax.set_xlabel('E$_{iso}$ [erg]')
plt.grid()
plt.savefig('./Plots/Eiso.png')

fig, ax = plt.subplots()
plt.hist(durationPointing,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
ax.set_yscale('log')
ax.set_ylabel('#')
ax.set_xlabel('Duration [s]')
plt.grid()
plt.savefig('./Plots/durationPointing.png')

fig, ax = plt.subplots()
plt.hist(foundFirst,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
plt.hist(scheduledPointings,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
plt.hist(possiblePointings,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
plt.hist(foundTimes,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
#plt.hist((found3[found2>-1]+1)-(found2[found2>-1]+1),bin, histtype='step', stacked=True,  fill=False,label='3D-2D')
ax.set_ylabel('#')
ax.set_yscale('log')
ax.set_xlabel('Time interval')
plt.grid()
plt.savefig('./Plots/aboutPointings.png')


fig, ax = plt.subplots()
plt.hist(pGW,bin, histtype='step', stacked=True,  fill=False,color='dodgerblue')
ax.set_ylabel('#')
ax.set_xlabel('Covered probability')
plt.grid()
plt.savefig('./Plots/durationPointing.png')


fig, ax = plt.subplots()
plt.plot(Area,foundFirst,'+',color='dodgerblue')
ax.set_ylabel('Number of pointing when source covered')
ax.set_xlabel('90% C.R. [deg$^2$]')
plt.grid()
plt.savefig('./Plots/Area_firstCovered.png')

fig, ax = plt.subplots()
plt.plot(foundFirst,durationPointing,'+',color='dodgerblue')
ax.set_xlabel('Number of pointing when source covered')
ax.set_ylabel('Duration [s]')
plt.grid()
plt.savefig('./Plots/Area_firstCovered.png')


fig, ax = plt.subplots()
plt.errorbar(np.arange(len(V2)),V2,xerr=0.5,fmt='o',linestyle='dotted' , color='dodgerblue')
plt.plot(np.arange(len(V2)),V2,'o',color='dodgerblue')
ax.set_ylabel('%')
ax.set_xlabel('E$_{iso}$ [erg]')
ax.set_xticklabels(['10$^{46}$','10$^{47}$','10$^{48}$','10$^{49}$','10$^{50}$','10$^{51}$','10$^{52}$','10$^{53}$'])
#ax.legend()
plt.grid()
plt.savefig('./Plots/Percentage_comparison.png')

'''
fig, ax = plt.subplots()

plt.errorbar(np.arange(len(V2)),V2,xerr=0.5,fmt='o',linestyle='dotted' , color='dodgerblue')
plt.plot(np.arange(len(V2)),V2,'o',color='dodgerblue')
#popt, pcov = curve_fit(func, np.arange(len(V2)), V2,p0=[100, 0.6,2])
#plt.plot(np.arange(0,len(V2),0.1), func(np.arange(0,len(V2),0.1), *popt), 'r-',label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
#plt.plot(np.arange(len(V3)),V3,color='dodgerblue',label='% found 2D')
ax.set_ylabel('%')
ax.set_xlabel('E$_{iso}$ [erg]')
ax.set_xticklabels(['10$^{46}$','10$^{47}$','10$^{48}$','10$^{49}$','10$^{50}$','10$^{51}$','10$^{52}$','10$^{53}$'])
#ax.legend()
plt.grid()
plt.savefig('./Plots/Percentage_comparison.png')



fig, ax = plt.subplots()
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1]),bin, histtype='step', stacked=True,  fill=False,label='2D')
#plt.hist(np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
plt.plot(luminosity2,Theta2,'+',color='dodgerblue')
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1])-np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
#ax.set_yscale('log')
ax.set_xscale('log')
print(luminosity2)
ax.set_xlabel('${theta}$ [deg]')
ax.set_xlabel('E$_{iso}$ [erg]')
plt.grid()
plt.savefig('./Plots/Eiso_theta.png')


bin=100
fig, ax = plt.subplots()
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1]),bin, histtype='step', stacked=True,  fill=False,label='2D')
#plt.hist(np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
plt.plot(luminosity2,TotalPointingsBarbara2,'+',color='dodgerblue')
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1])-np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(5*10e47,3*10e52)
ax.set_ylabel('Number of $\t{potential}$ scheduled observations')
ax.set_xlabel('E$_{iso}$ [erg]')
plt.grid()
plt.savefig('./Plots/Eiso_TotalPointings.png')


bin=100
fig, ax = plt.subplots()
#plt.plot(luminosity2,scheduledPointings2,'+',label='2D scheduledPointings')
plt.plot(Theta2,foundFirst2,'+',label='2D foundFirst')
#plt.plot(luminosity3,scheduledPointings3,'+',label='3D scheduledPointings',color='dodgerblue')
plt.plot(Theta3,foundFirst3,'+',label='3D foundFirst')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('Number of pointing when source covered')
ax.set_xlabel('${theta}$ [deg]')
ax.legend()
plt.grid()
plt.savefig('./Plots/FoundFirst_Theta.png')

bin=100
fig, ax = plt.subplots()
#plt.plot(luminosity2,scheduledPointings2,'+',label='2D scheduledPointings')
plt.plot(luminosity2,foundFirst2,'+',label='2D foundFirst')
#plt.plot(luminosity3,scheduledPointings3,'+',label='3D scheduledPointings',color='dodgerblue')
plt.plot(luminosity3,foundFirst3,'+',label='3D foundFirst')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('Number of pointing when source covered')
ax.set_xlabel('E$_{iso}$ [erg]')
ax.legend()
plt.grid()
plt.savefig('./Plots/FoundFirst_Lum.png')

bin=100
fig, ax = plt.subplots()
#plt.plot(luminosity2,scheduledPointings2,'+',label='2D')
plt.plot(luminosity2,foundTimes2,'+',label='2D',color='limegreen')
#plt.plot(luminosity3,scheduledPointings3,'+',label='3D')
plt.plot(luminosity3,foundTimes3,'+',label='3D',color='dodgerblue')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('Number of times that is source covered')
ax.set_xlabel('E$_{iso}$ [erg]')
ax.legend()
plt.grid()
plt.savefig('./Plots/TimesFoundPointings_Lum.png')



bin=100
fig, ax = plt.subplots()
#plt.plot(luminosity2,scheduledPointings2,'+',label='2D scheduledPointings')
x=[1,10,100,1000]
yerr=[statistics.stdev(foundFirst2[(foundFirst2>-1)&(Area2<3)]),statistics.stdev(foundFirst2[(foundFirst2>-1)&(Area2>3)&(Area2<31)]),statistics.stdev(foundFirst2[(foundFirst2>-1)&(Area2>31)&(Area2<316)]),statistics.stdev(foundFirst2[(foundFirst2>-1)&(Area2>316)])]
y=[np.median(foundFirst2[(foundFirst2>-1)&(Area2<3)]),np.median(foundFirst2[(foundFirst2>-1)&(Area2<31)]),np.median(foundFirst2[(foundFirst2>-1)&(Area2>31)&(Area2<316)]),np.median(foundFirst2[(foundFirst2>-1)&(Area2>316)])]
plt.plot(Area2[foundFirst2>-1],foundFirst2[foundFirst2>-1],'+',label='P$_{GW-in-FoV}$ scheduler')
plt.errorbar(x,y,yerr=yerr,fmt='o',label='Binned median$\pm \sigma$')
#plt.plot(luminosity3,scheduledPointings3,'+',label='3D scheduledPointings',color='dodgerblue')
#plt.plot(Area3,foundFirst3,'+',label='3D foundFirst')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('Number of pointing when source covered')
ax.set_xlabel('90% C.R. [deg$^2$]')
ax.legend()
plt.grid()
plt.savefig('./Plots/FoundFirst_Area.png')

bin=100
fig, ax = plt.subplots()
#plt.plot(luminosity2,scheduledPointings2,'+',label='2D')
plt.plot(Area2,foundTimes2,'+',label='2D',color='limegreen')
#plt.plot(luminosity3,scheduledPointings3,'+',label='3D')
plt.plot(Area3,foundTimes3,'+',label='3D',color='dodgerblue')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('Number of times that is source covered')
ax.set_xlabel('90% C.R. [deg$^2$]')
ax.legend()
plt.grid()
plt.savefig('./Plots/TimesFoundPointings_Area.png')



bin=100
fig, ax = plt.subplots()
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1]),bin, histtype='step', stacked=True,  fill=False,label='2D')
#plt.hist(np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
plt.plot(luminosity2,TimeFoundFirst2,'+',label='2D',color='limegreen')
plt.plot(luminosity3,TimeFoundFirst3,'+',label='3D',color='dodgerblue')
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1])-np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Covered at time')
ax.set_xlabel('E$_{iso}$ [erg]')
ax.legend()
plt.grid()
plt.savefig('./Plots/Times_Lum.png')

bin=100
fig, ax = plt.subplots()
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1]),bin, histtype='step', stacked=True,  fill=False,label='2D')
#plt.hist(np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
plt.hist(TimeFoundFirst2,100,density=True,histtype='step', stacked=True,  fill=False,color='dodgerblue')
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1])-np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Covered at time')
ax.set_xlabel('E$_{iso}$ [erg]')
ax.legend()
plt.grid()
plt.savefig('./Plots/Times_hist.png')

########################################################################################
########################################################################################
########################################################################################
#Import other files!
filenameObsTimes='/Users/mseglar/Documents/CurrentPhD/CTA/AllObsTimes.txt'
ID,duration= np.genfromtxt(filenameObsTimes, usecols=(0,3),unpack=True)
print(ID[0])
bin=100
fig, ax = plt.subplots()
#plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1]),bin, histtype='step', stacked=True,  fill=False,label='2D')
#plt.hist(np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
plt.hist(duration,100)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('#')
ax.set_xlabel('Duration [s]')
ax.legend()
plt.grid()
plt.savefig('./Plots/HistObservationsDurations.png')
'''
'''
fig, ax = plt.subplots()

plt.plot(0,max(duration[[ID==1.0]]),'+',color='green',label='Max')
plt.plot(0,np.mean(duration[[ID==1.0]]),'+',color='dodgerblue',label='Mean')
plt.plot(0,np.median(duration[[ID==1.0]]),'+',color='firebrick',label='Median')
for i in range(2,int(max(ID))):
    plt.plot(i,max(duration[[ID==i]]),'+',color='green')
    plt.plot(i,np.mean(duration[[ID==i]]),'+',color='dodgerblue')
    plt.plot(i,np.median(duration[[ID==i]]),'+',color='firebrick')
    #plt.hist(np.divide(found2[found2>-1]+1,pointings2[found2>-1])-np.divide(found3[found3>-1]+1,pointings3[found3>-1]),bin, histtype='step', stacked=True,  fill=False,label='3D')
ax.set_yscale('log')
ax.set_ylabel('Duration [s]')
ax.set_xlabel('Number of $\t{potential}$ scheduled observations')
ax.legend()
plt.grid()
plt.savefig('ICRCSimulations_OnAxis/Plots/MeanObservationsDurations_NumberObs.png')
'''

