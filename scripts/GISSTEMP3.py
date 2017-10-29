
import numpy as np
from pylab import plot, show, bar, legend, colors, axes, xlabel, ylabel
from pylab import title, savefig, axis, figure, semilogy, mean, exp, sqrt
from pylab import log, arctan, sum
import netCDF4
from scipy import stats
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits import basemap
from dTgiss_ann import dTgiss_ann
from pop import pop


# Get GISS dT data
years_giss,lats_giss,lons_giss,temp_mn,tmp=dTgiss_ann()

# Get population data
years,lats,lons,pop_tot,wt_pop=pop()
lons, lats = np.meshgrid(lons, lats) 
npops=len(years)

# Regrid temperature onto population grid
ntime=len(years_giss)
temp_pop_f=np.zeros(ntime)
temp_pop_i=np.zeros(ntime)
temp_area=np.zeros(ntime)
weight_area=np.cos(lats*np.pi/180.0)
# Calculate mean temperature exposure for human population
for itime in range(0, ntime):
    dTgiss=tmp[itime,:,:]
    dT = basemap.interp(dTgiss,lons_giss,lats_giss, \
                     lons,lats,checkbounds=False,masked=False,order=1)
    weight_pop_i=wt_pop[0,:,:]   
    weight_pop_f=wt_pop[npops-1,:,:]  
    temp_pop_i[itime]=sum(dT*weight_pop_i)/sum(weight_pop_i)   
    temp_pop_f[itime]=sum(dT*weight_pop_f)/sum(weight_pop_f)  
    temp_area[itime]=sum(dT*weight_area)/sum(weight_area)     
   
# dTgiss=rtmp[itime,:,:]
# dT = basemap.interp(dTgiss,lons_giss,lats_giss, \
#                     lons,lats,checkbounds=False,masked=False,order=1)



# Create map.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines()
# draw lat/lon grid lines every 30 degrees.
m.drawparallels(np.arange(-90.,91.,15))
m.drawmeridians(np.arange(-180.,181.,20.))
# find x,y of map projection grid.
x, y = m(lons, lats)
clevs = np.arange(-2,6,1)
cs = m.contourf(x,y,dT,clevs,cmap='gist_rainbow_r')
# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('K')
plt.title('Interpolated T Anomaly for '+str(years_giss[ntime-1]))
savefig("dT_annual-map-latest.pdf")
plt.show()

# Create arrays of latitude and longitude for dT data
lons_giss, lats_giss = np.meshgrid(lons_giss, lats_giss) 


yr_giss=np.asarray(years_giss)
plot(yr_giss,temp_mn,'k-',label='Area Av') 
plot(yr_giss,temp_pop_i,'g',label='Init Population Av')  
plot(yr_giss,temp_pop_f,'r',label='Final Population Av') 
#plot(years_giss,temp_area,'b',label='Area Av - regrid')  
xlab='Year'
ylab='Temperature Anomaly (K)'
xlabel(xlab,size=14)
ylabel(ylab,size=14)
legend(loc=2, fontsize=11)
savefig("dT_annual-timeseries-1880-2016.pdf")
figure()


jon=np.where(yr_giss >= 2000)[0]
print('Post 2000 years=',jon)

plot(yr_giss[jon],temp_mn[jon],'k-',label='Area Av') 
plot(yr_giss[jon],temp_pop_i[jon],'g',label='Init Population Av')  
plot(yr_giss[jon],temp_pop_f[jon],'r',label='Final Population Av') 
#plot(years_giss,temp_area,'b',label='Area Av - regrid')  
xlab='Year'
ylab='Temperature Anomaly (K)'
xlabel(xlab,size=14)
ylabel(ylab,size=14)
legend(loc=2, fontsize=11)
savefig("dT_annual-timeseries-2000-2016.pdf")
figure()

# Create map.
temp=tmp[0,:,:]
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines()
# draw lat/lon grid lines every 30 degrees.
m.drawparallels(np.arange(-90.,91.,15))
m.drawmeridians(np.arange(-180.,181.,20.))
# find x,y of map projection grid.
x, y = m(lons_giss, lats_giss)
clevs = np.arange(-2,6,1)
cs = m.contourf(x,y,temp,clevs,cmap='gist_rainbow_r')

# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('K')
plt.title('T Anomaly for '+str(years_giss[0]))
plt.show()

ntime=len(years)
    
plot(years,1.0e-9*pop_tot,'k*')    
xlab='Year'
ylab='Human Population (billions) '
xlabel(xlab,size=14)
ylabel(ylab,size=14)
# plot(years_giss,5*temp_mn,'r') 
savefig("Pop_annual-timeseries-1860-2005.pdf")
figure() 

# Create initial weight map.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines()
# draw lat/lon grid lines.
m.drawparallels(np.arange(-90.,91.,15))
m.drawmeridians(np.arange(-180.,181.,20.))
# find x,y of map projection grid.
x, y = m(lons, lats)
clevs = np.arange(1,100,10)
z=1.0e6*wt_pop[0,:,:]
# print('max(z) = ',max(z))
cs = m.contourf(x,y,z,clevs,cmap='gist_rainbow_r')

# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('Population weight (x$10^6$)')
plt.title('Population Weighting '+str(years[0]))
savefig("Pop_weight-map-1860.pdf")
plt.show()

# Create final weight map.
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
m.drawcoastlines()
# draw lat/lon grid lines.
m.drawparallels(np.arange(-90.,91.,15))
m.drawmeridians(np.arange(-180.,181.,20.))
# find x,y of map projection grid.
x, y = m(lons, lats)
clevs = np.arange(1,100,10)
z=1.0e6*wt_pop[ntime-1,:,:]
# print('max(z) = ',max(z))
cs = m.contourf(x,y,z,clevs,cmap='gist_rainbow_r')

# add colorbar.
cbar = m.colorbar(cs,location='bottom',pad="5%")
cbar.set_label('Population weight (x$10^6$)')
plt.title('Population Weighting '+str(years[ntime-1]))
savefig("Pop_weight-map-2005.pdf")
plt.show()

