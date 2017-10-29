def pop():
    import numpy as np
    from pylab import log, arctan, sum
    from netCDF4 import Dataset

## Change file pathway directory "dir" below for local usage ##
    dir=''
    file_in=Dataset(dir+'histsoc_population_0.5deg_1861-2005.nc4')

# Read in the data from the NetCDF files
    lats=file_in.variables['lat'][:].astype('float64')
    lons=file_in.variables['lon'][:].astype('float64')
    times=file_in.variables['time'][:].astype('float64')
    rpop=file_in.variables['number_of_people'][:,:,:].astype('float64')    # Time, lat, lon

# Set-up timseries variables
    ntime=rpop.shape[0]
    print('ntime = ',ntime)
    pop_tot=np.zeros(ntime)
    years=[0]*ntime
    
# Create arrays of latitude and longitude
#    lons, lats = np.meshgrid(lons, lats)

# Weights for population averages
    wt_pop=rpop

# Calculate and plot global mean anomalies
    for itime in range(0, ntime):

        years[itime]=1861+itime
        pop=rpop[itime,:,:]
        pop_tot[itime]=sum(pop)          
        wt_pop[itime,:,:]=pop/pop_tot[itime]
        
    return years,lats,lons,pop_tot,wt_pop;

    
