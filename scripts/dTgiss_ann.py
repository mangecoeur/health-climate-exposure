def dTgiss_ann():
    import numpy as np
    from pylab import sum
    from netCDF4 import Dataset
    
    ## Change file pathway directory "dir" below for local usage ##
    dir=''
    file_in=Dataset(dir+'gistemp1200_ERSSTv4.nc')

# Read in the data from the NetCDF files
    lats=file_in.variables['lat'][:].astype('float64')
    lons=file_in.variables['lon'][:].astype('float64')
    rtmp=file_in.variables['tempanomaly'][:,:,:].astype('float64')    # Time, lat, lon


# Set-up timseries variables
    ntime=rtmp.shape[0] 
    nyears=int(ntime/12)
    print('nyears = ',nyears)
    years=[0]*nyears
    temp_mn=np.zeros(nyears)
    tmp=np.zeros((nyears,rtmp.shape[1],rtmp.shape[2]))

# Create arrays of latitude and longitude
    lons_2d, lats_2d = np.meshgrid(lons, lats)

# Weighting factor for global means
    weight=np.cos(lats_2d*np.pi/180.0)

# Calculate and plot global mean anomalies
    for n in range(0, nyears):
        years[n] = 1880 + n
        i1=n*12
#       tmp[n,:,:]=rtmp[i1,:,:]/12.0
        for m in range(0,12):
            i2=i1+m
            tmp[n,:,:]=tmp[n,:,:]+rtmp[i2,:,:]/12.0 
            
        temp=tmp[n,:,:]
        temp_mn[n]=sum(temp*weight)/sum(weight)  
        print(years[n],temp_mn[n])
        
    print('Shape tmp = ',tmp.shape,', Shape weight = ',weight.shape)    
    return years,lats,lons,temp_mn,tmp;

    
