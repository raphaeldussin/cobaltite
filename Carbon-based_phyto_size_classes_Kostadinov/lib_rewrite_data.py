import os
import os.path
import numpy as np
import netCDF4 as nc

class rewrite_data():

	def __init__(self):
		return None

	def rewrite_obs_data(self,fileobs='',timevalue=0,fileout='./out.nc'):
		fid  = nc.Dataset(fileobs)
		var1 = fid.variables['C_biomass_total'][:]		
		lat  = fid.variables['Grid_Latitudes'][:]		
		lon  = fid.variables['Grid_Longitudes'][:]		
		nx,ny = var1.shape
		spval = -9999.

		var1out = var1.transpose()
		var1out[np.where(np.isnan(var1out))] = spval

		var1out = var1out[::-1,:]
		lat = lat[::-1]

		fout = nc.Dataset(fileout,'w',format='NETCDF3_64BIT_OFFSET')
		fout.description = ''
		#fout.whatever # global attribute
		# dimensions
		fout.createDimension('time', None)
		fout.createDimension('lat', ny)
		fout.createDimension('lon', nx)
		# variables
		time       = fout.createVariable('time', 'f8', ('time',))
		latitudes  = fout.createVariable('lat', 'f8', ('lat',))
		longitudes = fout.createVariable('lon', 'f8', ('lon',))
		outdata1   = fout.createVariable('C_biomass_total', 'f8', ('time','lat','lon',),fill_value=spval)

		outdata1.units = 'mg/m^3'
		outdata1.Size_range = '0.5-50 um (micrometers) in diameter'

		time.units = 'month'
		time.calendar = '1-12'

		# data
		time[0] = timevalue
		latitudes[:]    = lat
		longitudes[:]   = lon
		outdata1[0,:,:] = var1out
		# close
		fout.close()
