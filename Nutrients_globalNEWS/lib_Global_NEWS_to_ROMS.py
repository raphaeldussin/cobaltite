#!/usr/bin/env python

# import libs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pyroms
import numpy as np
import netCDF4 as nc
import matplotlib.cm as cm
import matplotlib.colors as cl
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
import pandas as pd
import scipy.ndimage as si

class news2roms():

	def __init__(self,grdname):
		self.grd = pyroms.grid.get_ROMS_grid(grdname)
		self.lonmin = self.grd.hgrid.lon_rho.min() - 360.
        	self.lonmax = self.grd.hgrid.lon_rho.max() - 360.
        	self.latmin = self.grd.hgrid.lat_rho.min()
        	self.latmax = self.grd.hgrid.lat_rho.max()
		return None

	def setup_map(self):
	        ''' set up the map '''
		plt.figure(figsize=[8.,8.])
		ax = plt.axes(projection=ccrs.PlateCarree())
		ax.set_extent([self.lonmin, self.lonmax, self.latmin, self.latmax])
		ax.coastlines()
		return ax

	def from_river_mouth_to_roms_cell(self,lon,lat):
		'''function to find closest ROMS ocean land point to river mouth'''
		# Convert latitude and longitude to
		# spherical coordinates in radians.
		degrees_to_radians = np.pi/180.0

		# phi = 90 - latitude
		phi1 = (90.0 - lat)*degrees_to_radians
		phi2 = (90.0 - self.grd.hgrid.lat_rho)*degrees_to_radians

		# theta = longitude
		theta1 = lon*degrees_to_radians
		theta2 = self.grd.hgrid.lon_rho*degrees_to_radians

		# Compute spherical distance from spherical coordinates.

		# For two locations in spherical coordinates
		# (1, theta, phi) and (1, theta, phi)
		# cosine( arc length ) =
		#    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
		# distance = rho * arc length

		cos = (np.sin(phi1)*np.sin(phi2)*np.cos(theta1 - theta2) +
		       np.cos(phi1)*np.cos(phi2))
		arc = np.arccos( cos )

		# Remember to multiply arc by the radius of the earth
		# in your favorite set of units to get length.
		arc[np.where(self.grd.hgrid.mask_rho == 0)] = 1.e36
            
		jcell,icell = np.unravel_index(arc.argmin(),self.grd.hgrid.lon_rho.shape)
		jcell=int(jcell)
		icell=int(icell)
		if self.grd.hgrid.mask_rho[jcell,icell] == 0:
			print 'Error : cell on land'
		if arc.min() > 0.01: # exceeded proximity threshold 0.01 * earth radius = 64km
			jcell=0 ; icell=0
		return jcell, icell

	def read_database(self,csvdir):
		''' read csv database into pandas '''
		# pandas read excel function not working, I needed to export to individual csv files
		# read the NEWS data from csv file
		news_basins       = csvdir + 'basins-Table1.csv'
		news_riverexports = csvdir + 'river_exports-Table1.csv'
		news_hydro        = csvdir + 'hydrology-Table1.csv'

		df_news_basins = pd.read_csv(news_basins)
		df_news_exports = pd.read_csv(news_riverexports)
		df_news_hydro = pd.read_csv(news_hydro)

		self.newsdb = pd.concat([df_news_basins,df_news_exports,df_news_hydro],axis=1)

		# write discharge in m3/s from Qact in km3/yr
		self.newsdb['Discharge'] =self.newsdb['Qact'] * 1e+9 / 86400 / 365.25

		# Raw form of nutrient loads is Mg N yr-1, convert to moles N sec-1
		self.newsdb['DIN_conc'] = self.newsdb['Ld_DIN'] /86400/365*1e6/14 / self.newsdb['Discharge']
		self.newsdb['DIP_conc'] = self.newsdb['Ld_DIP'] /86400/365*1e6/31 / self.newsdb['Discharge']
		self.newsdb['DON_conc'] = self.newsdb['Ld_DON'] /86400/365*1e6/14 / self.newsdb['Discharge']
		self.newsdb['DOP_conc'] = self.newsdb['Ld_DOP'] /86400/365*1e6/31 / self.newsdb['Discharge']
		self.newsdb['DOC_conc'] = self.newsdb['Ld_DOC'] /86400/365*1e6/12 / self.newsdb['Discharge']
		self.newsdb['POC_conc'] = self.newsdb['Ld_POC'] /86400/365*1e6/12 / self.newsdb['Discharge']
		self.newsdb['TSS_conc'] = self.newsdb['Ld_TSS'] /86400/365*1e6    / self.newsdb['Discharge']
		
		self.newsdb['NO3_CONC']   =        self.newsdb['DIN_conc']
		self.newsdb['LDON_CONC']  = 0.3 *  self.newsdb['DON_conc']
		self.newsdb['SLDON_CONC'] = 0.35 * self.newsdb['DON_conc']
		self.newsdb['SRDON_CONC'] = 0.35 * self.newsdb['DON_conc']
		self.newsdb['PO4_CONC']   =        self.newsdb['DIP_conc']
		self.newsdb['LDOP_CONC']  = 0.3 *  self.newsdb['DOP_conc']
		self.newsdb['SLDOP_CONC'] = 0.35 * self.newsdb['DOP_conc']
		self.newsdb['SRDOP_CONC'] = 0.35 * self.newsdb['DOP_conc']
		self.newsdb['NDET_CONC']  = 0. *   self.newsdb['DIN_conc']
		self.newsdb['PDET_CONC']  = 0. *   self.newsdb['DIP_conc']

		# expand dataframe
		nmouth = self.newsdb['mouth_lon'].shape[0]
		jcell_roms = np.zeros(nmouth,dtype='i4')
		icell_roms = np.zeros(nmouth,dtype='i4')
		jcellpos = pd.DataFrame(jcell_roms, columns=['jcell_roms'])
		icellpos = pd.DataFrame(icell_roms, columns=['icell_roms'])
		self.newsdb = pd.concat([self.newsdb,jcellpos,icellpos],axis=1)
		return None

	def extract_domain(self):
		''' extract database on roms domain '''
		self.newsdb_domain = self.newsdb[(self.newsdb.mouth_lat >= self.latmin)& (self.newsdb.mouth_lat <= self.latmax) &
		(self.newsdb.mouth_lon >= self.lonmin) & (self.newsdb.mouth_lon <= self.lonmax)]

		# we want only rivers whose discharge > 10 m3/s
		self.newsdb_domain = self.newsdb_domain[(self.newsdb_domain.Discharge > 10)]

		self.newsdb_domain = self.newsdb_domain.sort_values('Discharge', ascending=True)
		return None

	def plot_river_mouth(self,db='global'):
		''' plots '''
		ax = self.setup_map()
		if db == 'global':
			ax.plot(self.newsdb['mouth_lon'],self.newsdb['mouth_lat'],'ko')
		elif db == 'domain':
			ax.plot(self.newsdb_domain['mouth_lon'],self.newsdb_domain['mouth_lat'],'ko')
		plt.show()
		return None

	def move_rivermouth2roms(self):
		''' move river mouth to closest ocean point '''
		nmouth_local = self.newsdb_domain['mouth_lon'].shape[0]
		idj = self.newsdb_domain.columns.get_loc("jcell_roms")
		idi = self.newsdb_domain.columns.get_loc("icell_roms")

		for km in np.arange(nmouth_local):
			jcell_roms, icell_roms = \
			self.from_river_mouth_to_roms_cell(self.newsdb_domain['mouth_lon'].values[km], self.newsdb_domain['mouth_lat'].values[km])

			dbindex = self.newsdb_domain.index[km]
			self.newsdb_domain.loc[dbindex,"jcell_roms"] = jcell_roms
			self.newsdb_domain.loc[dbindex,"icell_roms"] = icell_roms
		return None

	def create_rivers_input(self,field,debug=False):

		rivers_input = np.zeros(self.grd.hgrid.mask_rho.shape)
		nmouth_local = self.newsdb_domain['mouth_lon'].shape[0]

		for kriver in np.arange(nmouth_local):
			print ' working on ' , self.newsdb_domain['basinname'].values[kriver]
			#if self.newsdb_domain['basinname'].values[kriver] == 'Delaware':
			#	self.create_one_river(rivers_input,kriver,field,rspread=10,debug=True)
			#else:
			#	self.create_one_river(rivers_input,kriver,field,rspread=10)
			self.create_one_river(rivers_input,kriver,field,rspread=10)

		rivers_input = rivers_input * self.grd.hgrid.mask_rho

		if debug:
			ax = self.setup_map()
			rivers_input = np.ma.masked_values(rivers_input,0.)
			ax.pcolormesh(self.grd.hgrid.lon_rho-360,self.grd.hgrid.lat_rho,self.grd.hgrid.mask_rho,cmap=cm.binary_r,vmin=-3,vmax=1)
			C = ax.pcolormesh(self.grd.hgrid.lon_rho-360,self.grd.hgrid.lat_rho,rivers_input)
			plt.colorbar(C)
			plt.show()
		return rivers_input

	def create_one_river(self,datainout,kriver,field,rspread=10,debug=False):
		''' create nutrient inputs for one river '''
		jcenter = self.newsdb_domain['jcell_roms'].values[kriver]
		icenter = self.newsdb_domain['icell_roms'].values[kriver]
		nutrient_value = self.newsdb_domain[field].values[kriver]

		if ( jcenter == 0 ) and ( icenter == 0 ):
			print 'out of domain, nothing to do'
			pass
		else:
			mask_river = np.zeros(self.grd.hgrid.mask_rho.shape)
			mask_river[jcenter,icenter] = 1 

			imin   = icenter - rspread
			jmin   = jcenter - rspread 
			imax   = icenter + rspread + 1
			jmax   = jcenter + rspread + 1

			mask_zoom_in_before = mask_river[jmin:jmax,imin:imax]
			lsm_zoom_in  = self.grd.hgrid.mask_rho[jmin:jmax,imin:imax]

			mask_zoom_in_after = self.mask_river(mask_zoom_in_before,lsm_zoom_in)
	
			if debug:
				plt.figure()
				plt.pcolormesh(mask_zoom_in)
				plt.title('before')
				plt.figure()
				plt.pcolormesh(lsm_zoom_in)
				plt.title('mask')
				plt.figure()
				plt.pcolormesh(mask_zoom_in)
				plt.title('after')
				plt.show()

			tmp = datainout[jmin:jmax,imin:imax].copy()
			tmp[np.where(mask_zoom_in_after == 1)] = nutrient_value
			datainout[jmin:jmax,imin:imax] = tmp.copy()
		return datainout

	def mask_river(self,data,mask):
		''' generate mask for one river '''
		data_old = data.copy()
		for kk in np.arange(1000):
			data_new = si.binary_dilation(data_old)
			data_new = data_new * mask
			if data_old.mean() == data_new.mean():
				print 'converged at iteration ', kk
				break
			data_old = data_new.copy()
		return data_new

	def create_all_nutrients(self,fileout):
		''' loop over all nutrients '''
		#list_nutrients = ['NO3_CONC','LDON_CONC','SLDON_CONC','SRDON_CONC','PO4_CONC','LDOP_CONC','SLDOP_CONC','SRDOP_CONC','NDET_CONC','PDET_CONC']

		fid = nc.Dataset(fileout,'w',format='NETCDF3_64BIT_OFFSET')
		fid.description = 'Global NEWS to ROMS grid (raphael.dussin@gmail.com)'
	        # dimensions
		ny, nx = self.grd.hgrid.mask_rho.shape
	        fid.createDimension('eta_rho', ny)
	        fid.createDimension('xi_rho', nx)
	        fid.createDimension('time', None)
	        # variables
	        latitudes  = fid.createVariable('lat_rho',    'f4', ('eta_rho','xi_rho',))
	        longitudes = fid.createVariable('lon_rho',    'f4', ('eta_rho','xi_rho',))
	        times      = fid.createVariable('time',       'f4', ('time',))
        	no3_conc   = fid.createVariable('NO3',   'f8', ('time','eta_rho','xi_rho',))
        	ldon_conc  = fid.createVariable('LDON',  'f8', ('time','eta_rho','xi_rho',))
        	sldon_conc = fid.createVariable('SLDON', 'f8', ('time','eta_rho','xi_rho',))
        	srdon_conc = fid.createVariable('SRDON', 'f8', ('time','eta_rho','xi_rho',))
        	po4_conc   = fid.createVariable('PO4',   'f8', ('time','eta_rho','xi_rho',))
        	ldop_conc  = fid.createVariable('LDOP',  'f8', ('time','eta_rho','xi_rho',))
        	sldop_conc = fid.createVariable('SLDOP', 'f8', ('time','eta_rho','xi_rho',))
        	srdop_conc = fid.createVariable('SRDOP', 'f8', ('time','eta_rho','xi_rho',))
        	ndet_conc  = fid.createVariable('NDET',  'f8', ('time','eta_rho','xi_rho',))
        	pdet_conc  = fid.createVariable('PDET',  'f8', ('time','eta_rho','xi_rho',))
	        # data
	        latitudes[:]      = self.grd.hgrid.lat_rho.copy()
	        longitudes[:]     = self.grd.hgrid.lon_rho.copy()
	        times[:]          = 0.
	        no3_conc[0,:,:]   = self.create_rivers_input('NO3_CONC')
	        ldon_conc[0,:,:]  = self.create_rivers_input('LDON_CONC')
	        sldon_conc[0,:,:] = self.create_rivers_input('SLDON_CONC')
	        srdon_conc[0,:,:] = self.create_rivers_input('SRDON_CONC')
	        po4_conc[0,:,:]   = self.create_rivers_input('PO4_CONC')
	        ldop_conc[0,:,:]  = self.create_rivers_input('LDOP_CONC')
	        sldop_conc[0,:,:] = self.create_rivers_input('SLDOP_CONC')
	        srdop_conc[0,:,:] = self.create_rivers_input('SRDOP_CONC')
	        ndet_conc[0,:,:]  = self.create_rivers_input('NDET_CONC')
	        pdet_conc[0,:,:]  = self.create_rivers_input('PDET_CONC')
	        # close
	        fid.close()
		return None




