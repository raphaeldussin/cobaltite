from brushcutter import *

#cfile = '/Volumes/ESM-ADMIN/data-nico/PSD_PFT_C_MO_climatology_SeaWiFS_06_test.nc'
cfile = '/Volumes/ESM-ADMIN/data-nico/out.nc'
romsgrd = '/Users/raphael/STORAGE/ROMS/GRIDS/CCS_7k_0-360_fred_grd.nc'

# ---------- define segments on MOM grid -----------------------
domain = obc_segment('domain', romsgrd,istart=0,iend=181,jstart=0,  jend=481,target_model='ROMS')

# ---------- define variables on each segment ------------------
C_biomass_domain = obc_variable(domain,'C_biomass_total',geometry='line',obctype='radiation')
#C_biomass2_domain = obc_variable(domain,'C_biomass_truc',geometry='line',obctype='radiation')

# ---------- interpolate T/S from SODA file -------------------
interp_t2s = C_biomass_domain.interpolate_from(cfile,'C_biomass_total',frame=0,depthname='',method='bilinear')
#C_biomass2_domain.interpolate_from(cfile,'C_biomass_truc',frame=0,depthname='',method='bilinear',interpolator=interp_t2s)

# ---------- list segments and variables to be written -------
list_segments = [domain]

list_variables = [C_biomass_domain] 

list_vectvariables = []

#----------- time --------------------------------------------
time = C_biomass_domain.timesrc

# ---------- write to file -----------------------------------
lib_ioncdf.write_obc_file(list_segments,list_variables,list_vectvariables,time,output='test_CCS.nc')
