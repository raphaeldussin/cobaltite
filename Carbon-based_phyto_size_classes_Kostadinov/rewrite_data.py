from lib_rewrite_data import rewrite_data

dirdatain=''
dirdataout=''

june = rewrite_data()
june.rewrite_obs_data(fileobs=dirdatain + 'PSD_PFT_C_MO_climatology_SeaWiFS_06.nc',fileout=dirdataout + 'PSD_PFT_C_MO_climatology_SeaWiFS_06.nc')
