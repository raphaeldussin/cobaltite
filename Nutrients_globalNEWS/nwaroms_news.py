import lib_Global_NEWS_to_ROMS as lgn

# pick a ROMS domain
grdname='NWA'
csvdir = '/Users/raphael/WORK/GlobalNEWS/export_csv/'

# setup domain
nwanews = lgn.news2roms(grdname)
# read database from csv files
nwanews.read_database(csvdir)
# extract rivers on domain
nwanews.extract_domain()
# check plot
#nwanews.plot_river_mouth(db='domain')
# move rivers to closest ocean point
nwanews.move_rivermouth2roms()
# expand to coastal area and write nutrient concentration to file
nwanews.create_all_nutrients('./NEWS_nutrients_NWA_RD-04-10-2017.nc')

