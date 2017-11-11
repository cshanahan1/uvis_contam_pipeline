from make_uvis_contam_plot import *
from set_paths_and_params import *
import glob
import os.path
import os

def main_wrapper_make_uvis_contam_plots(paths):
	
	data_dir = paths['data']
	filters = paths['filters']
	objects = paths['objects']
	plot_output_dir = paths['plot_output_dir']+'/'
	
	#create stats.dat file
	with open(plot_output_dir+'/stats.dat','w') as f:
		f.write('#obj,amp,filter,slope_mjd,b,stdev_mjd,slope_yr,stdev_yr\n')
		
	for obj in objects:
		print obj
		paths = [data_dir+'/data/{0}/{1}/flt_cleans'.format(obj, f) for f in filters]
		for path in paths:
			ifile = glob.glob(path+'/*.dat')
			if len(ifile)>0:
			    print ifile[0]
			    filename =  os.path.basename(ifile[0])
			    #filtername = filename.replace('_photcat.dat','')
			    filtername = filename.replace('_photcat.dat','')

			    #make plots
			    main_fluxDif_MJD(ifile[0], filtername, obj,plot_output_dir,wrapper=True)
			else:
				print 'no photcat. moving on' 
			
if __name__ == '__main__':

	main_wrapper_make_uvis_contam_plots(set_paths_and_params())