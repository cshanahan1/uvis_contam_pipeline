#! /usr/bin/env python


""" 
This script is a wrapper that runs the python routine LAcosmic - imported from 'lacosmic.py'
-on all filters (set in set_paths.py) and objects specified. A file sorting function is also
called from 'lacosmic.py' to sort all the clean/mask files for inspection/photometry.

Usage:

	After specifying paths/etc in 'set_paths.py', running this script will perform cosmic ray
	removal on all 'flt' files for each object/filter specified.
	
	Parameters to set in '__main__':
	
		objects: list of objects you wish to run this over
		create_pngs: If True, diagnostic pngs will be created of clean/masks. slows down 
					 excecution, but important for tweaking parameters in the param dictionary
		temp_for_photometry: If true, only .cleans from the 'temp_for_photometry' directory - presumably
						   the files that haven't been processed. If not, the .cleans from the main data
						   directory's 'flt_cleans' subdirectory will be processed.
		
		These parameters must be set as lists equal in length to the number of objects. This is so multiprocessing's 
		'Pool' can use the zipped arguements.
		
	
	See the documentation in the main script where this wrapper imports functions from - 'lacosmic.py' - 
	for more information.
		
		
	Author:
	
	 C Shanahan, Jul. 2016
	
"""

import glob
from lacosmic import *
from set_paths_and_params import *
from lacosmic_param_dictionary import lacosmic_param_dictionary 
from astropy.io import fits as pyfits
from multiprocessing import Pool

def lacosmic_on_filters(args):
	obj_list=set_paths_and_params()['objects']
	filter=args[0]
	create_pngs=args[1]
	temp_folder=args[2]
	data_dir=set_paths_and_params()['data']+'/data'
	#filters=glob.glob(data_dir+'/'+object+'/*')
	#filters_list = [path.replace(data_dir+'/'+object+'/','') for path in filters]
	dirs_list=[data_dir+'/'+obj+'/'+filter for obj in obj_list]
	if temp_folder:
		dirs_list=[d + '/temp_for_photometry' for d in dirs_list]
	param_dict = lacosmic_param_dictionary()
	for dir in dirs_list:
		if 'GRW70' in dir:
			object = 'GRW70'
		if 'GD153' in dir:
			object = 'GD153'
		fits_list = glob.glob(str(dir)+'/'+'*flt.fits')
		if len(fits_list) == 0 : 
			print 'no temp directory for / empty temp directory' + dir
			continue 
		print 'entering filter directory ' + dir
	
		if filter not in param_dict.keys():
			print "Filter {0} not in Param Dictionary. Using default values.".format(filter)
			sigclip, objlim  = 5.0,2
			sigfrac,niter,sigclip = 0.3, 3, 9.5
			sigclip_pf=9.5
		else:
			sigclip = param_dict[filter][0]
			sigfrac = param_dict[filter][1]
			objlim = param_dict[filter][2]
			niter = param_dict[filter][3]
			sigclip_pf = param_dict[filter][4]
		
		
		for fits in fits_list:
		
			if temp_folder:
				origin=data_dir+'/'+str(object)+'/'+str(filter)+'/temp_for_photometry/'
				cleanfile = fits.replace('_flt.fits','_flt.clean.fits')
				folder_cleanfile=cleanfile.replace('/temp_for_photometry/','/temp_for_photometry/flt_cleans/')
				
			
				if os.path.isfile(cleanfile):
					print cleanfile + ' this already exists in working dir! moving on'
					continue
				if os.path.isfile(folder_cleanfile):
					print folder_cleanfile + ' this already exists in flt_cleans, moving on' 
					continue
				
			else:
				origin=data_dir+'/'+str(object)+'/'+str(filter)
				
			print 'working on ' , fits
			run_lacosmic(fits, sigclip, sigfrac, objlim, niter, sigclip_pf)
			
		#now, sort the files!	
		sort_files(origin,temp_folder,create_pngs)
			

def main_run_lacosmic_on_filters():
	data_dir=set_paths_and_params()['data']
	#make sure to set the correct object!
	objects = set_paths_and_params()['objects']
	filters_list=set_paths_and_params()['filters']
	create_pngs=[False]*len(filters_list)
	temp_for_photometry=[True]*len(filters_list)
	
	args=zip(filters_list,create_pngs,temp_for_photometry)
	
	p = Pool(8)
	
	p.map(lacosmic_on_filters,args)
	

		
	for filt in filters_list:
		lacosmic_on_filters((filt,False,True))
	
if __name__ == '__main__':
	main_run_lacosmic_on_filters()

	
	




