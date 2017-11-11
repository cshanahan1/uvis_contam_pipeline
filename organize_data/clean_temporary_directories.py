"""Cleans up the temporary holding directories used to keep new files seperate for processing. 
Moves all the files to the permenant storage location and deletes the temp_for_lacosmic directories.

Used in the WFC3/UVIS contamination and stability monitoring program.  

Author:

	C. Shanahan, August 2016
	
Use:
	Set the filters you wish to run this script over, and the appropriate location of the data
	directory that contains the object in set_paths (or manually here..).
	
	
History:

	Aug. 2016:
	
		*Original script
	
"""



import sys
from set_paths_and_params import *
import glob
from os.path import isdir
import os
from astropy.io import fits
import shutil

def clean_temp_directories(paths):

	objects = paths['objects']
	filters = paths['filters']
	data_dir = paths['data']+'/data'
	

	for object in objects:
	
		object_dir=data_dir+'/'+object+'/'
		flt_dirs=[object_dir+f for f in filters]
	
		for dir in flt_dirs:
			#first move all the data files out of the temp directory
			temp_dir=dir+'/temp_for_photometry/'
			flts=glob.glob(temp_dir+'*flt.fits')
			for flt in flts:
				print 'moving ' + flt + ' to ' + dir
				shutil.move(flt,dir)
			
			###now clean out the /temp_for_photometry/flt_cleans
			stuff_in_flt_cleans=glob.glob(temp_dir+'flt_cleans/*.fits')

					
			#make a permanent flt_cleans in the filter directory if it doesn't already exist
			if not isdir(str(dir)+'/flt_cleans'):
				os.makedirs(dir+'/flt_cleans')
				print 'making directory ' + dir+'/flt_cleans'
			
			for thing in stuff_in_flt_cleans:
				print 'moving ' + thing + ' to ' + dir+'/flt_cleans'
				if '.dat' not in thing: #don't overwrite photcat just in case !
					shutil.move(thing,dir+'/flt_cleans')
				else:
					os.remove(thing)
					print 'removing temp catalog in temp /flt_cleans'

		
			###same process for flt_masks
			stuff_in_flt_masks=glob.glob(temp_dir+'flt_masks/*')
		
			if not isdir(str(dir)+'/flt_masks'):
				os.makedirs(dir+'/flt_masks')
				print 'making directory ' + dir+'/flt_masks'
		
			for thing in stuff_in_flt_masks:
				print 'moving ' + thing + ' to ' + dir+'/flt_masks'
				shutil.move(thing,dir+'/flt_masks')
			

		
if __name__ == '__main__':

	clean_temp_directories(set_paths_and_params())
		
		
		



