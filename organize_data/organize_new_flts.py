#! /usr/bin/env python

"""
Organize newly obtained data into the correct file structure, sorted
into  directories for the object and filter. 

Used in the WFC3/UVIS contamination and stability monitoring program.  


Author:

	C. Shanahan, August 2016
	
Use:
	In 'set_paths.py', specify the 'new_data' directory where all the new flts are located.
	Also specify the 'data' directory in this configuration file: this is where this script will
	sort the files into the correct object/filter subdirectories
	
	If you're using this script for the standard purpose of adding new data to the exsisting
	UVIS contam analysis, you will want to set the 'temp_for_lacosmic' parameter to True.
	This will seperate the new data from the old data that has already been cleaned/analyzed and simply
	add it to the analysis. You may switch this off if you wish to redo the analysis on the entire data set. 
	
History:

	Aug. 2016:
	
		* Original script
	
"""

import sys
from set_paths import *
import glob
import os
from astropy.io import fits
import shutil
from set_paths_and_params import *


def organize_flts(path_to_flts,data_directory,photometry_temp_directory=True):

	fits_list=glob.glob(path_to_flts+'/'+'*flt.fits')
	
	stripped_filenames=[f[-18:] for f in fits_list]

	object_names=[]
	file_names=[]
	filter_names=[]
	
	#open up each fits and extract filter and object name
	for i, ifile in enumerate(fits_list):
		file_names.append(ifile)
		hdr=fits.open(ifile)[0].header
		filterr=hdr['FILTER']
		filter_names.append(filterr)
		object_name=hdr['TARGNAME']
		print filterr + ', ' + object_name
		if '70' in object_name:
			object_names.append('GRW70')
		elif '153' in object_name:
			object_names.append('GD153')
		else:
			object_names.append(object_name)
			
	#create a set of unique filter and object names
	objects=set(object_names)
	filters=set(filter_names)
	
	
	#make the directories for sorting by object/filter if they dont exist
	for object in objects:
		for filter in filters:
			new_path=data_directory+'/'+object+'/'+filter
			if not os.path.isdir(new_path):
				print 'making directory ' + new_path
   				os.makedirs(new_path)
   			if not os.path.isdir(new_path + '/flt_cleans'):
   			    print 'making directory ' + new_path + '/flt_cleans'
   			    os.makedirs(new_path + '/flt_cleans')
   			if not os.path.isdir(new_path + '/flt_masks'):
   			    print 'making directory ' + new_path + '/flt_cleans'
   			    os.makedirs(new_path + '/flt_masks')
   				
   				
	if photometry_temp_directory:
		for object in objects:
			for filter in filters:
				new_path=data_directory+'/'+object+'/'+filter+'/temp_for_photometry/'
				if not os.path.isdir(new_path):
					print 'making directory ' + new_path
					os.makedirs(new_path)
					
		for i, file in enumerate(file_names):
			object=object_names[i]
			filter=filter_names[i]
			new_temp=data_directory+'/'+object+'/'+filter+'/temp_for_photometry/'
			if not os.path.isfile(new_temp+stripped_filenames[i]):
				shutil.move(file,new_temp)
				print 'moving ' + stripped_filenames[i] + ' to temp location in ' + object + ', ' + filter 
			else: 
				print stripped_filenames[i]+' is already in the temp directory for ' + object + ', ' + filter +'. Deleting from newData'
				os.remove(file)
					
					
	else:
		for i, file in enumerate(file_names):
			print 'moving ' + stripped_filenames[i]
			object=object_names[i]
			filter=filter_names[i]
			shutil.move(file,data_directory+'/'+object+'/'+filter+'/')


	
					
		
if __name__ == '__main__':

	paths = set_paths_and_params()
	path_to_flts = paths['data']+'/new_flts'
	data_directory = paths['data']+'/data'
	organize_flts(path_to_flts,data_directory)
	
