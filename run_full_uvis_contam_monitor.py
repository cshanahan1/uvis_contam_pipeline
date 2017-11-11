#! /usr/bin/env python

"""Runs all components of the  UVIS contamination monitor code. Steps include, in order:

1. Querys the QuickLook database for files from all contam proposals by proposal ID,
	checks them against files that have already been analyzed or retrieved. 
	If the files are new, it also checks the CALWF3 version. If the pipeline version 
	matches the desired version (specified in set_paths_and_params.py), the FLTs are 
	copied directly into the new data directory. If the version is not correct, then the 
	root names of the FLTs are written to a file for manual retrieval from MAST, where 
	files are processed with the most up to date pipeline/reference files.

2. 	Runs organize_data.organize_new_flts.

	Organizes the new data based on object / filter, puts the new _flt.fits in a 
	'temp_for_photometry' subdirectory within each object / fillter directory to keep new, 
	unprocessed data seperate from data that has already been analyzed. 
	
3.  Runs cosmic_ray_removal.run_lacosmic_on_filters

	Runs the cosmic ray rejection wrapper on the all new data, which has been organized by 
	object and by filter into temporary directories for processing. The CR rejection 
	returns cleaned _flt.fits (_flt.clean.fits) and the associated mask files for each
	observation. The masks are kept for manual inspection in case something goes 
	awry, but are not essential to the analysis pipeline. 
	
4. Performs point source photometry on the new, cleaned FLTs. See docstrings in 
   ptsrc_photometry.py for details on the process. Photometry catalogs for each 
   object/filter are appended with fluxes/exposure times / etc for the newly 
   processed data. 
   
5. Temporary directories are cleaned out, with the fully processed data being moved to
	the permenant data storage directory. See organize_data.clean_temporary_directories.
   
5. Creates plots, and appends newly processed data to the stats catalogs which have slopes
	and standard deviations for the full data set.
	

Authors
-------
	C. Shanahan, October 2016 
	

Use
-------
	>> python run_full_uvis_contam_monitor.py
	
	
Outputs
-------
	If this is the first time you are running this with all the unprocessed FLTs from
	contam monitor programs, this will create all the cleaned FLTs and the masks, 
	photometry catalogs, statistics tables, and plots for all files.
	
	If this is run periodically, it will create the cleaned FLTs and the masks, 
	photometry catalogs, statistics tables, and plots for all NEW files that do not exist
	in the analysis already.
	
	
Notes
-------
	* Make sure set_paths_and_params.py has been updated with the correct paths to the 
	local data directories etc. This script runs in one command, no args, from the root 
	directory of this code base. All code depends on the parameters set in this file.
	  
	* While this monitor is nearly totally automated by running this script, there are a 
	   few things that will need to be manually updated every once in a while.
	   
	   	1. the 'programs' list in this script under all the imports. These are all the 
	   		proposal IDs that have contam monitor data, so the most recent one will need
	   		to be added every cycle. 
	   		
	   	2. If new files were not processed with the correct specified version of the 
	   		calibration pipeline, they will need to be manually requested from MAST. 
	   		

"""
import time 
import os
import numpy as np
import shutil
import sys

import logging
from pyql.logging.logging_functions import configure_logging
from pyql.logging.logging_functions import log_info
from pyql.logging.logging_functions import log_fail

from set_paths_and_params import *
from organize_data.get_and_check_new_data import *
from organize_data.organize_new_flts import *
from organize_data.run_calwf3 import *
from organize_data.clean_temporary_directories import *
from cosmic_ray_removal.run_lacosmic_on_filters import *
from ptsrc_photometry.wrapper_ptsrc_photometry_photutils import *
from plotting_and_statistics.wrapper_make_uvis_contam_plots import *




def object_names():

	"""
	Contains a list of every possible iteration of how the object names might be spelled
	add to this list if you find another!
	
	Returns: 
		list of object names that the database recognizes for GD153 and GRW70
	
	"""
	return ['GRW+70D5824','GD-153','GD153']
		
#--------------------------------------------------------------------------------#

		
def main_uvis_contam_monitor(paths,calibrate_raws = True):

	#######################################################
	#		  query the QL database for new data          #
	#######################################################

	#unpack paths & parameters - from set_paths_and_params.py
	data_dir = paths['data']+'/data'
	new_data_dir = paths['data']+'/new_flts'
	raw_dir = paths['data']+'/new_raws'
	filters=paths['filters']
	objects=paths['objects']
	proposal_ids = paths['prop_ids']
	cal_vers = paths['cal_ver']


	#Iterate over proposal IDS to search for unprocessed files in the QL database
	
	print 'Querying the QL database for FLTs from contam monitor programs.'
	
	for proposal in proposal_ids:
	
		flts=query_ql(proposal,object_names())
		print len(flts), ' total observations in proposal ', proposal
		
		if len(flts)==0:
			continue
			
		#check all files in QL directories for proposal against files in data directory
		
		new_flts=check_new_files_against_existing_files(flts,data_dir,new_data_dir,raw_dir)
		
		
		if len(new_flts)==0:
			print 'no new files found for proposal {0}'.format(proposal)
			continue
			
		else:
			print 'found '+str(len(new_flts))+' new files for proposal {0}'.format(proposal)
		
			cal_version_flts=check_cal_version(new_flts, cal_vers)
			
			new_flts=cal_version_flts[0] #files calibrated with correct version of calwf3
			
			print len(new_flts) , ' of the new files were calibrated with the correct version of the pipeline. adding files to new data directory.'
				
			for new_flt in new_flts:
				shutil.copy(new_flt, new_data_dir)
				
			wrong_version_flts=cal_version_flts[1] 
			wrong_version_raws=[f.replace('flt','raw') for f in wrong_version_flts] 
			
			print len(wrong_version_flts) , ' of the new flts found were not calibrated with the correct version of calwf3 specified in set_paths_and_params.py. Copying the raw files to new_raws for calibration.'

			for raw in wrong_version_raws:
				shutil.copy(raw, raw_dir)
			
	new_flts=glob.glob(new_data_dir+'/*flt.fits')[0:10]
	new_raws = glob.glob(raw_dir+'/*flt.fits')
	
	#######################################################
	#		  calibrate the raw files if there are any    #
	#######################################################
	
	"""if calibrate_raws:
	    print 'Calibrating RAW files.' 
	    cwd = os.getcwd()
	    os.chdir(raw_dir)
	    print os.getcwd()
	    os.system('python ' + 'run_calwf3.py')
	    os.chdir(cwd)"""

	#######################################################
	#		  now organize the newly fetched data         #
	#######################################################
	
	if len(new_flts) > 0:
		#put into temp holding areas for lacosmic / photometry
		print 'organizing new FLTs into proper star and filter subdirectories.'
		organize_flts(new_data_dir,data_dir,photometry_temp_directory=True)
	if len(new_raws) > 0:
		print 'organizing new FLTs into proper star and filter subdirectories.'
		organize_flts(raw_dir,data_dir,photometry_temp_directory=True)
		
		print "done organizing new FLTs."
	else:
		print 'no new FLTs to organize.'
	

	#######################################################
	#		  Run LAcosmic on the new data                #
	#######################################################
	
	print 'Cleaning cosmic rays from new images.'
	main_run_lacosmic_on_filters()
	print 'Done cleaning cosmic rays from new images.'
	
	#######################################################
	#		         Run the photometry                   #
	#######################################################
	print 'Beginning photometry on new images.'
	wrapper_ptsrc_photometry_photutils(set_paths_and_params())
	#main_photom()
	print 'Photometry complete.'
	
	#######################################################
	#		         Clean up temp directories            #
	#######################################################
	
	print 'Cleaning up temporary directories.'
	
	clean_temp_directories(paths)
	
	
	#######################################################
	#		         Generate plots                       #
	#######################################################
	#print 'Generating plots and stats file...'
	main_wrapper_make_uvis_contam_plots(paths)
	
	
	#######################################################
	#		         Make stats catalogs                  #
	#######################################################
	

		
	
if __name__=='__main__':

	paths=set_paths_and_params()	
	main_uvis_contam_monitor(paths)	



		
			
	 
	
		
		
		


		
		
		