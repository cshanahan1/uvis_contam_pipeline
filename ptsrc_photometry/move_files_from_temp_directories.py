"""
	Organizes the newly ingested files in 'temp_for_lacosmic' directories back to the
	main date storage directories after they have been processed.
	
	Usage:
		Imported by the ptsrc_photometry code to organize files.
		
"""

from set_paths import *
import shutil
import glob

from os.path import isdir

import os

data_dir=set_paths()['data']

def move_files_from_temp_directories(objects,filters):

	for object in objects:
	
		object_dir=data_dir+'/'+object+'/'
		flt_dirs=[object_dir+f for f in filters]
	
		for dir in flt_dirs:
			#first move all the data files out of the temp directory
			temp_dir=dir+'/temp_for_lacosmic/'
			flts=glob.glob(temp_dir+'*flt.fits')
			
			for i, flt in enumerate(flts):
				print dir+'/'+flt[-18:]
				if os.path.isfile(dir+'/'+flt[-18:]) == False:
					print 'moving ' + flt + ' to ' + dir
					shutil.move(flt,dir)
				else:
					print 'flt was already in permenant location. deleting'
					os.remove(flt)
			
			###now clean out the /temp_for_lacosmic/flt_cleans
			stuff_in_flt_cleans=glob.glob(temp_dir+'flt_cleans/*')
			

		
			#make a permanent flt_cleans in the filter directory if it doesn't already exist
			if not isdir(str(dir)+'/flt_cleans'):
				os.makedirs(dir+'/flt_cleans')
				print 'making directory ' + dir+'/flt_cleans'
			
			for i,thing in enumerate(stuff_in_flt_cleans):
				
				if '.dat' not in thing: #don't overwrite photcat just in case !
					print 'moving ' + thing + ' to ' + dir+'/flt_cleans'
					if os.path.isfile(dir+'/flt_cleans/'+thing[-23:]) == False:
						shutil.move(thing,dir+'/flt_cleans')
						print 'moving clean'
					else:	
						print 'clean was already in permenent location. deleting'
						os.remove(thing)
				#else:
					#os.remove(thing)
					#print 'removing temp catalog in temp /flt_cleans'


		
			###same process for flt_masks
			stuff_in_flt_masks=glob.glob(temp_dir+'flt_masks/*')	
		
			if not isdir(str(dir)+'/flt_masks'):
				os.makedirs(dir+'/flt_masks')
				print 'making directory ' + dir+'/flt_masks'
		
			for i, thing in enumerate(stuff_in_flt_masks):
				print 'moving ' + thing + ' to ' + dir+'/flt_masks'
				print dir+'/flt_masks/'+thing[-22:]
				if os.path.isfile(dir+'/flt_masks/'+thing[-23:])==False:
					shutil.move(thing,dir+'/flt_masks')
					print 'moving mask'
				else:
					print 'mask was already in permenent location. deleting'
					os.remove(thing)

	



if __name__=='__main__':
	directory_with_flts=''
	directory_with_cleans=''
	directory_with_masks=''
	flts = glob.glob(directory_with_flts+'*flt.fits')
	cleans=glob.glob(directory_with_cleans+'*')
	cleans=[clean for clean in cleans if '.dat' not in clean] #prevent from accidentally copying
															 #some temporary phot catalog
	masks=glob.glob(directory_with_masks+'*mask.fits')
	
	