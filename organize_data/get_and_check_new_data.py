from astropy.io import fits
import glob
import os.path

from pyql.database.ql_database_interface import session
from pyql.database.ql_database_interface import Master
from pyql.database.ql_database_interface import UVIS_flt_0
from pyql.database.ql_database_interface import UVIS_flt_1

def query_ql(proposal_id, object_names, file_type='flt'):

	"""
	Querys the quicklook database by proposal ID for the full paths to files. Default file 
	type is _flt.fits but you can specify _raw.fits if you wish to calibrate the data.
	
	Parameters
	----------
	list_new_files : list of str
		A list of all new file paths that you wish to check against existing data.
		
	paths: dictionary 
		The dictionary from set_paths.py that specifys all paths/parameters.
		
	Returns
	-------
	
	new_unique_files : list
		List of all files in list_new_files that aren't already in the data directories.
		
	"""
	objs=[]
	query_results=[]	
	roots=[]
	
	for object in object_names:
		results = session.query(Master.dir,UVIS_flt_0.targname,UVIS_flt_0.filter,\
								Master.rootname).join(UVIS_flt_0).\
								filter(UVIS_flt_0.detector=='UVIS').\
								filter(UVIS_flt_0.proposid==proposal_id).\
								filter(UVIS_flt_0.targname==object).\
								filter(UVIS_flt_0.filter != 'G280' ).all()
		rootnames=[obs.rootname for obs in results]
		for root in rootnames:
			roots.append(root)
		for path in results:
			query_results.append(path[0]) #because the query returns a tuple
		objsss=[obs.targname for obs in results]
		filters=[obs.filter for obs in results]

		for i, ob in enumerate(objsss):
			objs.append(ob)

	
	# all unique directories
	query_results=set(query_results)	
		
	paths=[]
	for dir in query_results:
		files_in_dir = glob.glob(dir + '/*{}.fits'.format(file_type))
		for filee in files_in_dir:
			if filee[-18:-9] in roots: 
				paths.append(filee)
	
	return paths
	
def check_new_files_against_existing_files(list_new_files,data_dir,new_flts_dir,raw_dir):

	"""
	Given a list of files, this function will check the names of these files against 
	all files that are alreadu in the UVIS contam monitor data directories. This function
	will return a list of all the files that don't already exist in these directories 
	so you can ensure you are not copying over any duplicate files
	
	Parameters
	----------
	list_new_files : list of str
		A list of all new file paths that you wish to check against existing data.
		
	paths: dictionary 
		The dictionary from set_paths.py that specifys all paths/parameters.
		
	Returns
	-------
	
	new_unique_files : list
		List of all files in list_new_files that aren't already in the data directories.
		
	"""


	existing_flts=[] 
	
	#####check for existing files in the data directory and in the temp directories #####
	
	
	objects = [os.path.basename(item) for item in glob.glob(data_dir+'/*')]
	
	existing_flts = []
	existing_temp_flts= []
	existing_bad_flts = []
	
	for object in objects:

		filters=[os.path.basename(i) for i in glob.glob(data_dir+'/'+object+'/*')]
						  
		for filterr in filters:
			
			files_directory = '{}/{}/{}'.format(data_dir,object,filterr)
			files = glob.glob(files_directory+'/*flt.fits') 

			temp_files_directory = files_directory+'/temp_for_photometry'
			temp_files = glob.glob(data_dir+'/' + object + '/'+filterr + \
								   '/temp_for_photometry/*flt.fits')
			
			bad_files_directory = files_directory+'/flt_cleans/failed_source_finding'
			bad_files_temp_directory = temp_files_directory+'/flt_cleans/failed_source_finding'
			bad_files = glob.glob(bad_files_directory + '/*/*.fits')
			bad_files = bad_files + glob.glob(bad_files_temp_directory + '/*/*.fits')
			bad_files = bad_files + glob.glob(data_dir + '/bad_files/*/*')
			
			if len(files) > 0: # if you found files in data directory
				for filee in files:
					existing_flts.append(filee)
					
			if len(temp_files) > 0: # if you found files in temp directory
				for filee in temp_files:
					existing_temp_flts.append(filee)
					
			if len(bad_files) > 0: # if you found files in the bad_files directory
				for filee in bad_files:
					existing_bad_flts.append(filee)
			

	#now look for existing flts in the 'new_flt' and 'new_raws' directory 
	
	files_in_new_data = glob.glob(new_flts_dir+'/*.fits')+glob.glob(raw_dir+'/*.fits')
		
	if len(files_in_new_data) > 0:
		for filee in files_in_new_data:
			existing_flts.append(filee)
	
	
	all_existing_files = existing_flts + existing_temp_flts +existing_bad_flts
	
	existing_file_rootnames = [item[-18:-9] for item in all_existing_files]
	new_file_rootnames = [item[-18:-9] for item in list_new_files]
		
	new_uniques=set(new_file_rootnames)-set(existing_file_rootnames)

	new_unique_flts=[]
	for new_file in list_new_files:
		for unique in new_uniques:
			if unique in new_file:
				new_unique_flts.append(new_file)

	return new_unique_flts
	
def check_cal_version(list_of_flts,cal_versions):
	"""
	Checks if files have been calibrated with the correct version of calwfc3.
	
	 Returns:
	 	A 2D list, the first being the flts that are calibrated with the correct
	 	version of the pipeline and the second list of flts that aren't calibrated 
	 	with the correct version.
	 	
	 """
	
	good_flts = []
	bad_flts=[]
	for flt in list_of_flts:
		hdulist=fits.open(flt)
		hdr=hdulist[0].header
		hdulist.close
		print hdr['CAL_VER']
		if hdr['CAL_VER'] in cal_versions:

	
			good_flts.append(flt)
		else:
			bad_flts.append(flt)
		
	return [good_flts,bad_flts]
	