"""

This a setup script that should be carefully checked and run before running any of 
the other scripts in the contamination monitor code! 

This contains a dictionary where all necessary paths are specified. If any of these 
directories are missing, they will be created when this script is executed. A copy of this
file containing the paths dictionary will be copied to every code subdirectory for easy
access (without requiring sys.path.append...). 


Author
------
    C. Shanahan, August 2016 
    
Usage
-----

    Specify the location of the:
        Directory containing all the code: "scripts"
        Location of the new .flts (or .raws) : "new_flts"
        Location of permenant data directory for processed data: "data"
        Location where files that fail source finding/plotting will be moved:"bad_files"
        Location of the pixel area maps: "pam"
        List of proposal IDs for UVIS contam/zeropoint programs: "prop_ids"
        The subset (or full set) of filters you wish to run the code over: "filters"
        The subset (or full set) of objects you wish to run te code over: "objects"
        
    The set_up_dirs() function will only run if it is your first time running the code.
    The main data storage directory is set in the 'set_paths' dictionary, and if the 'data'
    and 'new_flt' subdirectories do not exist in that location, they will be created.


    
 """
 
import glob
import os
import shutil

""""filters:['F475W','F218W','F225W', 'F336W', 'F606W','F275W',\
              'F438W','F814W','F280N','F343N','F373N','F390M','F390W','F395N','F410M',\
              'F467M','F469N','F502N','F547M','F555W','F850LP'],
              
              "filters":['F475W','F218W','F225W', 'F336W', 'F606W','F275W',\
              'F438W','F814W','F280N','F343N','F373N','F390M','F390W','F395N','F410M',\
              'F467M','F469N','F502N','F547M','F555W','F850LP'],
              
              "prop_ids":['14021', '13089', '13088', '14018', '13575', '13574', '12699', \
              '14384', '11903', '14815', '11798', '12698', '11907', '14382', '12334', \
              '12333', '11450', '11426'],
"""

def set_paths_and_params():

    """Specify full paths to various output directories here, as well as the full set of 
        filters and stars that the code will be run for. """
        
    #for ease if all you want all your output directories to be subdirectories
    #in one location 
    base_path = "/grp/hst/wfc3r/cshanahan/UVIS_contam_monitor"
        
    paths_and_params = { "scripts":base_path,
              "data":base_path+"/data",
              "pam_dir":"/grp/hst/wfc3t/cshanahan/",
              "plot_output_dir":base_path+"/plots",
               "filters":['F475W','F218W','F225W', 'F336W', 'F606W','F275W',\
              'F438W','F814W','F280N','F343N','F373N','F390M','F390W','F395N','F410M',\
              'F467M','F469N','F502N','F547M','F555W','F850LP'],
              "objects":['GRW70','GD153'],
              "prop_ids":['14815'],
              "cal_ver" : ['3.5','3.4'],
              "aperture_sizes" : [2,3,8,10,12,15,20,25,30,50,60]
              }
              
    return paths_and_params
    
def set_up_dirs():

    """This sets up the necessary 'data', 'new_flts', and 'plot_output_dir"
         subdirectories if they do not exist yet (first time running this code)"""
    
    paths=set_paths_and_params()
  
    data_dir = paths['data']
    
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    data_subdirs = ['/data','/new_flts','/new_raws']
    
    paths = [data_dir+subdir for subdir in data_subdirs] + [paths[key] for key in paths]
    for item in paths:
        if type(item) != list:
            if not os.path.isdir(item):
                print 'Making directory ' + item
                os.mkdir(item)
                
def copy_calwf3_script():

    paths=set_paths_and_params()
    
    print 'copying the latest version of run_calwf3.py to new_raws directory'
    shutil.copy(paths['scripts']+'/organize_data/run_calwf3.py',paths['data']+'/new_raws')
    
    


if __name__ == '__main__':

    origin=set_paths_and_params()['scripts']
    
    subdirs=['organize_data','cosmic_ray_removal','ptsrc_photometry',\
            'plotting_and_statistics']
    dirs = [origin + '/' + subdir for subdir in subdirs]
    for dir in dirs: #copy this script to each of the sub directories with code
        print 'copying set_paths to ' + dir
        shutil.copy('set_paths_and_params.py',dir)
        
    #create 'data' and 'new_flt' subdirectories if they don't exist
    set_up_dirs() 
    
    #copy the latest version of run_calwf3.py to new_raws directory
    copy_calwf3_script()
    
        
