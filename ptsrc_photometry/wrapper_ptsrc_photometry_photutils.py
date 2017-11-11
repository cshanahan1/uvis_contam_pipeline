import glob
from ptsrc_photometry_photutils import *
from set_paths_and_params import *


def wrapper_ptsrc_photometry_photutils(paths, PAMcorr = True,temp_dir = True):

    objects = paths['objects']
    filters = paths['filters']
    data_dir = paths['data']+'/data'
    
    
    if temp_dir:
        strr = 'temp_for_photometry'
    else:
        strr = ''
            
    for ob in objects:
        for filterr in filters:
            dir = data_dir+'/{0}/{1}/{2}/flt_cleans'.format(ob,filterr,strr)
            phot_file_path = data_dir+'/{0}/{1}/flt_cleans/{1}_photcat.dat'.format(ob,filterr,strr)
            
            ifiles = glob.glob(dir+'/*clean.fits')
            
            if len(ifiles) > 0:           
                print 'beginning photometry on {0}, {1}'.format(ob,filterr)
                main_ptsrc_photometry(ifiles,phot_file_path, paths, PAMcorr)
            else:
                print 'no files. moving on to next filter.'
                continue
            
            
if __name__ == '__main__':

    paths = set_paths_and_params()
    
    wrapper_ptsrc_photometry_photutils(paths,temp_dir = True)
    
    