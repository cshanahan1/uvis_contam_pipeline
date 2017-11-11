"""Monte-Carlo simulation to find error in flux vs time
slopes. 

Results to be published with UVIS contam and stability
ISR (2014).

Author:

    C.M. Gosmeyer, July 2014
    C. Shanahan, Oct 2016
    
Use:

    cd to directory where you want files generated.

    >>> python montecarlo_fluxvstime.py
    
Outputs:

    For each filter,
    
        PNG file. `<filter>_<amp>_<max_iterations>_histo.png`
            Histogram of the slopes.
        text file. `<filter>_<amp>_<max_iterations>_points.png`
            Randomly selected indices from `deltaflux_ls` array.
        text file. `<filter>_<amp>_<max_iterations>_slopes.png`
            Slope and slope error for each iteration.
        text file. `<filter>_<amp>_<max_iterations>_stats.png`
            Statistics of random sample: average, standard dev,
            median, max, min, and number of bins in histogram.

Notes:

    Dependent on file structure pre-new flat build.
"""

# import module that sets formatting parameters
from matplotlib import rc
# change default font for all plot text to LaTeX font; also change font size
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'size': 14})
# allow TeX commands to be used in text strings
rc('text', usetex=True)
    

import numpy as np
import matplotlib.pyplot as pyplot
import glob
import os
import shutil
from os import remove
from astropy.io import ascii
from random import randint, seed
import pandas as pd 
from set_paths_and_params import *


#-------------------------------------------------------------------------------#  
def calculate_divisor_uvis(divisor_type, amp_letter, amp, flux, expt, mjd, MJD_min, MJD_max):

    cutrange = cutout_mjd_range(amp, mjd, MJD_min, MJD_max)
    if len(cutrange) < 3:
        return 1.0
    
    MJD=mjd[cutrange]
    flux=flux[cutrange]/expt[cutrange]
    divisor = np.median(flux[MJD<=(min(MJD)+50)])
    
    return divisor  
    
def calculuate_deltaflux_errors(flux, expt, merr, divisor, divisor_type):

    # Constant through iteration 
    flux_ref = divisor
    if divisor_type == 'firstimage' or len(flux) < 3:
        # The ref image is Fref = Flux[0] / expt[0]
        # The ref image error is 
        # sig_ref = sqrt[  [(dfref/dflux[0])Fref * sigflux[0] ]^2 + [(dfref/dexpt[0])Fref * sigexpt[0] ]^2 ]
        #         = sqrt[  [1./expt[0] * sigflux[0] ]^2 + [Flux[0] / (expt[0])^2  * sigexpt[0] ]^2 ]
        # We'll suppose that error in exposure time is negligibly small and 
        # the error in flux [in e-] is megerr[0] * Flux[0] / log(e). 
        # Therefore, error in a single reference image reduces to:
        sig_ref = (1./expt[0]) * (merr[0]*flux[0]/1.0857 )
        
    elif divisor_type == 'avfirst3images':    
            
        sig_ref_0 = merr[0]*flux[0]/1.0857
        sig_ref_1 = merr[1]*flux[1]/1.0857
        sig_ref_2 = merr[2]*flux[2]/1.0857
        
        sig_ref = 0.33333333 * (np.sqrt( (sig_ref_0 / expt[0])**2 + \
                 (sig_ref_1 / expt[1])**2 + (sig_ref_2 / expt[2])**2  ))
        
    elif divisor_type == 'median':
        # The ref image is the median of the set. 
        # Assume the standard deviation of the entire set provides the
        # approximate error of the median. 
        stdev_fluxrate = np.std(flux/expt)
        sig_ref = stdev_fluxrate / flux_ref
    
    
    flux_meas = flux/expt

    sig_meas =  merr * (flux / 1.0857) 

    sigma_deltaflux = np.sqrt( ((1.0/flux_ref)**2) * ((sig_meas/expt)**2 + \
                      (flux_meas/flux_ref)**2 * sig_ref**2)) 

    return sigma_deltaflux

def cutout_mjd_range(amp_ls, mjd_ls, mjd_min, mjd_max):
    """Cuts out all indices that fall within given MJD range for
    a given amplifier.
    
    Parameters:
        amp_ls : array of indices
            Amp array for specific letter read from 
            ``<filter>_photcat.dat``.
        mjd_ls : array of floats
            MJDs read from ``<filter>_photcat.dat``.
        mjd_min : float
            Minimum value of MJD to plot after.
        mjd_max : float
            Maximum value of MJD to plot before.
        
    Returns:
        cutrange : list of indices
            Indices corresponding to desired mjd range.
    
    Outputs:
        nothing
    """
    cutrange = []
   
    if mjd_min <= 54900.:
        mjd_min = 55000   

    for index in amp_ls:
        if mjd_ls[index] >= mjd_min and mjd_ls[index] <= mjd_max:
            cutrange.append(index)
            
    return cutrange



def load_data(path_to_photcat,filterr):
    """Loads data from ``<filter>_photcat.dat`` file.
    
    Parameters:
        path_to_photcat : string
            Path to ``<filter>_photcat.dat`` file. 
        detector : string
            Either 'IR' or 'UVIS'.
    
    Returns:
        ampA_ls : array of ints
            Indices corresponding to amp A.
        ampC_ls : array of ints
            Indices corresponding to amp C.            
        mjd_ls : array of floats
            MJD dates of observations.
        flux_ls : array of floats
            Fluxes of observations at aperture 10 pixels.  
        merr_ls : array of floats
            Magnitude errors of observations at aperture 10 pixels.
        expt_ls : array of floats
            Exposure times of observations.          
        imagefilename_ls : array of strings
            Names of the FLT files.
        filtername : string
            Name of the filter.
    
    Outputs:
        nothing
    """
    detector='UVIS'
    pathname = glob.glob(path_to_photcat + filterr +'_photcat.dat')[0]  
    print 'looking for photcat'
    if path_to_photcat == '':
        photcat_filename = glob.glob(filterr+'*_photcat.dat')[0]  
    else: 
        photcat_filename = (pathname.split('/'))[len(pathname.split('/'))-1]
    print photcat_filename
    
    if photcat_filename[5:6].upper() == 'LP':
        filtername = photcat_filename[0:6]
    elif photcat_filename[0:1].upper() == 'G':
        filtername = photcat_filename[0:4]
    else:
        filtername = photcat_filename[0:5]

    if detector.upper() == 'UVIS':    
        # Read in photcat file
            
        df = pd.read_csv(pathname)
        mjd_ls = np.array(df['mjd'].tolist())
        chip_ls = np.array(df['chip'].tolist()) 
        AXIS1, AXIS2 = np.array(df['axis1'].tolist()) ,np.array(df['axis2'].tolist()) 
        XC, YC = np.array(df['xc'].tolist()) ,np.array(df['yc'].tolist()) 
        BACK = np.array(df['back'].tolist())
        expt_ls, flux_ls = np.array(df['exptime'].tolist()) ,np.array(df['f10'].tolist()) 
        merr_ls = np.array(df['m10err'].tolist()) 
        imagefilename_ls = np.array(df['#filename'].tolist())
        shut_ls, amp_ls= np.array(df['shutter'].tolist()),np.array(df['amp'].tolist()) 
    
        # List of indexes for amps.
        ampA_ls = np.where(amp_ls == 'A')[0]
        ampC_ls = np.where(amp_ls == 'C')[0]

        return ampA_ls, ampC_ls, mjd_ls, flux_ls, merr_ls, expt_ls, imagefilename_ls, filtername
        
   


#-------------------------------------------------------------------------------#  

def calculate_flux(amp_ls, mjd_ls, flux_ls, merr_ls, expt_ls, imagefilename_ls, \
                    mjd_min, mjd_max, divisor_type, amp_letter, filtername, \
                    detector, remove_dips=False):
    """Calculates the flux and flux error for the desired MJD range 
    and amplifier subset.
    
    Option to skip the flux dips in blue filters F218W, F225W, and F275W.
    
    Parameters:
        amp_ls : array of ints
            Indices corresponding to given amplifier.            
        mjd_ls : array of floats
            MJD dates of observations.
        flux_ls : array of floats
            Fluxes of observations at aperture 10 pixels.  
        merr_ls : array of floats
            Magnitude errors of observations at aperture 10 pixels.
        expt_ls : array of floats
            Exposure times of observations.          
        imagefilename_ls : array of strings
            Names of the FLT files.
        mjd_min : float
            Minimum MJD desired. 
        mjd_max : float
            Maximum MJD desired.
        divisor_type : string
            What the fluxes should all be divided by -- what decides 
            where 'zero' on the plot sits.
            'firstimage', 'avfirst3images', and 'median' presently 
            only supported options.
        amp_letter : string
            Letter of the given amplifier. 
        filtername : string
            Name of the filter.
        detector : string
            Either 'IR' or 'UVIS'.
        remove_dips : {True, False}
            Off by default. Switch on if you want to skip the flux dips
            in blue filters F218W, F225W, and F275W.
    
    Returns:
        mjd_cutrange_ls : array of floats
            MJDs in the given min, max, and amp limits.
        deltaflux_ls : array of floats
            Flux differences for the given MJD and amp limits.
        deltaflux_err_ls : array of floats
            Error in flux differences for the given MJD and amp limits.
        rootname_ls : array of strings
             Rootname of files for the given MJD and amp limits.
        
    Outputs:
        nothing
    """
    if detector.upper() == 'UVIS':
        divisor = calculate_divisor_uvis(divisor_type, amp_letter, amp_ls, flux_ls, \
                                         expt_ls, mjd_ls, mjd_min, mjd_max)

        if remove_dips == True and mjd_min == 55000. and filtername in ['F218W', 'F225W', 'F275W']:
            print "ENTERING REMOVE DIPS"
            cutrange = remove_flux_dips(filtername, amp_ls, mjd_ls, mjd_min, \
                                        mjd_max, flux_ls, expt_ls, divisor)
        
        else:
            cutrange = cutout_mjd_range(amp_ls, mjd_ls, mjd_min, mjd_max)
            
        print cutrange
    
        mjd_cutrange_ls = mjd_ls[cutrange]    
        deltaflux_ls = ((flux_ls[cutrange]/expt_ls[cutrange])/divisor - 1.0)*100.
        deltaflux_err_ls = calculuate_deltaflux_errors(flux_ls[cutrange], \
                                                       expt_ls[cutrange], \
                                                       merr_ls[cutrange], divisor, \
                                                       divisor_type)

        imagefilename_cutrange_ls = imagefilename_ls[cutrange]
        rootname_ls = []
        for filename in imagefilename_ls:
            rootname_ls.append(filename.split('_flt.clean.fits')[0])
        rootname_ls = np.array(rootname_ls)
    

    
    return mjd_cutrange_ls, deltaflux_ls, deltaflux_err_ls, rootname_ls


#-------------------------------------------------------------------------------#  

def run_montecarlo(amp_ls, mjd_cutrange_ls, deltaflux_ls, \
                   deltaflux_err_ls, rootname_ls, filtername, max_iterations, \
                   random_ls_size, amp_letter, dest):
    """Runs the monte carlo simulation to find error in flux difference
    vs mjd slopes. 
    
    i.e., randomly selects a subset of the indices available in the given
    amp and MJD limits, therefore returning random flux-MJD pairs, from
    which a slope is calculated. This is repeated for `max_iterations`. 
    A histogram and text files containing statistical information are 
    generated at the end.
    
    Parameters:
        amp_ls : array of ints
            Indices corresponding to given amplifier.            
        mjd_cutrange_ls : array of floats
            MJDs in the given min, max, and amp limits.
        deltaflux_ls : array of floats
            Flux differences for the given MJD and amp limits.
        deltaflux_err_ls : array of floats
            Error in flux differences for the given MJD and amp limits.
        rootname_ls : array of strings
            Rootname of files for the given MJD and amp limits.
        filtername : string
            Name of filter.
        max_iterations : int
            Number of iterations to run simulation.    
        random_ls_size : int
            Number of indices from the whole dataset (selected from
            `deltaflux_ls`) to randomly draw.
        amp_letter : string
            Letter of the given amplifier. 
        dest : string
            Path to where you want output files generated.
    
    Returns:
        nothing
    
    Outputs:
        PNG file. `<filter>_<amp>_<max_iterations>_histo.png`
            Histogram of the slopes.
        text file. `<filter>_<amp>_<max_iterations>_points.png`
            Randomly selected indices from `deltaflux_ls` array.
        text file. `<filter>_<amp>_<max_iterations>_slopes.png`
            Slope and slope error for each iteration.
        text file. `<filter>_<amp>_<max_iterations>_stats.png`
            Statistics of random sample: average, standard dev,
            median, max, min, and number of bins in histogram.
    """

    slope_data_ls = []
    points_ls = []
    count = 0
    
    while count < max_iterations:
        
        # Create list of indexes of the size of the flux list
        flux_index_ls = range(0, len(deltaflux_ls))
            
        # Randomly select subset of observations via the flux list
        flux_random_ls = random_select_subset_without_repeats(flux_index_ls, random_ls_size)
        
        # Save the random indices to a list.
        # Can then use these indexes on a amp-truncated list of fluxes
        # and MJDs to recall which observations were used. 
        points_ls.append(flux_random_ls)
        
        # Cut the randomly selected observations from flux and time lists
        mjd_random_ls = mjd_cutrange_ls[flux_random_ls]
        deltaflux_random_ls = deltaflux_ls[flux_random_ls]
        deltaflux_err_random_ls = deltaflux_err_ls[flux_random_ls]
        rootname_random_ls = rootname_ls[flux_random_ls]
        
        #run_tcreate(deltaflux_random_ls, mjd_random_ls, deltaflux_err_random_ls, rootname_random_ls, dest) 
        #run_tlinear(filtername, amp_letter, outfilemark='_montecarlo_'+str(count), dest=dest)                    
        
        # Read in *hdr.tlinear.dat file
        #tlinear_data = ascii.read(dest + filtername+'_'+amp_letter+'_montecarlo_'+str(count)+'_tlinear.hdr.dat', names=['col1','col2','col3'])
        #slope = float(tlinear_data['col3'][5])
        slope,b=np.polyfit(mjd_random_ls,deltaflux_random_ls,1,w=deltaflux_err_random_ls)
        p,slope_err,_1,_1,_1=np.polyfit(mjd_random_ls,deltaflux_random_ls,1,full=True,w=deltaflux_err_random_ls)
        #slope_err = float(tlinear_data['col3'][7])
        #remove(dest + filtername+'_'+amp_letter+'_montecarlo_'+str(count)+'_tlinear.hdr.dat')
        #remove(dest + filtername+'_'+amp_letter+'_montecarlo_'+str(count)+'_tlinear.dat')
        
        # Store the iteration, slope, and slope err
        slope_data_ls.append([int(count), slope, slope_err, slope*365., slope_err*365.])
        print '-----------------------'
        print count#, slope, slope_err, slope*365., slope_err*365.
                        
        count += 1

    # Write slopes, errors, etc. to file
    print "Writing slopes to file..."
    ascii.write(np.array(slope_data_ls), \
                dest + filtername+'_'+amp_letter+'_'+str(max_iterations)+'_slopes.dat', \
                names=['#Count', 'Slope', 'Slope_Err', 'Slope/Year', 'Slope_Err/Year'])
    
    # Write random points into file
    print "Writing random indices to file..."
    ascii.write(np.array(points_ls), dest + filtername+'_'+amp_letter+\
                                     '_'+str(max_iterations)+'_points.dat')
    
    # Do some statistics: average, std dev, & median slope 
    slopes = []
    slopesyear = []
    for slope_data in slope_data_ls:
        slopes.append(slope_data[1])
        slopesyear.append(slope_data[3])
    
    print "---Slopes:"
    print slopes
    
    print "---Slopes/Year:"
    print slopesyear
    
    av_slopes = np.average(slopes)
    print "---Average Slope:"
    print av_slopes
    
    av_slopesyear = np.average(slopesyear)
    print "---Average Slope/Year:"
    print av_slopesyear
    
    stddev_slopes = np.std(slopes)
    print "---Std Dev Slope:"
    print stddev_slopes
    
    stddev_slopesyear = np.std(slopesyear)
    print "---Std Dev Slope/Year:"
    print stddev_slopesyear
    
    median_slopes = np.median(slopes)
    print "---Median Slope:"
    print median_slopes
    
    median_slopesyear = np.median(slopesyear)
    print "---Median Slope/Year:"
    print median_slopesyear
    
    
    # Create slope histogram.
    print "Creating histogram of slopes..."
    min_slope, max_slope = np.min(slopes), np.max(slopes)
    print "---Max, Min Slope:"
    print min_slope, max_slope
    binsize = .00002
    num_bins = abs(np.floor((max_slope - min_slope) / binsize))
    fig, ax = pyplot.subplots()
    
    print "---Number of bins:"
    print num_bins
    ax.set_xlabel(r'Slope')
    ax.set_title(r'{0} , {1} , {2}'.format(filtername,amp_letter,str(max_iterations)))
    
    n, bins, patches = ax.hist(slopes, num_bins, fc='royalblue', ec='white')
    for b in patches[::2]: b.set_facecolor('blue')
    
    fig.savefig(dest + filtername+'_'+amp_letter+'_'+str(max_iterations)+'_slopes_histo.png')
    pyplot.close()
    
    # Write slope statistics to file.
    print "Writing slope statistics to file..."
    stats = [[av_slopes], [stddev_slopes], [median_slopes], [max_slope], [min_slope], [num_bins]]
    ascii.write(stats, dest + filtername+'_'+amp_letter+'_'+str(max_iterations)+'_slopes_stats.dat', \
                names=['#Average', 'Std_Dev', 'Median', 'Max', 'Min', 'Num_Bins'])
                
                
    # Create slope/year histogram.
    print "Creating histogram of slopes/year..."
    min_slopeyear, max_slopeyear = np.min(slopesyear), np.max(slopesyear)
    print "---Max, Min Slope/Year:"
    print min_slopeyear, max_slopeyear
    binsize = .008
    num_bins = abs(np.floor((max_slopeyear - min_slopeyear) / binsize))
    fig, ax = pyplot.subplots()
    
    print "---Number of bins:"
    print num_bins
    ax.set_xlabel('Slope/Year')
    ax.set_title(filtername + ', ' + amp_letter + ', ' + str(max_iterations))
    
    n, bins, patches = ax.hist(slopesyear, num_bins, fc='royalblue', ec='white')
    for b in patches[::2]: b.set_facecolor('blue')
    
    fig.savefig(dest + filtername+'_'+amp_letter+'_'+str(max_iterations)+'_slopes-year_histo.png')
    pyplot.close()
    
    # Write slope/year statistics to file.
    print "Writing slope/year statistics to file..."
    stats = [[av_slopesyear], [stddev_slopesyear], [median_slopesyear], [max_slopeyear], \
            [min_slopeyear], [num_bins]]
    ascii.write(stats, dest + filtername+'_'+amp_letter+'_'+str(max_iterations)+'_slopes-year_stats.dat', \
                names=['#Average', 'Std_Dev', 'Median', 'Max', 'Min', 'Num_Bins'])


#-------------------------------------------------------------------------------#  

def random_select_subset_without_repeats(data_ls, size):
    """Pseudo-randomly selects a subset of elements
    (for a total number set by `size`) from a list,
    without repeats.  
    
    Parameters:
        data_ls : array
            Data from which want to randomly select elements.
        size : int
            How large you desire `random_ls` to be.
            Must be smaller than length of `data_ls`. 
        
    Returns:
        random_ls : list 
            List of randomly selected elements from `data_ls`.
            
    Outputs:
        nothing
    
    Notes:
        Adapted from
        http://code.activestate.com/recipes/59883-random-selection-of-elements-in-a-list-with-no-rep/
    """
    random_ls = []
    
    print data_ls
    iterator = 0
    while iterator < size:
        index = randint(0, len(data_ls)-2-iterator)
        elem = data_ls[index]
        data_ls[index] = data_ls[len(data_ls)-2-iterator]
        random_ls.append(elem)
        iterator += 1
    
    return random_ls    


#-------------------------------------------------------------------------------#  

def move_files(dir_name):
    """Creates directory `dir_name` and moves the histograms 
    and corresponding text files to it.
    
    Parameters:
        dir_name : string
            Name of directory to be created.
    
    Returns:
        nothing
        
    Outputs:
        nothing
    """
    file_ls = glob.glob('*dat')   
    plot_ls = glob.glob('*png')
    file_ls += plot_ls
    os.makedirs(dir_name)
    for file in file_ls:
        shutil.move(str(file), dir_name)

#-------------------------------------------------------------------------------#  
def single_filter(path_to_photcat, max_iterations, detector, mjd_min, \
                  mjd_max, remove_dips, dest,filterr):
    """Performs the monte carlo function suite over both
    amps A and C, separately, for a single filter.
    
    Parameters:
        path_to_photcat : string
            Path to the ``<filter>_photcat.dat`` file.
        max_iterations : int
            Number of iterations to run simulation.
        mjd_min : float
            Minimum MJD from which to select data subset. 
        mjd_max : float
            Maximum MJD from which to select data subset.   
        remove_dips : {True, False}
            Off by default. Switch on if you want to skip the flux dips
            in blue filters F218W, F225W, and F275W.
        dest : string
            Path to where you want output files generated.
                
    Returns:
        nothing
        
    Outputs:
        See :func:`run_montecarlo`.
    """
    if detector.upper() == 'UVIS':
        ampA_ls, ampC_ls, mjd_ls, flux_ls, merr_ls, expt_ls, imagefilename_ls, filtername = \
            load_data(path_to_photcat,filterr)
    
        mjd_cutrange_lsA, deltaflux_lsA, deltaflux_err_lsA, rootname_lsA = \
            calculate_flux(ampA_ls, mjd_ls, flux_ls, merr_ls, expt_ls, imagefilename_ls, \
                           mjd_min, mjd_max, 'avfirst3images', 'A', filtername, \
                           detector, remove_dips)
        mjd_cutrange_lsC, deltaflux_lsC, deltaflux_err_lsC, rootname_lsC = \
            calculate_flux(ampC_ls, mjd_ls, flux_ls, merr_ls, expt_ls, imagefilename_ls, \
                           mjd_min, mjd_max, 'avfirst3images', 'C', filtername, \
                           detector, remove_dips)

        random_ls_sizeA = int(len(deltaflux_lsA)/2)
        random_ls_sizeC = int(len(deltaflux_lsC)/2)

        run_montecarlo(ampA_ls, mjd_cutrange_lsA, deltaflux_lsA, \
                       deltaflux_err_lsA, rootname_lsA, \
                       filtername, max_iterations, random_ls_sizeA, 'A', output_dir)
        run_montecarlo(ampC_ls, mjd_cutrange_lsC, deltaflux_lsC, \
                       deltaflux_err_lsC, rootname_lsC, \
                       filtername, max_iterations, random_ls_sizeC, 'C', output_dir)


                   
#-------------------------------------------------------------------------------#  

                                     
#-------------------------------------------------------------------------------#  

if __name__ == '__main__':
    paths = set_paths_and_params()
    data_dir=paths['data']+'/data'
    max_iterations = 10000
    plot_output_dir = paths['plot_output_dir']
    detector = 'UVIS'
    stars = paths['objects']
    
    if not os.path.isdir(plot_output_dir + '/monte_output'):
        os.mkdir(plot_output_dir + '/monte_output')
    
    
    for star in stars:
		if not os.path.isdir(plot_output_dir + '/monte_output/{}'.format(star)):
			os.mkdir(plot_output_dir + '/monte_output/{}'.format(star))

		if star == 'GD153':
			filters=['F218W','F225W', 'F336W', 'F606W','F275W','F438W','F814W']
			photcats=[data_dir+'/GD153/'+filt+'/flt_cleans/' for filt in filters]
			output_dir = plot_output_dir + '/monte_output/{}/'.format(star)

		if star == 'GRW70':
			filters=['F218W','F225W', 'F336W', 'F606W','F275W','F438W','F814W','F343N',\
			'F373N','F390M','F390W','F395N','F410M','F467M','F547M','F555W','F850LP']
			photcats=[data_dir+'/GRW70/'+filt+'/flt_cleans/' for filt in filters]
			output_dir = plot_output_dir + '/monte_output/{}/'.format(star)

		for i, photcat in enumerate(photcats):
			print photcat
			if os.path.isdir(photcat.replace('flt_cleans/','flt_cleans')):
				print 'ok'
				single_filter(photcat, max_iterations, detector, 0.0, 58000, False, output_dir,filters[i])
