import numpy as np
def append_to_stats_file(obj,filtername,plot_output_dir,mjdA,difA,sigA,mjdC,difC,sigC):

    """
    Appends the stats.dat file with object name, filtername, MJD, and the slope and 
    standard deviation for the linear fit between flux % difference. File is placed in 
    plot_output_dir, which is set it the set_paths.py script. 
    
    This function is only invoked when this plotting script is being run over all filters 
    and objects in the wrapper script uvis_contam_plots_overFilters.py. When this script 
    is run on its own, the call to this function is switched off.
        
        Parameters:
            obj: string
                name of object (GRW70 or GD153)
            filtername : string
                name of filter
            plot_output_dir : string
                path to where the file is written               

        Returns: 
            nothing.
            
        Outputs: 
        
            'stats.dat' file 

        """

    if len(difA)>0:
        slopeA,bA = np.polyfit(mjdA, difA, 1,w=sigA)
        slopeA_yr = slopeA*365.
        stdA = np.std(difA)
        stdA_yr =stdA*365.
        with open(plot_output_dir+'/stats.dat','a') as f:
            f.write('{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(obj,'A',filtername,slopeA,\
                                                               bA,stdA,slopeA_yr,stdA_yr))
    if len(difC)>0:
        slopeC,bC = np.polyfit(mjdC, difC, 1,w=sigC)
        slopeC_yr = slopeC*365.
        stdC = np.std(difC)
        stdC_yr =stdC*365.
        with open(plot_output_dir+'/stats.dat','a') as f:
            f.write('{0},{1},{2},{3},{4},{5},{6},{7}\n'.format(obj,'C',filtername,slopeC,\
                                                               bC,stdC,slopeC_yr,stdC_yr))


    

