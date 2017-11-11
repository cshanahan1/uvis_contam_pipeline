#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from set_paths_and_params import * 
# import module that sets formatting parameters
from matplotlib import rc
# change default font for all plot text to LaTeX font; also change font size
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'size': 14})
# allow TeX commands to be used in text strings
rc('text', usetex=True)


""" 

	Makes plots of slope vs. filter. Seperate plot for each star, and different markers
	for quad A and quad C slopes.
	
	Reads in the 'stats.dat' with the linear fits to the data for all filters / chips.
	This file is created for all data by running the main plotting script 
	'uvis_contam_plots_overFilters.py'. 
	
	
"""

def choose_filters():
	""" Set which filters you wish to include in the plot
		Put in order of increasing wavelength"""
	
	filters = ['F218W', 'F225W', 'F275W', 'F336W', 'F438W','F606W','F814W']
		
	return filters 
	
def read_in_stats_file(ifile):
	""" Parses stats.dat file. Filters data to only include rows with the filters selected
		in 'choose_filters()'.
	
		Returns:
			pandas dataframe with slope, filter, amp (A or C), and star
	"""

	df = pd.read_csv(ifile)

	#clunky way of unpacking but numpy.loadtxt doesn't like dealing with 
	#a mix of strings and floats in a csv so this is easier
	filters=np.array(df['filter'].tolist())
	slopes=np.array(df['slope_yr'].tolist())
	amps=np.array(df['amp'].tolist())
	objs=np.array(df['#obj'].tolist())
	
	selected_filters=choose_filters() # set above, what filters to make plot with
	
	for i, filter in enumerate(filters):
		if filter not in selected_filters:
			slopes[i]=-999
			
	
	filters = filters[slopes!=-999]
	amps = amps[slopes!=-999]
	objs = objs[slopes!=-999]
	slopes = slopes[slopes!=-999]
	
	filter_indicies = []#used for plotting strings on x axis
	for filter in filters:
		for i , val in enumerate(selected_filters):
			if filter == val:
				filter_indicies.append(i)
	
	dat={'filters':filters,'amps':amps,'objs':objs,'slopes':slopes,\
		'filter_indicies':filter_indicies}
	df =  pd.DataFrame(data=dat,index=None)
	
	return df

def read_montecarlo(obj):
	""" Reads in the montecarlo data table for either GD153 and GRW70.
	
		make sure the header for that file is:
		
		Filter , Amp , slope_yr , mean_mc , stdev.
		
		The existing script that generates the table outputs a different header to be pasted into 
		latex so swap it with that.
		
		Returns:
			pandas dataframe with filters, chips, and the MC standard deviations
		 
		  """
	plot_output_dir = set_paths_and_params()['plot_output_dir']
	datfile=plot_output_dir+'/monte_output/'+obj+'_monteCarloTable.txt' #path to where montecarlo table lives
	
	
	df = pd.read_csv(datfile, delimiter = ' & ')
	print df
	
	filters=np.array(df['Filter'].tolist())
	amps=np.array(df['Amp'].tolist())
	stdevs=np.array(df['stdev'].tolist())
	
	dat = {'filters':filters,'amps':amps,'stdevs':stdevs}
	
	df =  pd.DataFrame(data=dat,index=None)
	
	return df
	
def return_matching_stdev(mc_data,amp,filterr):

	stdevs = np.array(mc_data['stdevs'].tolist())
	amps = np.array(mc_data['amps'].tolist())
	filters = np.array(mc_data['filters'].tolist())
	
	print stdevs[(amps==amp)&(filters==filterr)][0]
	return stdevs[(amps==amp)&(filters==filterr)][0]

def make_plots(df, obj):	

	df = df[df['objs']==obj]

	slopes=df['slopes'].values
	amps=df['amps'].values
	filters = df['filters'].values 
	filter_indicies=df['filter_indicies'].values
	
	slopes_chip_A = slopes[amps=='A']
	filters_chip_A = filters[amps=='A']
	print filters_chip_A
	filter_indicies_chipA=filter_indicies[amps=='A']
	
	print 'chip A' 
	for i, filt in enumerate(filters_chip_A):
		print filt, ' ' , slopes_chip_A[i]

	slopes_chip_C = slopes[amps=='C']
	filters_chip_C =filters[amps=='C']
	print filters_chip_C
	filter_indicies_chipC=filter_indicies[amps=='C']
	
	print 'chip C'
	for i, filt in enumerate(filters_chip_C):
		print filt, ' ' ,slopes_chip_C[i]
	
	#now find errorbars
	errs_a=[]
	errs_c=[]
	mc_data=read_montecarlo(obj)
	
	for val in filters_chip_A:
		print val
		print return_matching_stdev(mc_data,'A',val)
		errs_a.append(return_matching_stdev(mc_data,'A',val))
		
	for i, val in enumerate(filters_chip_A):
		errs_c.append(return_matching_stdev(mc_data,'A',val))
		
	#plt.scatter(filter_indicies_chipA,slopes_chip_A, c='r',marker='^', s=40,label = 'Quad. A',yerr=errs_a)
	#plt.scatter(filter_indicies_chipC,slopes_chip_C,marker='^', s=40, c='b',label = 'Quad. C',yerr_errs_c)
	
	plt.errorbar(filter_indicies_chipA,slopes_chip_A,yerr=errs_a,fmt='r^',label = 'Quad. A')
	plt.errorbar(filter_indicies_chipC,slopes_chip_C,yerr=errs_c,fmt='b^',label = 'Quad. C')
	

	selected_filters=choose_filters()
	plt.xticks(np.arange(len(selected_filters)),selected_filters)
	plt.ylabel(r'$\Delta$ Flux [\% per year]',fontsize=20)
	plt.ylim(-0.7,0.2)
	plt.xlim(-0.5,7)
	plt.xlabel(r'Filter',fontsize=20)
	plt.text(0.02,0.07,obj,fontsize=20)
	plt.legend(loc='best')
	plt.savefig(set_paths_and_params()['plot_output_dir']+'/'+obj+'_slope_vs_filter.png')
	plt.show()

def main_make_plot_slope_vs_filter(ifile):

	df = read_in_stats_file(ifile)
	
	make_plots(df, 'GRW70')
	make_plots(df, 'GD153')
	
	
if __name__=='__main__':
	plot_output_dir = set_paths_and_params()['plot_output_dir'] #where the stats.dat file is
	main_make_plot_slope_vs_filter(plot_output_dir+'/stats.dat')
	
	
	
	
