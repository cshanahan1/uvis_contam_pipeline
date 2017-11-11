"""Reads in all the stats files from the monte-carlo simulations and outputs the info to
a text file. This is to make generating a latex table for the ISR easier so I don't 
have to open each file!"""

import glob
import numpy as np
from set_paths_and_params import *

plot_output_dir = set_paths_and_params()['plot_output_dir']+'/'

#read in my stats file
#obj,amp,filter,slope_withWeights,b_withWeights,slope_noWeights,b_noWeights,stdev
my_objs=[]
my_amps=[]
my_filters=[]
my_slopes=[]



with open(set_paths_and_params()['plot_output_dir']+'/stats.dat','r') as f:
	lines=f.readlines()[1:]
	for line in lines:
		line=line.split(',')
		my_objs.append(line[0])
		my_amps.append(line[1])
		my_filters.append(line[2])
		my_slopes.append(line[3])
		

stars=['GD153','GRW70']

objects=[]
filters=[]
amps=[]
means=[]
stdevs=[]


for star in stars:
	print 'reading in files for ' + star 
	path=plot_output_dir+'monte_output/{0}'.format(star)
	stats_files_list=glob.glob(path+'/*slopes_stats.dat')
	print stats_files_list

	for file in stats_files_list:
		filee=file.replace(path,'')
		filee=filee.replace('_100000_slopes_stats.dat','')
		filter=filee.split('_')[0].replace('/','')
		if filter == 'F850L':
			filter='F850LP' #hacky fix its late leave me alone
		amp=filee.split('_')[1]
		
		objects.append(star)
		amps.append(amp)
		filters.append(filter)
		
 	
		with open(file,'r') as ff:
			dat= ff.readlines()[1].split()
			means.append(dat[0])
			stdevs.append(dat[1])

			
print len(objects)	
print means		


with open(plot_output_dir+'monte_output/GD153_monteCarloTable.txt', 'w') as ff:
	ff.write('Filter & Amp & slope_yr & mean_mc & stdev\n')
		
with open(plot_output_dir+'/monte_output/GRW70_monteCarloTable.txt', 'w') as ff:
	ff.write('Filter & Amp & slope_yr & mean_mc & stdev\n')

for i, star in enumerate(objects):
	
	for j, obj in enumerate(my_objs):
		if (star==obj) & (my_filters[j]==filters[i]) & (my_amps[j]==amps[i]):
			print obj
			print amps[i]
			print filters[i]
			my_slope=float(float(my_slopes[j]) * 365. )
			my_slope='%.4f'%(my_slope)
			print 'my slope ' , my_slope
			print 'mc mean  ' , float(means[i])*365
			print '-----------------------'
	filter=filters[i]
	amp=amps[i]
	mean='%.4f'%(float(means[i])*365)
	stdev='%.4f'%(float(stdevs[i])*365)
			
	if star=='GD153':	
		with open(plot_output_dir+'monte_output/GD153_monteCarloTable.txt', 'a') as ff:
			if amp=='A':
				ff.write('{0} & A & {1} & {2} & {3}\n'.format(filter,my_slope,mean,stdev))
			if amp=='C':
				ff.write('{0} & C & {1} & {2} & {3}\n'.format(filter,my_slope,mean,stdev))
				
	if star=='GRW70':	
		with open(plot_output_dir+'monte_output/GRW70_monteCarloTable.txt', 'a') as ff:
			if amp=='A':
				ff.write('{0} & A & {1} & {2} & {3}\n'.format(filter,my_slope,mean,stdev))
			if amp=='C':
				ff.write('{0} & C & {1} & {2} & {3}\n'.format(filter,my_slope,mean,stdev))
		


