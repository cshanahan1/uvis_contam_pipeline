import numpy as np
import os.path

# import module that sets formatting parameters
from matplotlib import rc
# change default font for all plot text to LaTeX font; also change font size
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'size': 14})
# allow TeX commands to be used in text strings
rc('text', usetex=True)
	
import sys
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import *
from generate_stats_file import *


def load_phot_data(ifile):

    """ Reads in and parses photometry catalog.
        
        Returns:
        
            Pandas dataframe with important values from file (mjd, chip, flux in different 
            radii, etc... All acessable by name.
            
    """
    
    filtername = (os.path.basename(ifile)).replace('_photcat.dat', '')
    
    
    df = pd.read_csv(ifile,index_col = False)
    amp = np.array(df['amp'].tolist())
    back = np.array(df['back'].tolist())
    chip= np.array(df['chip'].tolist())
    expt= np.array(df['exptime'].tolist())
    f10 = np.array(df['f10'].tolist())
    ferr10 = np.array(df['f10err'].tolist())
    m10 = np.array(df['m10'].tolist())
    merr10 = np.array(df['m10err'].tolist())
    imagefilename = np.array(df['#filename'].tolist())
    mjd = np.array(df['mjd'].tolist())
    xc = np.array(df['xc'].tolist())
    yc = np.array(df['yc'].tolist())
    axis1= np.array(df['axis1'].tolist())
    axis2 = np.array(df['axis2'].tolist())
    shut= np.array(df['shutter'].tolist())

    #Measure the radius from center of the image.
    rads = np.sqrt((xc-axis1/2.0)**2 + (yc-axis2/2.0)**2)
    d = {'amp':amp, 'back':back, 'chip':chip, 'expt':expt, 'filtername': filtername, \
    'f10':f10, 'ferr10':ferr10,'imagefilename':imagefilename, 'mjd':mjd, 'rads':rads, \
    'shut':shut,'m10':m10, 'merr10':merr10}
    
    df=pd.DataFrame(data=d)
    #print df['ferr10']
    return df

def generate_divisor(MJD,flux,len_div=50.):

	"""
	
	Calculates the baseline average flux value to which all subsequent measurements are
	compared to measure trends in throughput. This baseline value is defined as the median 
	of all images taken within the first 50 days since the first observation (which should 
	capture all images on the first ~2 visits).
	
	Returns:
		
		median flux value of earliest observations.
		
	"""
		
	#print 'normalizing to ' +str(len(flux[MJD<=(min(MJD)+len_div)])) +  ' images'
	return np.median(flux[MJD<=(min(MJD)+len_div)])

		#return float(flux[MJD==min(MJD)])

def plot_lines(objname,mjdA,difA,mjdC,difC,mjdB,difB,mjdD,difD,sigA,sigC,filtername):

	fig=plt.figure()
	fig.subplots_adjust(hspace=0.1)
	ax1 = fig.add_subplot(1,1,1) #UVIS chip 1, amp A	
	ax1.set_xlim([54800,58500])
	
	if len(mjdA) > 0 :	

		m,b = polyfit(mjdA, difA, 1,w=sigA)
		x=np.linspace(mjdA.min(),mjdA.max(),10)
		ax1.plot(x,m*x+b,c='k',lw=3,label='Chip 1')	
		ax1.errorbar(mjdA,difA,yerr=sigA,linestyle="None",c='k',alpha=0.4)
		ax1.scatter(mjdA, difA, label = 'Chip 1',c='k',s=100, alpha = 0.4)
		ax1.text(0.01,0.06,'Chip 1 Slope , intercept= ' + str(m)[0:20] +', '+ str(b),transform=ax1.transAxes)
		ax1.text(0.4,0.9,'GRW70, '+filtername,transform=ax1.transAxes)
		
	if len(mjdC) > 0 :	
		m,b = polyfit(mjdC, difC, 1,w=sigC)
		ax1.text(0.01,0.02,'Chip 2 Slope = ' + str(m)[0:20] +', '+ str(b),transform=ax1.transAxes)
		x=np.linspace(mjdC.min(),mjdC.max(),10)
		ax1.plot(x,m*x+b,c='r',lw=3, label='Chip 2')
		
	plt.legend(loc='best')
	plt.xlabel(r'MJD',fontsize=20)
	plt.ylabel(r'$\Delta$ Flux [\%]',fontsize=20)
		
	plt.savefig(filtername+'lines.png')
	plt.show()
	

def plotting(objname,mjdA,difA,mjdC,difC,mjdB,difB,mjdD,difD,sigA,sigC,filtername,plot_output_dir):

	""" 
	
	Creates plot of % change in flux over MJD for a single object / filter. Each chip is 
	plotted seperatley. Scatter plot of % changes in flux with the linear fit (weighted by
	the individual measurement errors). 
	
	Parameters:
		objname : string
			name of star 
			
		mjdA,mjdB,mjdC,mjdD: arrays of floats
			array of MJD in each quadrant (A, B , C, D)
			
		difA, difB, difC, difD: arrays of floats
			array of % change in flux in each quadrant
			
		sigA, sigC: arrays of floats
			array of errors associated with % change in flux in quads A and C
			
		filtername: string
			name of filter
			
		plot_output_dir: string
			plot output directory
			
		Outputs:
			plot of % change in flux over MJD with linear fit, in the designated
			plot output directory

	
	"""
	fig = plt.figure()
	fig.subplots_adjust(hspace=0.1)
	
	if len(mjdA) > 0:
		ax1 = fig.add_subplot(2,1,1) #UVIS chip 1, amp A
		ax1.set_xlim([54800,58500])
		ax1.set_ylim((-3.5,3.5))
		ax1.scatter(mjdA, difA, label = 'amp A',c='r',s=80,edgecolor='k',zorder=3)
		ax1.errorbar(mjdA,difA,yerr=sigA,linestyle="None",c='r',markeredgecolor='r',zorder=2,capsize=3)
		ax1.scatter(mjdB,difB,label='amp B',c='k',s=80,marker='x',edgecolor='k',zorder=6)
		ax1.axhline(0,c='k',ls='dashed')
		ax1.legend(loc=1)
		m,b = polyfit(mjdA, difA, 1,w=sigA)
		#print m
		ax1.text(0.01,0.042,'m = ' + str(np.round(np.array(m*365.),3))+' % / MJD',transform=ax1.transAxes)
		stdev=np.std(difA)
		ax1.text(0.01,0.12,r'$\sigma$ = '+str(np.round(np.array(stdev),3)),transform=ax1.transAxes)
		x=np.linspace(mjdA.min(),mjdA.max(),10)
		ax1.plot(x,m*x+b,c='k',lw=3)
		
		
	plt.ylabel(r'$\Delta$ Flux [\%]',fontsize=20)
	ax1.text(0.45,0.9,'UVIS1',transform=ax1.transAxes)
	ax1.text(0,1.0,str(filtername),transform=ax1.transAxes)

	
	if len(mjdC) > 0:
		ax2= fig.add_subplot(2,1,2,sharex=ax1)
		ax2.set_xlim([54800,58500])
		ax2.set_ylim((-3.5,3.5))
		ax2.errorbar(mjdC,difC,yerr=sigC,linestyle="None",c='b',zorder=1,capsize=3)
		ax2.scatter(mjdC, difC, label = 'amp C',c='b',s=80,edgecolor='k',zorder=4)
		ax2.scatter(mjdD,difD,label='amp D',c='k',s=80,marker='x',edgecolor='k',zorder=5)
		
		ax2.axhline(0,c='k',ls='dashed')
		m,b = polyfit(mjdC, difC, 1,w=sigC)
		#print m
		ax2.text(0.01,0.042,'m = ' + str(np.round(np.array(m*365.),3))+' % / MJD',transform=ax2.transAxes)
		stdev=np.std(difC)
		ax2.text(0.01,0.12,r'$\sigma$ = '+str(np.round(np.array(stdev),3)),transform=ax2.transAxes)
		x=np.linspace(mjdC.min(),mjdC.max(),10)
		ax2.plot(x,m*x+b,c='k',lw=3)
		ax2.legend(loc=1)
		
		
		ax2.text(0.45,0.9,'UVIS2',transform=ax2.transAxes)

	fig.suptitle(objname,fontsize=20)
	plt.xlabel(r'MJD',fontsize=20)
	plt.ylabel(r'$\Delta$ Flux [\%]',fontsize=20)
	
	plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
	
	plt.savefig(plot_output_dir+str(objname)+'_'+str(filtername)+'.png',dpi=1000)
	
	plt.close()

def main_fluxDif_MJD(ifile,filtername,objname, plot_output_dir,wrapper):

	""" Main function that calls other functions to make plots and updates stats catalog.
	
		Parameters:
			ifile: string
				Path to photometry catalog.
			filtername: string
				Name of filter
			objname: string
				Name of star
			plot_output_dir: string
				Location where plots / stats catalog will be output
			wrapper: bool
				T if this function is called from wrapper script, F otherwise
			
		
	"""
	
	bad_guys=[]
	#get the data
	data = load_phot_data(ifile) #dataframe
	#select for chip and amp type
	
	dataA=data[(data['amp']=='A')]
	dataC=data[(data['amp']=='C')]
	
	dataB=data[(data['amp']=='B')]
	dataD=data[(data['amp']=='D')]
	

	####fluxes for each amp (chip), then find errors in flux from the magnitude errors####
	fluxA = np.array(dataA['f10'].tolist())/np.array(dataA['expt'].tolist())
	fluxC = np.array(dataC['f10'].tolist())/np.array(dataC['expt'].tolist())
	if len(fluxA)==0:
		mjdA,difA,sigA,divA=[],[],[],1
	else:

		mjdA=np.array(dataA['mjd'].tolist())
		
		#calculate divisor
		divA=generate_divisor(mjdA,fluxA)
		
		#calculate percent difference
		difA=np.array(((fluxA-divA)/divA)*100.)
		
		#flux error
		sigA=np.array(dataA['ferr10'].tolist())/np.array(dataA['expt'].tolist())
		sigA=np.array((((fluxA+sigA)-divA)/divA)*100.)-difA

		fileA=dataA['imagefilename']
		
		for i, file in enumerate(fileA):
			if np.abs(difA[i])>3.5:
				bad_guys.append((file,difA[i]))
				
		mjdA=mjdA[np.abs(difA)<3.5]
		sigA=sigA[np.abs(difA)<3.5]
		difA=difA[np.abs(difA)<3.5]
		
		
	if len(fluxC)==0:
		mjdC,difC,sigC,divC=[],[],[],1
		
	else:
		m10C=np.array(dataC['m10'].tolist())
		mErr_C=np.array(dataC['merr10'].tolist())
		#add mag errors to mag, convert that to flux, subtract central value of flux for flux error
		msig_C=m10C+mErr_C
		zeroPoint_C=m10C+2.5*np.log10(fluxC)
		sigC=np.abs(np.sqrt(10**((2*zeroPoint_C)/2.5)*(1/(2.5**2))*((np.log(10))**2)*(10**((-2.0*m10C)/2.5))*(mErr_C**2)))

		mjdC=np.array(dataC['mjd'].tolist())
		
		#calculate divisor
		divC=generate_divisor(mjdC,fluxC)
		#divC=divA #normalize to chipA
		
		#calculate percent difference
		difC=np.array(((fluxC-divC)/divC)*100.)
		
		#flux error
		sigC=np.array(dataC['ferr10'].tolist())/np.array(dataC['expt'].tolist())
		sigC=np.array((((fluxC+sigC)-divC)/divC)*100.)-difC

		fileC =dataC['imagefilename']
		
		for i, file in enumerate(fileC):
			if np.abs(difC[i])>3.5:
				bad_guys.append((file,difC[i]))

				
		mjdC=mjdC[np.abs(difC)<3.5]
		difC=difC[np.abs(difC)<3.5]
		sigC=sigC[np.abs(difC)<3.5]

	#add in data from amps B and D just for plotting not to fit the slopes to
	fluxB = np.array(dataB['f10'].tolist())/np.array(dataB['expt'].tolist())
	fluxD= np.array(dataD['f10'].tolist())/np.array(dataD['expt'].tolist())
	mjdB,mjdD =dataB['mjd'].tolist(),dataD['mjd'].tolist()

	
	#calculate percent difference
	if len(fluxB)>0:
		
		
		divB = np.median(fluxB[0:3])
		
		difB=np.array(((fluxB-divB)/divB)*100.)	
		
	else:
		difB=[]
		
	if len(fluxD)>0:
		print 'D'
		divD = np.median(fluxD[0:3])
		difD=np.array(((fluxD-divD)/divD)*100.)	
	
	else:
		difD=[]
	
	
	#find horrible outliers....throw out and examine
	
	fileB,fileD=dataB['imagefilename'],dataD['imagefilename']

	plotting(objname,mjdA,difA,mjdC,difC,mjdB,difB,mjdD,difD,sigA,sigC,filtername,plot_output_dir)
			
	
	if wrapper: #if this function is being called in the plotting wrapper over all filters,
				#make the stats.dat file
		append_to_stats_file(objname,filtername,plot_output_dir,mjdA,difA,sigA,mjdC,difC,sigC)
		print 'updating stats.dat file in output directory'

if __name__ == '__main__':

	###############
	#set these!!!!
	ifile='../data/data/GRW70/F606W/flt_cleans/F606W_photcat.dat'
	filtername='F606W'
	objname='GRW70'
	plot_output_dir=''
	###############
	
	
	main_fluxDif_MJD(ifile, filtername, objname,plot_output_dir,wrapper=False)