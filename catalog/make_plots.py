#! /usr/bin/env python
# Colm Coughlan 20.11.2015
# Dublin Institute for Advanced Studies


'''
L1551 IRS 5 field at 610 MHz
T Tau field
DG Tau field
'''

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import scipy.stats

myfontsize = 15
plt.rcParams.update({'font.size': myfontsize})

def sensitivity_limit(flux_325, sens_limit_upper, freq_high):
	freq_low = 322.666
	return(np.log(sens_limit_upper/flux_325)/np.log(freq_high/freq_low))

parser = argparse.ArgumentParser(description='Colm Coughlan. Dublin Institute for Advanced Studies. Make a survey plots from make_catalog.py output.')
parser.add_argument('stem1', type=str, help='Results stem one.')
parser.add_argument('stem2', type=str, help='Results stem two.')
parser.add_argument('stem3', type=str, help='Results stem three.')
parser.add_argument('output', type=str, help='Stem for output files.')

args = parser.parse_args()

stems = [args.stem1, args.stem2, args.stem3]
name = 'fullcat'

#########################################################################################################
#
# Residuals histogram
#
#########################################################################################################

# Freq 1
df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_freq1.fluxes.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
	
print('Making RMS histogram with '+str(len(df))+' datapoints.')
#data = np.log10(df['Resid_Isl_rms'].values)
data = df['E_Peak_flux'].values*(10**3)
median_rms1 = np.median(data)/(10**6)
print('Median peak uncertainty at 323MHz = '+str(median_rms1*1000.0)+' mJy/Beam')
data = np.log10(data)
mean, sigma = np.mean(data), np.std(data)
#data = data[np.where(data<(mean+1.0*sigma))]
xmin = np.min(data)
xmax = np.max(data)
#xmax = 2.0*np.median(data)
step = 20

hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
hist = hist/float(len(data))
cdf = np.cumsum(hist)
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
plt.plot(bin_edges[:-1],cdf)
plt.xlabel(r'$\mathrm{323\,MHz\,Peak\,flux\,uncertainty\,(log_{10}\mu JyBeam^{-1})}$')
plt.ylabel(r'$\mathrm{Fraction of sources}$')
plt.savefig(args.output+name+'_err_freq1.eps')
plt.clf()

# Freq 2
df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_freq2.fluxes.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
	
print('Making RMS histogram with '+str(len(df))+' datapoints.')
#data = np.log10(df['Resid_Isl_rms'].values)
data = df['E_Peak_flux'].values*(10**3)
median_rms2 = np.median(data)/(10**6)
print('Median peak uncertainty at 608MHz = '+str(median_rms2*1000.0)+' mJy/Beam')

data = np.log10(data)
mean, sigma = np.mean(data), np.std(data)
#data = data[np.where(data<(mean+1.0*sigma))]
xmin = np.min(data)
xmax = np.max(data)
step = 20

hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
hist = hist/float(len(data))
cdf = np.cumsum(hist)
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
plt.plot(bin_edges[:-1],cdf)
plt.xlabel(r'$\mathrm{608\,MHz\,Peak\,flux\,uncertainty\,(log_{10}\mu JyBeam^{-1})}$')
plt.ylabel(r'$\mathrm{Fraction of sources}$')
plt.savefig(args.output+name+'_err_freq2.eps')
plt.clf()

#########################################################################################################
#
# Flux histogram
#
#########################################################################################################
# freq1

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_freq1.fluxes.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
print('Making flux histograms with '+str(len(df))+' datapoints.')

data = df['Total_flux'].values/1000.0
data = np.log10(data)
xmin = np.min(data)
xmax = np.max(data)
step = 30
hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
plt.xlabel(r'$\mathrm{323\,MHz\,log_{10}(Integrated\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{N}$')
plt.savefig(args.output+name+'_total_flux_freq1.eps')
plt.clf()

data = df['Peak_flux'].values/1000.0
data = np.log10(data)
xmin = np.min(data)
xmax = np.max(data)
step = 30
hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
plt.xlabel(r'$\mathrm{323\,MHz\,log_{10}(Peak\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{N}$')
plt.savefig(args.output+name+'_peak_flux_freq1.eps')
plt.clf()

# freq2
df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_freq2.fluxes.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
print('Making flux histograms with '+str(len(df))+' datapoints.')

data = df['Total_flux'].values/1000.0
data = np.log10(data)
xmin = np.min(data)
xmax = np.max(data)
step = 30
hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
plt.xlabel(r'$\mathrm{608\,MHz\,log_{10}(Integrated\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{N}$')
plt.savefig(args.output+name+'_total_flux_freq2.eps')
plt.clf()

data = df['Peak_flux'].values/1000.0
data = np.log10(data)
xmin = np.min(data)
xmax = np.max(data)
step = 30
hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
plt.xlabel(r'$\mathrm{608\,MHz\,log_{10}(Peak\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{N}$')
plt.savefig(args.output+name+'_peak_flux_freq2.eps')
plt.clf()

# SI histogram

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'.spx.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
	
df = pd.concat(frame_list)

data = df['SI'].values
median, sigma = np.median(data), np.std(data)
xmin = np.min(data)
xmax = np.max(data)
step = 30

hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
print('Median alpha_GMRT = '+str(median)+' +/- '+str(sigma))
plt.xlabel(r'$\mathrm{\alpha_{GMRT}}$')
plt.ylabel(r'$\mathrm{N}$')
plt.savefig(args.output+name+'_alpha_gmrt.eps')
plt.clf()

# NVSS SI histogram

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'.nvss_spx.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
	
df = pd.concat(frame_list)

data = df['NVSS_SI'].values
#weights = df['NVSS_E_SI'].values
#weights = 1.0/np.power(weights,2.0)
median, sigma = np.median(data), np.std(data)
step = 30

data_both = df.loc[df['NVSS_SI_FP']=='B','NVSS_SI'].values
xmin1 = np.min(data_both)
xmax1 = np.max(data_both)

data_single = df.loc[df['NVSS_SI_FP']!='B','NVSS_SI'].values
xmin2 = np.min(data_single)
xmax2 = np.max(data_single)

xmin = np.min([xmin1, xmin2])
xmax = np.max([xmax1, xmax2])
hist1, bin_edges = np.histogram(data_both, bins=np.linspace(xmin, xmax, step))
hist2, bin_edges = np.histogram(data_single, bins=np.linspace(xmin, xmax, step))

plt.bar(bin_edges[:-1], hist1, width = np.abs((xmax-xmin)/step),color='red')
plt.bar(bin_edges[:-1], hist2, width = np.abs((xmax-xmin)/step), bottom = hist1)

#hist, bin_edges = np.histogram(data, bins=np.linspace(xmin, xmax, step))
#plt.bar(bin_edges[:-1], hist, width = np.abs((xmax-xmin)/step))
plt.xlim(min(bin_edges), max(bin_edges))
print('Median alpha_GMRT_NVSS = '+str(median)+' +\ '+str(sigma))
plt.xlabel(r'$\mathrm{\alpha_{GMRT-NVSS}}$')
plt.ylabel(r'$\mathrm{N}$')
plt.legend([r'$\mathrm{3\,point}$',r'$\mathrm{2\,point}$'])
plt.savefig(args.output+name+'_alpha_nvss.eps')
plt.clf()


# Source Position plot (combined frequencies. set xaxis)

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_positions.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)

df = pd.concat(frame_list)

plt.scatter(df['RA'].values, df['DEC'].values)
plt.gca().invert_xaxis()

plt.xlabel(r'$\mathrm{RA\,(deg)}$')
plt.ylabel(r'$\mathrm{DEC\,(deg)}$')
plt.ylim(16,28)
plt.xlim(61,73)
plt.savefig(args.output+name+'_positions.eps')
plt.clf()

# spectral index vs flux plots (GMRT and NVSS)

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
df.reset_index(inplace=True, drop=True)


flux = df.loc[np.all([df['Match']==True,df['Freq'].values<4E8],axis=0),'Total_flux'].values
spindx = df.loc[np.all([df['Match']==True,df['Freq'].values<4E8],axis=0),'SI']
min,max = np.min(flux),np.max(flux)
smin, smax = sensitivity_limit(min, 0.2, 607.664), sensitivity_limit(max, 0.2, 607.664)# Flux in mJy
min,max = np.log(min),np.log(max)

plt.xlabel(r'$\mathrm{323\,MHz\,log_{10}(Integrated\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{\alpha}$')
plt.ylim(-6,3)
plt.scatter(np.log(flux), spindx)
plt.plot([min, max], [smin,smax],'--r')
plt.savefig(args.output+name+'_alpha_gmrt_evo.eps')
plt.clf()

flux = df.loc[np.all([df['Match']==True,df['Freq'].values<4E8],axis=0),'Total_flux']
spindx = df.loc[np.all([df['Match']==True,df['Freq'].values<4E8],axis=0),'E_SI']

plt.xlabel(r'$\mathrm{323\,MHz\,log_{10}(Integrated\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{\Delta \alpha}$')
plt.scatter(np.log(flux), spindx)
plt.savefig(args.output+name+'_err_alpha_gmrt_evo.eps')
plt.clf()

flux = df.loc[np.all([df['NVSS']==True,df['Freq'].values<4E8],axis=0),'Total_flux']
spindx = df.loc[np.all([df['NVSS']==True,df['Freq'].values<4E8],axis=0),'NVSS_SI']
min,max = np.min(flux),np.max(flux)
smin, smax = sensitivity_limit(min, 1.5, 1400), sensitivity_limit(max, 1.5, 1400)# Flux in mJy
min,max = np.log(min),np.log(max)

plt.xlabel(r'$\mathrm{323\,MHz\,log_{10}(Integrated\,Flux\,(Jy))}$')
plt.ylabel(r'$\mathrm{\alpha_{GMRT-NVSS}}$')
plt.ylim(-2,1)
plt.scatter(np.log(flux), spindx)
plt.plot([min, max], [smin,smax],'--r')
plt.savefig(args.output+name+'_alpha_nvss_evo.eps')
plt.clf()


# Differential source counts - combined 3 fields, separate frequencies


pbeam1 = 0.8*81.0/60.0	# in degrees, low freq
pbeam2 = 0.8*43.0/60.0
area1 = np.power(pbeam1*(np.pi/180.0),2.0) * np.pi / (4.0*np.log(2))	# areas in steradiens
area2 = np.power(pbeam2*(np.pi/180.0),2.0) * np.pi / (4.0*np.log(2))

# load in data, plus peak fluxes for scaling source counts of weaker sources

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_freq1.fluxes.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
data = df['Total_flux'].values/1000.0
peakdata = df['Peak_flux'].values/1000.0

df = pd.DataFrame()
frame_list=[]
for i in range(len(stems)):
	tdf = pd.read_csv(stems[i]+'_freq2.fluxes.csv',skiprows=0,delimiter=',', engine='python')
	frame_list.append(tdf)
df = pd.concat(frame_list)
data2 = df['Total_flux'].values/1000.0
peakdata2 = df['Peak_flux'].values/1000.0


xmin = np.min(data)
xmax = 0.5#np.max(data)
xmin2 = np.min(data2)
xmax2 = 0.5*np.max(data2)
step = 20

counts, bin_edges, bin_num = scipy.stats.binned_statistic(x=data,values=data, statistic='count', bins=np.logspace(np.log10(xmin), np.log10(xmax), num = step, base = 10))
means, bin_edges, bin_num = scipy.stats.binned_statistic(x=data,values=data, statistic='mean', bins=np.logspace(np.log10(xmin), np.log10(xmax), num = step, base = 10))
peakmeans, bin_edges, bin_num = scipy.stats.binned_statistic(x=data,values=peakdata, statistic='mean', bins=np.logspace(np.log10(xmin), np.log10(xmax), num = step, base = 10))
stds, bin_edges, bin_num = scipy.stats.binned_statistic(x=data,values=data, statistic=np.std, bins=np.logspace(np.log10(xmin), np.log10(xmax), num = step, base = 10))
counts2, bin_edges2, bin_num2 = scipy.stats.binned_statistic(x=data2,values=data2, statistic='count', bins=np.logspace(np.log10(xmin2), np.log10(xmax2), num = step, base = 10))
means2, bin_edges2, bin_num2 = scipy.stats.binned_statistic(x=data2,values=data2, statistic='mean', bins=np.logspace(np.log10(xmin2), np.log10(xmax2), num = step, base = 10))
peakmeans2, bin_edges2, bin_num2 = scipy.stats.binned_statistic(x=data2,values=peakdata2, statistic='mean', bins=np.logspace(np.log10(xmin2), np.log10(xmax2), num = step, base = 10))
stds2, bin_edges2, bin_num2 = scipy.stats.binned_statistic(x=data2,values=data2, statistic=np.std, bins=np.logspace(np.log10(xmin2), np.log10(xmax2), num = step, base = 10))


# Correct low significance counts for sensitivity effects (scale by median rms / central rms)

central_rms1 = 0.000127
central_rms2 = 0.000067

corr1 = median_rms1/central_rms1
corr2 = median_rms2/central_rms2

print('Correcting 323 by '+str(corr1))
print('Correcting 608 by '+str(corr2))
for i in range(len(peakmeans)):
	if(peakmeans[i] < (7.0 * central_rms1)):
		counts[i] = counts[i] * corr1 * ( (6.0 * central_rms1)/peakmeans[i] )
		print 'Correcting 323 MHz bin'
		
for i in range(len(peakmeans2)):
	if(peakmeans2[i] < (7.0 * central_rms2)):
		counts2[i] = counts2[i] * corr2 * ( (6.0 * central_rms2)/peakmeans2[i] )
		print 'Correcting 608 MHz bin'

dstep=np.zeros((len(bin_edges)-1,1))
dstep2=np.zeros((len(bin_edges)-1,1))
for i in range(step-1):
	dstep[i] = bin_edges[i+1] - bin_edges[i]
	dstep2[i] = bin_edges2[i+1] - bin_edges2[i]

diff_sc = np.multiply( np.power(means, 2.5), (np.divide(counts, dstep.T) /area1))[0]
diff_sc2 = np.multiply( np.power(means2, 2.5), (np.divide(counts2, dstep2.T) /area2))[0]

#error = diff_sc*np.divide(stds,means)
#error2 = diff_sc2*np.divide(stds2,means2)

error = np.divide(diff_sc, np.sqrt(counts))
error2 = np.divide(diff_sc2, np.sqrt(counts2))

fig = plt.errorbar(bin_edges[:-1], diff_sc, yerr = error, linestyle='None',label='323 MHz')
fig = plt.errorbar(bin_edges2[:-1], diff_sc2, yerr = error2, linestyle='None',label='608 MHz')
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
plt.xlabel(r'$\mathrm{Flux\,(Jy)}$')
plt.ylabel(r'$\mathrm{Euclidian\,normalised\,differential\,source\,count}$')
plt.legend(loc=2)
plt.savefig(args.output+name+'_dsc.eps')
plt.clf()