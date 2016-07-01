#! /usr/bin/env python
# Colm Coughlan 20.11.2015
# Dublin Institute for Advanced Studies

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

myfontsize = 15
plt.rcParams.update({'font.size': myfontsize})

parser = argparse.ArgumentParser(description='Colm Coughlan. Dublin Institute for Advanced Studies. Make offset plots from make_catalog.py output.')
parser.add_argument('stem1', type=str, help='Results stem one.')
parser.add_argument('stem2', type=str, help='Results stem two.')
parser.add_argument('stem3', type=str, help='Results stem three.')
parser.add_argument('output', type=str, help='Stem for output files.')

args = parser.parse_args()

stems = [args.stem1, args.stem2, args.stem3]

# NVSS offset scatter plots for each field and frequency

frame_list=[[],[]]

for i in stems:
	for j in range(2):
		if j==0:
			freq = '325 MHz'
		else:
			freq = '610 MHz'
		name = i.split('/')[-1]
		if name == 'DGTau':
			pr_name = 'DG Tau'
		if name == 'TTau':
			pr_name = 'T Tau'
		if name == 'L1551':
			pr_name = 'L1551 N1'
			
		df = pd.read_csv(i+'_freq'+str(j+1)+'.nvss_offset.csv',skiprows=0,delimiter=',', engine='python')
		frame_list[j].append(df)
		
		median_ra, sigma_ra = np.median(3600.0*df['NVSS_RA_offset'].values), np.std(3600.0*df['NVSS_RA_offset'].values)
		median_dec, sigma_dec = np.median(3600.0*df['NVSS_DEC_offset'].values), np.std(3600.0*df['NVSS_DEC_offset'].values)
		plt.scatter(3600.0*df['NVSS_RA_offset'].values, 3600.0*df['NVSS_DEC_offset'].values)

		plt.xlabel(r'$\mathrm{RA\,Offset\,(arcsec)}$')
		plt.ylabel(r'$\mathrm{DEC\,Offset\,(arcsec)}$')
		plt.title(r'$\mathrm{'+pr_name+'\,'+freq+'}$')
		print('Running '+pr_name+' '+freq+' with '+str(len(df))+' points.')
		print('\t Median RA offset = %.2f +/- %.2f arcsec.'%(median_ra,sigma_ra))
		print('\t Median DEC offset = %.2f +/- %.2f arcsec.'%(median_dec,sigma_dec))
		plt.savefig(args.output+name+'_freq'+str(j+1)+'.eps')
		plt.clf()
		
for j in range(2):
	df = pd.concat(frame_list[j])
	print 'Printing '+str(len(df))+' sources at freq '+str(j+1)
	if j==0:
		freq = '325 MHz'
	else:
		freq = '610 MHz'
	median_ra, sigma_ra = np.median(3600.0*df['NVSS_RA_offset'].values), np.std(3600.0*df['NVSS_RA_offset'].values)
	median_dec, sigma_dec = np.median(3600.0*df['NVSS_DEC_offset'].values), np.std(3600.0*df['NVSS_DEC_offset'].values)
	plt.scatter(3600.0*df['NVSS_RA_offset'].values, 3600.0*df['NVSS_DEC_offset'].values)
	plt.ylim(-6,6)
	plt.xlim(-6,6)
	plt.xlabel(r'$\mathrm{RA\,Offset\,(arcsec)}$')
	plt.ylabel(r'$\mathrm{DEC\,Offset\,(arcsec)}$')
	print('Running entire '+freq+' survey with '+str(len(df))+' points.')
	print('\t Median RA offset = %.2f +/- %.2f arcsec.'%(median_ra,sigma_ra))
	print('\t Median DEC offset = %.2f +/- %.2f arcsec.'%(median_dec,sigma_dec))
	plt.savefig(args.output+'combined_freq'+str(j+1)+'.eps')
	plt.clf()