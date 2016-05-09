# Colm Coughlan
# DIAS

import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import sys

t = 1.0e4		# electron temperature
theta2 = 17.4	#arcsec	squared
theta = np.sqrt(theta2)	# geometric mean of source
flux = 1.88	# model flux in mJy at turnover frequency
d = 0.14	# distance in kpc


# optical depth

def tau_func(nu, t, em,alpha):
#	return( (3.014e-2) * (t**-1.5) * np.power(nu, alpha) * ( np.log(np.divide(4.955e-2,nu)) +1.5*np.log(t) ) * em )
	return( (8.235e-2) * np.power(t,-1.35) * np.power(nu,alpha) * em)	# altenhoff 1960 approximation
	
# flux density from temperature and optical depth (normalising to set flux scale to Jy)

def sed_func(nu, em,alpha1, norm1, alpha2, norm2):
	tau = tau_func(nu, t, em,alpha1)
	flux_density1 = norm1 * (1.0 - np.exp(-tau)) * np.power(nu/0.149, 2.0) # get flux density in Jy
	flux_density2 = norm2 * np.power((nu/100.0),alpha2)
	flux_density = flux_density1 + flux_density2
	return( flux_density )
	
# flux density from temperature and optical depth (normalising to set flux scale to Jy)
# should be identical to function above, except returns the two models separately

def sed_func_plotonly(nu, em,alpha1, norm1, alpha2, norm2):
	tau = tau_func(nu, t, em,alpha1)
	flux_density1 = norm1 * (1.0 - np.exp(-tau)) * np.power(nu/0.149, 2.0) # get flux density in Jy
	flux_density2 = norm2 * np.power((nu/100.0),alpha2)
	return( flux_density1, flux_density2 )

if(len(sys.argv)!=2):
	print("\tError: Takes 1 argument.")
	print("\tUseage: sedfit <datafile")
	sys.exit()
else:
	inputname = str(sys.argv[1])


# read data
df = pd.read_csv(inputname,delimiter='\t')

# perform fit, identify values and errors

res, cov = curve_fit(sed_func, df['Freq'].values, df['Flux'].values, sigma = df['Error'].values, maxfev = 999999,p0=[56397.0,-1.8,1.0,2.0,10.0]) # guesses that work for 900 GHz norm: p0=[56397.0,-1.8,1.0,10.0,1.0]
em, alpha1, norm1, alpha2, norm2 = res[0], res[1], res[2], res[3], res[4]
cov = np.sqrt(cov)
em_err, alpha1_err, norm1_err, alpha2_err, norm2_err = cov[0][0], cov[1][1], cov[2][2], cov[3][3], cov[4][4]

# calculate parameters from EM

v0 = np.power(8.235e-2, -1.0/alpha1)*np.power(t,1.35/alpha1)*np.power(em, -1.0/alpha1)
#v0 = 0.3045 * np.power(t,-0.643) * np.power(em, 0.476)	# this is the canonical formula (it is inconsistent to use this with a value of EM derived from a difference alpha)
mion = (3.4e-5) * np.power((t/1.0e4),0.175) * np.power(v0, 0.05) * np.power(flux,0.5) * np.power(d, 2.5) * np.power(theta,1.5)
ne = (7.2e3) * np.power((t/1.0e4),0.175) * np.power(v0, 0.05) * np.power(flux,0.5) * np.power(d, -0.5) * np.power(theta,-1.5)

print 'em = '+str(em)+' +/- '+str(em_err)
print 'Alpha1 = '+str(alpha1)+' +/- '+str(alpha1_err)
print 'norm1 = '+str(norm1*1000.0)+' +/- '+str(norm1_err*1000.0) + ' mJy'
print 'Alpha2 = '+str(alpha2)+' +/- '+str(alpha2_err)
print 'norm2 = '+str(norm2*1000.0)+' +/- '+str(norm2_err*1000.0) + ' mJy'
print('Derived parameters:')
print('\t v0(GHz)   =   '+str(v0))
print('\t mion(solar masses)   =   '+str(mion))
print('\t n_e(cm-3)   =   '+str(ne))

# Generate plot of derived model

model_freqs = np.logspace(-3,3,num=1000, base=10.0)
model_fluxes1, model_fluxes2 = sed_func_plotonly(model_freqs, em, alpha1, norm1, alpha2, norm2)
model_taus = tau_func(model_freqs, t, em, alpha1)

plt.figure(0)
plt.loglog(model_freqs, model_fluxes1,'k--',label=r'$\mathrm{Free-Free}$')
plt.loglog(model_freqs, model_fluxes2,'k:',label=r'$\mathrm{Power\,law}$')
plt.loglog(model_freqs, model_fluxes1 + model_fluxes2,'k-',label=r'$\mathrm{Combined}$')
plt.errorbar(df['Freq'].values, df['Flux'].values, yerr = df['Error'].values, linestyle='None', fmt='ko', ecolor='black')
plt.xlabel(r'$\nu\,\mathrm{(GHz)}$')
plt.ylabel(r'$S_{\nu}\,\mathrm{(Jy)}$')
plt.title(r'$\mathrm{T\,Tau\,Spectral\,Energy\,Distribution}$')
plt.legend(loc=2)

plt.figure(1)
plt.loglog(model_freqs, model_taus)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Tau')
plt.title('T Tau opacity vs. frequency')

plt.show()
