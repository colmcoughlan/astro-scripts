import numpy as np
from scipy.optimize import curve_fit
from scipy import constants
import pandas as pd
import matplotlib.pyplot as plt
import sys

Boltzmann = 1.38064852e-23	# scipy constants v.0.9 does not have this, but constants.Boltzmann should work in modern versions
t = 1.0e4
theta2 = 17.4	#arcsec	squared



# should not need to change anything below here



#theta2 = theta2 * np.power( (np.pi/(3600.0*180.0)) ,2)
omega = np.pi * theta2 / (4.0 * np.log(2))

print 't = '+str(t)
print 'omega = '+str(omega)

def tau_func(nu, t, em,alpha):
#	return( (3.014e-2) * (t**-1.5) * np.power(nu, alpha) * ( np.log(np.divide(4.955e-2,nu)) +1.5*np.log(t) ) * em )
	return( (8.235e-2)*np.power(t,-1.35)*np.power(nu,alpha)*em)

def sed_func(nu, em,alpha1, norm1, alpha2, norm2):
	# from Anna Scaife's pebbles paper, page 9
	# this tau needs nu in GHz, t in K, EM in cm-6 pc

	tau = tau_func(nu, t, em,alpha1)
	wavelength = (constants.c/(1.0e9*nu))*100.0	# get lambda in centimetres from freq in gigahertz
	flux_density1 = norm1 * (1.0 - np.exp(-tau)) * theta2 * t / (wavelength*wavelength)	# get flux density in Jy/Beam. This formula needs theta in arcsec, wavelength in cm
#	flux_density1 = norm1 * np.power((nu/0.149),alpha1)
	flux_density2 = norm2 * np.power((nu/900.0),alpha2)
	flux_density = flux_density1 + flux_density2
#	flux_density = (2.0 * Boltzmann * t /(wavelength*wavelength)) * (1.0 - np.exp(-tau)) * omega
	return( flux_density )

def sed_func_plotonly(nu, em,alpha1, norm1, alpha2, norm2):
	# from Anna Scaife's pebbles paper, page 9
	# this tau needs nu in GHz, t in K, EM in cm-6 pc

	tau = tau_func(nu, t, em,alpha1)
	wavelength = (constants.c/(1.0e9*nu))*100.0	# get lambda in centimetres from freq in gigahertz
	flux_density1 = norm1 * (1.0 - np.exp(-tau)) * theta2 * t / (wavelength*wavelength)	# get flux density in Jy/Beam
#	flux_density1 = norm1 * np.power((nu/0.149),alpha1)
	flux_density2 = norm2 * np.power((nu/900.0),alpha2)
	return( flux_density1, flux_density2 )

if(len(sys.argv)!=2):
	print("\tError: Takes 1 argument.")
	print("\tUseage: sedfit <datafile")
	sys.exit()
else:
	inputname = str(sys.argv[1])


df = pd.read_csv(inputname,delimiter='\t')

res, cov = curve_fit(sed_func, df['Freq'].values, df['Flux'].values, sigma = df['Error'].values, maxfev = 999999,p0=[56397.0,-1.8,0.001,10.0,2.0])
em, alpha1, norm1, alpha2, norm2 = res[0], res[1], res[2], res[3], res[4]
cov = np.sqrt(cov)
em_err, alpha1_err, norm1_err, alpha2_err, norm2_err = cov[0][0], cov[1][1], cov[2][2], cov[3][3], cov[4][4]


print 'em = '+str(em)+' +/- '+str(em_err)
print 'Alpha1 = '+str(alpha1)+' +/- '+str(alpha1_err)
print 'norm1 = '+str(norm1)+' +/- '+str(norm1_err)
print 'Alpha2 = '+str(alpha2)+' +/- '+str(alpha2_err)
print 'norm2 = '+str(norm2)+' +/- '+str(norm2_err)

model_freqs = np.logspace(-3,3,num=100, base=10.0)
model_fluxes1, model_fluxes2 = sed_func_plotonly(model_freqs, em, alpha1, norm1, alpha2, norm2)
model_taus = tau_func(model_freqs, t, em, alpha1)

plt.figure(0)
plt.loglog(model_freqs, model_fluxes1 + model_fluxes2,label='Combined model')
plt.loglog(model_freqs, model_fluxes1,label='Free-free model')
plt.loglog(model_freqs, model_fluxes2,label='Power law model')
plt.errorbar(df['Freq'].values, df['Flux'].values, yerr = df['Error'].values, linestyle='None', fmt='o',label='Data')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Flux (Jy)')
plt.title('T Tau spectral energy distribution')
plt.legend(loc=4)

plt.figure(1)
plt.loglog(model_freqs, model_taus)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Tau')
plt.title('T Tau opacity vs. frequency')

plt.show()
