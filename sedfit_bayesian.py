# Colm Coughlan
# DIAS

import numpy as np
import pymc as mc
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import sys
import corner

myfontsize = 15
plt.rcParams.update({'font.size': myfontsize})

theta_p1 = 7.6#8.2#7.6	# arcsec (deconvolved)
theta_p2 = 1.4#2.4#1.4
theta_p1_err = 2.3
theta_p2_err = 1.03
theta2 = theta_p1 * theta_p2	#arcsec	squared
omega = theta2 * np.pi / (4.0*np.log(2))	# solid angle
theta = np.sqrt(theta2)	# geometric mean of source
c_omega = 7.21586e-4 * omega	# in mJy
d = 0.14	# distance in kpc
d_err = 0.01	# ten parsec error from Kenyon 1994
flat_t = 1.0e4	# flat temperature that can be used for fitting (K)

# optical depth

def tau_func(nu, t, em,alpha):
#	return( (3.014e-2) * (t**-1.5) * np.power(nu, alpha) * ( np.log(np.divide(4.955e-2,nu)) +1.5*np.log(t) ) * em )
	return( (8.235e-2) * np.power(t,-1.35) * np.power(nu,alpha) * em)	# altenhoff 1960 approximation
	
# flux density from temperature and optical depth (normalising to set flux scale to mJy)

def sed_func(nu, em,alpha1, t, alpha2, norm2):
	tau = tau_func(nu, t, em,alpha1)
	flux_density1 = t * c_omega * (1.0 - np.exp(-tau)) * np.power(nu, 2.0) # get flux density in mJy
	flux_density2 = norm2 * np.power((nu/100.0),alpha2)	# normalised power law
	total_flux_density = flux_density1 + flux_density2
	
	if(not plot_mode):
		return( total_flux_density )
	else:
		return( flux_density1, flux_density2 )

if(len(sys.argv)!=2):
	print("\tError: Takes 1 argument.")
	print("\tUseage: sedfit <datafile")
	sys.exit()
else:
	inputname = str(sys.argv[1])


# read data, convert to mJy
df = pd.read_csv(inputname,delimiter='\t')
df['Flux'] = np.multiply((df['Flux'].values),1000.0)
df['Error'] = np.multiply((df['Error'].values),1000.0)

# perform fit, identify values and errors

plot_mode = False
res, cov = curve_fit(sed_func, df['Freq'].values, df['Flux'].values, sigma = df['Error'].values, maxfev = 999999,p0=[56397.0,-1.8,1.0e4,2.0,10.0])
em, alpha1, t, alpha2, norm2 = res[0], res[1], res[2], res[3], res[4]
cov = np.sqrt(cov)
em_err, alpha1_err, t_err, alpha2_err, norm2_err = cov[0][0], cov[1][1], cov[2][2], cov[3][3], cov[4][4]

# calculate parameters from EM

#flux = 1.88	# model flux in mJy at turnover frequency
tau = 1.0
v0 = np.power(8.235e-2, -1.0/alpha1)*np.power(t,1.35/alpha1)*np.power(em, -1.0/alpha1)
#v0 = 0.3045 * np.power(t,-0.643) * np.power(em, 0.476)	# this is the canonical formula (it is inconsistent to use this with a value of EM derived from a difference alpha)
#mion = (3.4e-5) * np.power((t/1.0e4),0.175) * np.power(v0, 0.05) * np.power(flux,0.5) * np.power(d, 2.5) * np.power(theta,1.5)
#ne = (7.2e3) * np.power((t/1.0e4),0.175) * np.power(v0, 0.05) * np.power(flux,0.5) * np.power(d, -0.5) * np.power(theta,-1.5)
ne = np.sqrt(em / (d*theta*(7.1e-3)))	# derive ne from em and theta
mion = (1.47500387e-9) * ne * np.power(theta*d,3.0)

# find error in these parameters

u = np.power(v0, alpha1)
omega_err = omega * np.sqrt(np.power(theta_p1_err/theta_p1,2.0) + np.power(theta_p2_err/theta_p2,2.0))
theta_err = np.power(2.0*theta_p1*theta_p2,-0.5)*np.sqrt(np.power(theta_p2*theta_p1_err,2.0) + np.power(theta_p1*theta_p2_err,2.0))
tau_err = tau * np.sqrt( np.power(em_err/em,2.0) + np.power( (alpha1*alpha1_err)/v0 ,2.0) + np.power( (1.35*t_err)/t,2.0) )
#v0_err = np.sqrt( np.power(np.power(u,(1.0/alpha1)-1.0)/alpha1,2.0) * ( np.power( (alpha1_err * np.log(u))/(alpha1*u), 2.0) + (np.power(v0,2.0*alpha1)*( np.power(em_err/em,2.0) + np.power((1.35*t_err)/t,2.0) ) ) ))
v0_err = (v0/alpha1) * np.sqrt(np.power( (1.35*t_err)/t,2.0) + np.power(em_err/em,2.0) + np.power((np.log(-1.0/alpha1)*alpha1_err)/alpha1,2.0))
#flux_err = flux * (np.power((tau*tau_err)/(np.exp(tau)-1.0),2.0) + np.power(t_err/t,2.0) + np.power(omega_err/omega,2.0))
#ne_err = ne * np.sqrt( np.power( (0.175*t_err)/t,2.0) + np.power( (0.5*flux_err)/flux,2.0) + np.power( (-0.5*d_err)/d,2.0) + np.power( (-1.5*theta_err)/theta,2.0))
ne_err = ne * np.sqrt( np.power( (-0.5*d_err)/d,2.0) + np.power( (-0.5*theta_err)/theta,2.0) + np.power( (0.5*em_err)/em,2.0) )
#mion_err = mion * np.sqrt( np.power( (0.175*t_err)/t,2.0) + np.power( (0.5*flux_err)/flux,2.0) + np.power( (2.5*d_err)/d,2.0) + np.power( (1.5*theta_err)/theta,2.0))
mion_err = mion * np.sqrt( np.power(ne_err/ne, 2.0) + np.power((3.0*theta_err)/theta, 2.0) + np.power((3.0*d_err)/d, 2.0) )

print('\n\n')
print('Fitted parameters:')
print '\t Emission measure = '+str(em)+' +/- '+str(em_err)+' pc cm^-6'
print '\t Alpha (free-free) = '+str(alpha1)+' +/- '+str(alpha1_err)
print '\t Temperature = '+str(t)+' +/- '+str(t_err) + ' K'
print '\t Alpha (power law) = '+str(alpha2)+' +/- '+str(alpha2_err)
print '\t K (power law) = '+str(norm2)+' +/- '+str(norm2_err) + ' mJy'
print('\n')
print('Derived parameters:')
print('\t v0(GHz)   =   '+str(v0)+' +/- '+str(v0_err))
print('\t mion(solar masses)   =   '+str(mion)+' +/- '+str(mion_err))
print('\t n_e(cm-3)   =   '+str(ne)+' +/- '+str(ne_err))
print('\n\n')

# Generate plot of derived model

model_freqs = np.logspace(-3,3,num=1000, base=10.0)
plot_mode = True
model_fluxes1, model_fluxes2 = sed_func(model_freqs, em, alpha1, t, alpha2, norm2)
model_taus = tau_func(model_freqs, t, em, alpha1)

plt.loglog(model_freqs, model_fluxes1/1000.0,'k--',label=r'$\mathrm{Free-Free}$')
plt.loglog(model_freqs, model_fluxes2/1000.0,'k:',label=r'$\mathrm{Power\,law}$')
plt.loglog(model_freqs, model_fluxes1/1000.0 + model_fluxes2/1000.0,'k-',label=r'$\mathrm{Combined}$')
plt.errorbar(df['Freq'].values, df['Flux'].values/1000.0, yerr = df['Error'].values/1000.0, linestyle='None', fmt='ko', ecolor='black')
plt.xlabel(r'$\nu\,\mathrm{(GHz)}$')
plt.ylabel(r'$S_{\nu}\,\mathrm{(Jy)}$')
plt.title(r'$\mathrm{NLSQ\,T\,Tau\,Spectral\,Energy\,Distribution}$')
plt.xlim(0.05,1000.0)
plt.ylim(0.0003,15.0)
plt.legend(loc=2)

plt.savefig('SED_NLSQ.png', format='png')
plt.clf()

# Now Bayesian mode:

def power_law_model(x,f,err):
	a = mc.Uniform('a', lower=0.0, upper=1e3)
	b = mc.Uniform('b', lower=-100.0, upper=100.0)

	# same function for bayesian analysis
	@mc.deterministic(plot=False)
	def power_law(nu = x, a = a, b = b):
		return( a*np.power(x,b) )

	y = mc.Normal('power_law', mu=power_law, tau=1.0/np.power(err,2.0), value=f, observed=True)
	
	return( locals() )
	
def joint_model(x,f,err):
	em = mc.Uniform('EM', lower=1.0e4, upper=1.0e7)
	alpha1 = mc.Uniform('Alpha', lower=-2.1, upper=0.0)	# this "alpha" is referred to as eta in the paper
	t = mc.Uniform('T', lower=3.0e3, upper=2.0e4)
	alpha2 = mc.Uniform('Alpha_dust', lower=0.0, upper=4.0)
	norm2 = mc.Uniform('K', lower=10.0, upper=40.0)

	# same function for bayesian analysis
	@mc.deterministic(plot=False)
	def joint_model_flux(nu = x, em = em,alpha1 = alpha1, t = t, alpha2 = alpha2, norm2 = norm2):
		tau = tau_func(nu, t, em,alpha1)
		flux_density1 = t * c_omega * (1.0 - np.exp(-tau)) * np.power(nu, 2.0) # get flux density in mJy
		flux_density2 = norm2 * np.power((nu/100.0),alpha2)	# normalised power law
		total_flux_density = flux_density1 + flux_density2
		return( total_flux_density )

	y = mc.Normal('joint_model', mu=joint_model_flux, tau=1.0/np.power(err,2.0), value=f, observed=True)
	
	return( locals() )


#map = mc.MAP(model(df['Freq'].values, df['Flux'].values, df['Error'].values))
#map.fit()

#sampler = mc.MCMC( power_law_model(df['Freq'].values, df['Flux'].values, df['Error'].values) )
sampler = mc.MCMC( joint_model(df['Freq'].values, df['Flux'].values, df['Error'].values) )

sampler.sample(iter = 800000, burn= 200000)


# power law plotting
'''
sampler.a.summary()
sampler.b.summary()

plt.figure(2)
model_fluxes = sampler.stats()['a']['mean']*np.power(model_freqs,sampler.stats()['b']['mean'])
plt.loglog(model_freqs, model_fluxes)
plt.errorbar(df['Freq'].values, df['Flux'].values, yerr = df['Error'].values, linestyle='None', fmt='ko', ecolor='black')
'''

sampler.em.summary()
sampler.alpha1.summary()
sampler.t.summary()
sampler.alpha2.summary()
sampler.norm2.summary()
sampler.write_csv("fit_statistics.csv", variables=['EM', 'Alpha', 'T', 'Alpha_dust', 'K'])

# calculate parameters from EM

em, alpha1, t, alpha2, norm2 = sampler.stats()['EM']['mean'], sampler.stats()['Alpha']['mean'], sampler.stats()['T']['mean'], sampler.stats()['Alpha_dust']['mean'], sampler.stats()['K']['mean']
em_err, alpha1_err, t_err, alpha2_err, norm2_err = sampler.stats()['EM']['standard deviation'], sampler.stats()['Alpha']['standard deviation'], sampler.stats()['T']['standard deviation'], sampler.stats()['Alpha_dust']['standard deviation'], sampler.stats()['K']['standard deviation']

tau = 1.0
v0 = np.power(8.235e-2, -1.0/alpha1)*np.power(t,1.35/alpha1)*np.power(em, -1.0/alpha1)
ne = np.sqrt(em / (d*theta*(7.1e-3)))	# derive ne from em and theta
mion = (1.47500387e-9) * ne * np.power(theta*d,3.0)

# find error in these parameters

u = np.power(v0, alpha1)
omega_err = omega * np.sqrt(np.power(theta_p1_err/theta_p1,2.0) + np.power(theta_p2_err/theta_p2,2.0))
theta_err = np.power(2.0*theta_p1*theta_p2,-0.5)*np.sqrt(np.power(theta_p2*theta_p1_err,2.0) + np.power(theta_p1*theta_p2_err,2.0))
tau_err = tau * np.sqrt( np.power(em_err/em,2.0) + np.power( (alpha1*alpha1_err)/v0 ,2.0) + np.power( (1.35*t_err)/t,2.0) )
ne_err = ne * np.sqrt( np.power( (-0.5*d_err)/d,2.0) + np.power( (-0.5*theta_err)/theta,2.0) + np.power( (0.5*em_err)/em,2.0) )
mion_err = mion * np.sqrt( np.power(ne_err/ne, 2.0) + np.power((3.0*theta_err)/theta, 2.0) + np.power((3.0*d_err)/d, 2.0) )
v0_err = (v0/alpha1) * np.sqrt(np.power( (1.35*t_err)/t,2.0) + np.power(em_err/em,2.0) + np.power((np.log(-1.0/alpha1)*alpha1_err)/alpha1,2.0))

print('\n\n')
print('Bayesian parameters:')
print '\t Emission measure = '+str(em)+' +/- '+str(em_err)+' pc cm^-6'
print '\t Alpha (free-free) = '+str(alpha1)+' +/- '+str(alpha1_err)
print '\t Temperature = '+str(t)+' +/- '+str(t_err) + ' K'
print '\t Alpha (power law) = '+str(alpha2)+' +/- '+str(alpha2_err)
print '\t K (power law) = '+str(norm2)+' +/- '+str(norm2_err) + ' mJy'
print('\n')
print('Derived parameters:')
print('\t v0(GHz)   =   '+str(v0)+' +/- '+str(v0_err))
print('\t mion(solar masses)   =   '+str(mion)+' +/- '+str(mion_err))
print('\t n_e(cm-3)   =   '+str(ne)+' +/- '+str(ne_err))
print('\n\n')


# make new plots

model_fluxes1, model_fluxes2 = sed_func(model_freqs, em, alpha1, t, alpha2, norm2)
model_taus = tau_func(model_freqs, t, em, alpha1)

plt.loglog(model_freqs, model_taus)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Tau')
plt.title('T Tau opacity vs. frequency')

plt.savefig('Opacity_vs_Frequency.png', format='png')
plt.clf()

plt.loglog(model_freqs, model_fluxes1/1000.0,'k--',label=r'$\mathrm{Free-Free}$')
plt.loglog(model_freqs, model_fluxes2/1000.0,'k:',label=r'$\mathrm{Power\,law}$')
plt.loglog(model_freqs, model_fluxes1/1000.0 + model_fluxes2/1000.0,'k-',label=r'$\mathrm{Combined}$')
plt.errorbar(df['Freq'].values, df['Flux'].values/1000.0, yerr = df['Error'].values/1000.0, linestyle='None', fmt='ko', ecolor='black')
plt.xlabel(r'$\nu\,\mathrm{(GHz)}$', fontsize=myfontsize)
plt.ylabel(r'$S_{\nu}\,\mathrm{(Jy)}$', fontsize=myfontsize)
plt.title(r'$\mathrm{T\,Tau\,Spectral\,Energy\,Distribution}$', fontsize=myfontsize)
plt.xlim(0.05,1000.0)
plt.ylim(0.0003,15.0)
plt.legend(loc=2)

plt.savefig('SED_Bayes.png', format='png')
plt.clf()

samples = np.array([sampler.em.trace(),sampler.alpha1.trace(),sampler.t.trace(),sampler.alpha2.trace(),sampler.norm2.trace()]).T
corner.corner(samples,labels=[r'$\mathrm{EM\,(pc\,cm^{-6})}$',r'$\eta$',r'$\mathrm{T\,(K)}$',r'$\alpha_{dust}$',r'$\mathrm{K_{100\,GHz}\,(mJy)}$'], truths=[em, alpha1, t, alpha2, norm2])

plt.savefig('corner_plot.png', format='png')
plt.clf()

#mc.Matplot.plot(sampler)
#plt.show()
