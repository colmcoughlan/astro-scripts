#!/usr/bin/python3
# -*- coding: utf-8 -*-
# python3 ttau_density.py
#
# exec(open('ttau_gas_parameters_new_eq.py').read())
#
# ... ttau_gas_parameters           v0.2                           20.Apr.2016
#
# parameter estimation of parameters of an ionised gas
# using the equations in Coughlan et al. 2016, ApJ
#
# electron density
#
# rewritten to 
# T_e = (5.84685 x 10^-4 x D x theta x n_e^2 x nu^eta)**(1/1.35)
# with
# n_e   / cm⁻3     electron density
# T_e   / K        electron temperature
# nu    / GHz      frequency observed
# this is the turnover frequency nu = 157 MHz
# with free-free spectral index eta = -1.83
# D     / kpc      distance to source
# theta / arcsec   size of source
# 
#
# followed by 
# ionised gas mass
#
# rewritten to 
# T_e = (3.355893054 x 10^15 x M_ion^2 / (D^5 x theta^5) x nu^eta)**(1/1.35)
# with
# M_ion / M_sun    ionised gas mass

from matplotlib import rc
from matplotlib.patches import Ellipse
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
myfontsize = 15
#rc('font',**{'family':'serif','serif':['Times']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator_x = MultipleLocator(1)
majorFormatter_x = FormatStrFormatter('%d')
minorLocator_x = MultipleLocator(0.5)
majorLocator_y = MultipleLocator(10000)
majorFormatter_y = FormatStrFormatter('%d')
minorLocator_y = MultipleLocator(5000)
fig, ax = plt.subplots()

# used values
# nu = 157 MHz
nu = 0.157
# eta = -1.83
eta = -1.83
# S_nu = 1.75 mJy
S_nu = 1.75
# D = 140 pc
D = 0.14
# n_e = x cm^-3
n_e = [5e3,7.5e3,1e4,1.5e4]
# M_ion = x M_sun
M_ion = [1e-7,1e-6,4e-6,1e-5]
# theta_min = 0.5 arcsec, theta_max = 5 arcsec
theta_min = 0.5
theta_max = 10.0
# limits in temperture to be plotted
T_e_min = 1000
T_e_max = 10000

theta_list = np.zeros(100)
temp_list1  = np.zeros([len(n_e),len(theta_list)])
temp_list2  = np.zeros([len(n_e),len(theta_list)])

# create a list for theta
theta_steps = 100
delta_theta = (theta_max - theta_min) / theta_steps

for i in range(theta_steps):
   theta_list[i] = theta_min + i * delta_theta 

#
# electron density
#
const = (5.84685e-4 * D * nu**eta)

for j in range(len(n_e)):
   for i in range(theta_steps):
      temp           = (const * n_e[j]**2 * theta_list[i])**(1/1.35)
      temp_list1[j,i] = temp

ax.xaxis.set_major_locator(majorLocator_x)
ax.xaxis.set_major_formatter(majorFormatter_x)
ax.yaxis.set_major_locator(majorLocator_y)
ax.yaxis.set_major_formatter(majorFormatter_y)

# for the minor ticks, use no labels; default NullFormatter
ax.xaxis.set_minor_locator(minorLocator_x)
ax.yaxis.set_minor_locator(minorLocator_y)

#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
plt.text(8.8, 15500, r'5$\times$10$^3$', fontsize=14)
plt.text(8.0, 28000, r'7.5$\times$10$^3$', fontsize=14)
plt.text(6.9, 39000, r'1.0$\times$10$^4$', fontsize=14)
plt.text(3.5, 43000, r'n$_e$ = 1.5$\times$10$^4$', fontsize=14)
#plt.xticks([0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6], \\
# [r$0$,r$1$,r$2$,r$3$,r$4$,r$5$])
plt.xticks([0,1,2,3,4,5,6,7,8,9,10], [r'$0$',r'$1$',r'$2$',r'$3$',r'$4$',r'$5$',r'$6$',r'$7$',r'$8$',r'$9$',r'$10$'], fontsize=15)
#plt.yticks([])
#plt.xlabel(r'{$\theta$ [arcsec]}', fontsize=18)
#plt.ylabel(r'{T$_e$ [K]}', fontsize=18)
plt.axis([0,10,0,50000])

plt.xlabel(r'$\mathrm{\theta}\,\mathrm{(arcsec)}$', fontsize=myfontsize)
plt.ylabel(r'$\mathrm{T_{e}}\,\mathrm{(K)}$', fontsize=myfontsize)

for i in range(len(n_e)):
#   line, = plt.plot(theta_list,temp_list1[i], 'k-')
   line, = plt.plot(theta_list,temp_list1[i], 'k--')

#plt.savefig('tex_demo')
###plt.show()
#fig = plt.figure()
#fig.savefig('ttau.png')
#
# ionised gas mass
#
const = (2.687420543e14 / D**5 * nu**eta)

for j in range(len(n_e)):
   for i in range(theta_steps):
      temp           = (const * M_ion[j]**2 / theta_list[i]**5)**(1/1.35)
      temp_list2[j,i] = temp

plt.text(8.8, 11000, r'1$\times$10$^{-5}$', fontsize=14)
plt.text(8.6, 3000, r'4$\times$10$^{-6}$', fontsize=14)
plt.text(6.0, 2000, r'1$\times$10$^{-6}$', fontsize=14)
plt.text(2.2, 2200, r'M$_{ion}$ = 10$^{-7}$', fontsize=14)

for i in range(len(n_e)):
   line, = plt.plot(theta_list,temp_list2[i], 'k-')
#   line, = plt.plot(theta_list,temp_list2[i], 'k--')


# Colm: Add ellipse indicating region

ellipse = Ellipse(xy = (3.26, 14339.1562498), width = 2.0* 2.59, height = 2.0 * 3245.33287323)
ellipse.set_alpha(0.5)
ellipse.set_facecolor('gray')
ax = plt.gca()
ax.add_patch(ellipse)


plt.savefig('ttau_gas_density_mass.png')
#plt.show()
#fig = plt.figure()
#fig.savefig('ttau.png')
