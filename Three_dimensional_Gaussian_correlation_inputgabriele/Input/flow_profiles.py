# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 10:45:28 2018

@author: Gabriele
"""

import scipy as sp
import matplotlib.pyplot as plt

data = sp.io.loadmat('PL_input.mat')
ypos = data['ypos']
velmean = float(data['velmean'])
velgrad = data['velgrad']
delta99 = float(data['delta99'])
ufluc = data['ufluc']
Lam_input = data['Lambda']

fig, ax = plt.subplots()
ax.plot(ypos/delta99, velgrad, 'kx')
ax.set_xlim([0., 1.])
ax.set_xlabel(r'$x_2/\delta$', fontsize=18)
ax.set_ylabel(r'$\partial U_1  / \partial x_2$', fontsize=18)
ax.grid(True)
fig.savefig("CFD_U.eps")
fig.savefig("CFD_U.png") 

fig, ax = plt.subplots()
ax.plot(ypos/delta99, ufluc, 'kx')
ax.set_xlim([0., 1.])
ax.set_xlabel(r'$x_2/\delta$', fontsize=18)
ax.set_ylabel(r'$\hat{u}_2$', fontsize=18)
ax.grid(True)
fig.savefig("CFD_u2.eps")
fig.savefig("CFD_u2.png") 

fig, ax = plt.subplots()
ax.plot(ypos/delta99, 1.4*Lam_input/delta99, 'kx')
ax.set_xlim([0., 1.])
ax.set_ylim([0., 0.2])
ax.set_xlabel(r'$x_2/\delta$', fontsize=18)
ax.set_ylabel(r'$\Lambda/\delta$', fontsize=18)
ax.grid(True)
fig.savefig("CFD_Lambda.eps")
fig.savefig("CFD_Lambda.png") 