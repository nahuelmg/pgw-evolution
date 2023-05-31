#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 18:42:59 2021

@author: nahue
"""
from matplotlib import pyplot as plt
import numpy as np


kc_vec = np.loadtxt('kc_vec_damping_1')
spectrum_damping_1 = np.loadtxt('spectrum_100_damping_1')
spectrum_damping_2 = np.loadtxt('spectrum_100_damping_2')
spectrum_damping_3 = np.loadtxt('spectrum_100_damping_3')
spectrum_no_damping = np.loadtxt('spectrum_100_no_damping')
spectrum_damping_4 = np.loadtxt('spectrum_100_damping_4')
spectrum_damping_5 = np.loadtxt('spectrum_100_damping_5')


plt.figure()
plt.figure(dpi=300)
plt.loglog(kc_vec/(2.*np.pi),spectrum_no_damping,label='No damping')
plt.loglog(kc_vec/(2.*np.pi),spectrum_damping_3,label='Damping DTT ($K=cont$)')
plt.loglog(kc_vec/(2.*np.pi),spectrum_damping_4,label='Damping DTT ($K=cont$) 4')
plt.loglog(kc_vec/(2.*np.pi),spectrum_damping_1,label='Damping $K\propto \delta$')
plt.loglog(kc_vec/(2.*np.pi),spectrum_damping_5,label='Damping DTT ($K=cont$) 5')
plt.legend()
##plt.ylim(1e-41,10**(-35.5))
plt.grid(True)
plt.xlabel('f  [ Hz ]')
plt.ylabel('$\Omega_{GW}(f,\eta_0)$')
#plt.title('pGW spectrum today $Neff={}$'.format(Neff_today))
plt.savefig('pGW_spectrum')