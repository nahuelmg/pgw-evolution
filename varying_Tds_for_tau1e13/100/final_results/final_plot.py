from matplotlib import pyplot as plt
import numpy as np

plt.style.use('default')

#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Helvetica']

#plt.rcParams['font.family'] = ['serif']
#plt.rcParams['font.serif'] = ['Times New Roman']
#plt.rcParams['text.usetex'] = False
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 12


#%%

# Carga de todos los datos

temperature_Tds_0_1 = np.loadtxt('temperature_Tds_0_1', unpack = True)
temperature_Tds_0_25 = np.loadtxt('temperature_Tds_0_25', unpack = True)
temperature_Tds_0_07 = np.loadtxt('temperature_Tds_0_07', unpack = True)
temperature_Tds_0_002 = np.loadtxt('temperature_Tds_0_002', unpack = True)
temperature_Tds_0_125 = np.loadtxt('temperature_Tds_0_125', unpack = True)
temperature_Tds_0_15 = np.loadtxt('temperature_Tds_0_15', unpack = True)
temperature_Tds_0_225 = np.loadtxt('temperature_Tds_0_225', unpack = True)
temperature_Tds_0_2 = np.loadtxt('temperature_Tds_0_2', unpack = True)
temperature_Tds_0_175 = np.loadtxt('temperature_Tds_0_175', unpack = True)

kc_horizon_crossing_Tds_0_1 = np.loadtxt('kc_horizon_crossing_Tds_0_1', unpack = True)
kc_horizon_crossing_Tds_0_25 = np.loadtxt('kc_horizon_crossing_Tds_0_25', unpack = True)
kc_horizon_crossing_Tds_0_07 = np.loadtxt('kc_horizon_crossing_Tds_0_07', unpack = True)
kc_horizon_crossing_Tds_0_002 = np.loadtxt('kc_horizon_crossing_Tds_0_002', unpack = True)
kc_horizon_crossing_Tds_0_125 = np.loadtxt('kc_horizon_crossing_Tds_0_125', unpack = True)
kc_horizon_crossing_Tds_0_15 = np.loadtxt('kc_horizon_crossing_Tds_0_15', unpack = True)
kc_horizon_crossing_Tds_0_225 = np.loadtxt('kc_horizon_crossing_Tds_0_225', unpack = True)
kc_horizon_crossing_Tds_0_2 = np.loadtxt('kc_horizon_crossing_Tds_0_2', unpack = True)
kc_horizon_crossing_Tds_0_175 = np.loadtxt('kc_horizon_crossing_Tds_0_175', unpack = True)

geff_Tds_0_1 = np.loadtxt('geff_Tds_0_1', unpack = True)
geff_Tds_0_25 = np.loadtxt('geff_Tds_0_25', unpack = True)
geff_Tds_0_07 = np.loadtxt('geff_Tds_0_07', unpack = True)
geff_Tds_0_002 = np.loadtxt('geff_Tds_0_002', unpack = True)
geff_Tds_0_125 = np.loadtxt('geff_Tds_0_125', unpack = True)
geff_Tds_0_15 = np.loadtxt('geff_Tds_0_15', unpack = True)
geff_Tds_0_225 = np.loadtxt('geff_Tds_0_225', unpack = True)
geff_Tds_0_2 = np.loadtxt('geff_Tds_0_2', unpack = True)
geff_Tds_0_175 = np.loadtxt('geff_Tds_0_175', unpack = True)

f_s_Tds_0_1 = np.loadtxt('f_s_Tds_0_1', unpack = True)
f_s_Tds_0_25 = np.loadtxt('f_s_Tds_0_25', unpack = True)
f_s_Tds_0_07 = np.loadtxt('f_s_Tds_0_07', unpack = True)
f_s_Tds_0_002 = np.loadtxt('f_s_Tds_0_002', unpack = True)
f_s_Tds_0_125 = np.loadtxt('f_s_Tds_0_125', unpack = True)
f_s_Tds_0_15 = np.loadtxt('f_s_Tds_0_15', unpack = True)
f_s_Tds_0_225 = np.loadtxt('f_s_Tds_0_225', unpack = True)
f_s_Tds_0_2 = np.loadtxt('f_s_Tds_0_2', unpack = True)
f_s_Tds_0_175 = np.loadtxt('f_s_Tds_0_175', unpack = True)


Neff_today_Tds_0_1 = np.loadtxt('Neff_today_int_Tds_0_1_tau1e13', unpack = True)[0]
Neff_today_Tds_0_25 = np.loadtxt('Neff_today_int_Tds_0_25_tau1e13', unpack = True)[0]
Neff_today_Tds_0_07 = np.loadtxt('Neff_today_int_Tds_0_07_tau1e13', unpack = True)[0]
Neff_today_Tds_0_002 = np.loadtxt('Neff_today_int_Tds_0_002_tau1e13', unpack = True)[0]
Neff_today_Tds_0_125 = np.loadtxt('Neff_today_int_Tds_0_125_tau1e13', unpack = True)[0]
Neff_today_Tds_0_15 = np.loadtxt('Neff_today_int_Tds_0_15_tau1e13', unpack = True)[0]
Neff_today_Tds_0_225 = np.loadtxt('Neff_today_int_Tds_0_225_tau1e13', unpack = True)[0]
Neff_today_Tds_0_2 = np.loadtxt('Neff_today_int_Tds_0_2_tau1e13', unpack = True)[0]
Neff_today_Tds_0_175 = np.loadtxt('Neff_today_int_Tds_0_175_tau1e13', unpack = True)[0]

kc_vec = np.loadtxt('kc_vec_int_Tds_0_1_tau1e13', unpack = True) # para todos el mismo

spectrum_100_int_Tds_0_1_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_1_tau1e13', unpack = True)
spectrum_100_int_Tds_0_25_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_25_tau1e13', unpack = True)
spectrum_100_int_Tds_0_07_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_07_tau1e13', unpack = True)
spectrum_100_int_Tds_0_002_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_002_tau1e13', unpack = True)
spectrum_100_int_Tds_0_125_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_125_tau1e13', unpack = True)
spectrum_100_int_Tds_0_15_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_15_tau1e13', unpack = True)
spectrum_100_int_Tds_0_225_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_225_tau1e13', unpack = True)
spectrum_100_int_Tds_0_2_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_2_tau1e13', unpack = True)
spectrum_100_int_Tds_0_175_tau1e13 = np.loadtxt('spectrum_100_int_Tds_0_175_tau1e13', unpack = True)
spectrum_100_no_int_Tds_0_1 = np.loadtxt('spectrum_100_no_int_Tds_0_1', unpack = True)
spectrum_100_no_int_Tds_0_25 = np.loadtxt('spectrum_100_no_int_Tds_0_25', unpack = True)
spectrum_100_no_int_Tds_0_07 = np.loadtxt('spectrum_100_no_int_Tds_0_07', unpack = True)
spectrum_100_no_int_Tds_0_002 = np.loadtxt('spectrum_100_no_int_Tds_0_002', unpack = True)
spectrum_100_no_int_Tds_0_125 = np.loadtxt('spectrum_100_no_int_Tds_0_125', unpack = True)
spectrum_100_no_int_Tds_0_15 = np.loadtxt('spectrum_100_no_int_Tds_0_15', unpack = True)
spectrum_100_no_int_Tds_0_225 = np.loadtxt('spectrum_100_no_int_Tds_0_225', unpack = True)
spectrum_100_no_int_Tds_0_2 = np.loadtxt('spectrum_100_no_int_Tds_0_2', unpack = True)
spectrum_100_no_int_Tds_0_175 = np.loadtxt('spectrum_100_no_int_Tds_0_175', unpack = True)


#%%
#para grafico sin rdofs

freq_100_no_rdof_Tds_0_125 = kc_vec/(2*np.pi)
mask = np.logical_and(10**(-15)<freq_100_no_rdof_Tds_0_125,freq_100_no_rdof_Tds_0_125<10**(-14.8))


spec_1 = spectrum_100_no_int_Tds_0_125[0:12661].tolist()
spec_2 = spectrum_100_no_int_Tds_0_125[mask].tolist()

spec = []
spec.extend(spec_1)

for i in range(0,int((50000-12000)/666) ):
    spec.extend(spec_2)

spec = np.asarray(spec[0:50000])

plt.loglog(freq_100_no_rdof_Tds_0_125,spec)

#%%

# Análisis para Tds = 0.15 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_15[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_15_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_15,linewidth=2,label='Free evolution')
plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_15_tau1e13,linewidth=2,label='Viscous effects')
plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_15_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_15_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_15_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_15_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_15/(2*np.pi),f_s_Tds_0_15,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()


plt.figure()
plt.semilogx(temperature_Tds_0_15,f_s_Tds_0_15,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)

plt.figure()
plt.semilogx(temperature_Tds_0_15,geff_Tds_0_15)
plt.title('Neff={0:.2f}'.format(Neff_today_Tds_0_15))
plt.xlabel(r'$T\;\;\;{\rm [ GeV ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()

plt.figure()
plt.semilogx(kc_horizon_crossing_Tds_0_15/(2*np.pi),geff_Tds_0_15)
plt.title('Neff={0:.2f}'.format(Neff_today_Tds_0_15))
plt.axvline(1/(1e13*np.pi*2))
plt.xlabel(r'$f\;\;\;{\rm [ Hz ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()


#%%

# Análisis para Tds = 0.175 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_175[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_175_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_175,linewidth=2,label='Free evolution')
plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_175_tau1e13,linewidth=2,label='Viscous effects')
plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_175_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_175_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_175_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_175_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_175/(2*np.pi),f_s_Tds_0_175,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()


plt.figure()
plt.semilogx(temperature_Tds_0_175,f_s_Tds_0_175,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)

plt.figure()
plt.semilogx(temperature_Tds_0_175,geff_Tds_0_175)
plt.title('Neff={0:.2f}'.format(Neff_today_Tds_0_175))
plt.xlabel(r'$T\;\;\;{\rm [ GeV ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()

plt.figure()
plt.semilogx(kc_horizon_crossing_Tds_0_175/(2*np.pi),geff_Tds_0_175)
plt.title('Neff={0:.2f}'.format(Neff_today_Tds_0_175))
plt.axvline(1/(1e13*np.pi*2))
plt.xlabel(r'$f\;\;\;{\rm [ Hz ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()



#%%

# Análisis para Tds = 0.125 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_125[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_125_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,label='Free evolution')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,label='Viscous effects')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_125_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_125_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_125_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_125_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_125/(2*np.pi),f_s_Tds_0_125,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()


plt.figure()
plt.semilogx(temperature_Tds_0_125,f_s_Tds_0_125,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)



plt.figure()
plt.semilogx(temperature_Tds_0_125,geff_Tds_0_125)
plt.title(r'$N_{\rm eff}={}$')#.format(Neff_today_Tds_0_125))
plt.xlabel(r'$T\;\;\;{\rm [ GeV ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()

plt.figure()
plt.semilogx(kc_horizon_crossing_Tds_0_125/(2*np.pi),geff_Tds_0_125)
plt.title(r'$N_{\rm eff}={}$')#.format(Neff_today_Tds_0_125))
plt.axvline(1/(1e13*np.pi*2))
plt.xlabel(r'$f\;\;\;{\rm [ Hz ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()


#%%

# Análisis para Tds = 0.07 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_07[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_07_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_07,linewidth=2,label='Free evolution')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_07_tau1e13,linewidth=2,label='Viscous effects')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_07_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_07_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_07_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_07_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_07/(2*np.pi),f_s_Tds_0_07,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()

#plt.figure()
#plt.semilogx(temperature_Tds_0_07,geff_Tds_0_07)
#plt.title(r'$N_{\rm eff}={}$'.format(Neff_today_Tds_0_07))
#plt.xlabel(r'$T\;\;\;{\rm[ GeV ]}$')
#plt.ylabel(r'$g_*$')
#plt.legend()
#plt.grid()
#
#plt.figure()
#plt.semilogx(kc_horizon_crossing_Tds_0_07/(2*np.pi),geff_Tds_0_07)
#plt.title(r'$N_{\rm eff}={}$'.format(Neff_today_Tds_0_07))
#plt.axvline(1/(1e13*np.pi*2))
#plt.xlabel(r'$f\;\;\;{\rm[ Hz ]}$')
#plt.ylabel(r'$g_*$')
#plt.legend()
#plt.grid()

#%%

# Análisis para Tds = 0.002 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_002[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_002_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_002,linewidth=2,label='Free evolution')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_002_tau1e13,linewidth=2,label='Viscous effects')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_002_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_002_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_002_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_002_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_002/(2*np.pi),f_s_Tds_0_002,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()

plt.figure()
plt.semilogx(temperature_Tds_0_002,geff_Tds_0_002)
plt.title(r'$N_{\rm eff}={}$')#.format(Neff_today_Tds_0_002))
plt.xlabel(r'$T\;\;\;{\rm[ GeV ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()

plt.figure()
plt.semilogx(kc_horizon_crossing_Tds_0_002/(2*np.pi),geff_Tds_0_002)
plt.title(r'$N_{\rm eff}={}$')#.format(Neff_today_Tds_0_002))
plt.axvline(1/(1e13*np.pi*2))
plt.xlabel(r'$f\;\;\;{\rm[ Hz ]}$')
plt.ylabel(r'$g_*$')
plt.legend()
plt.grid()

#%%

# Análisis para Tds = 0.1 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_1[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_1_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_1,linewidth=2,label='Free evolution')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_1_tau1e13,linewidth=2,label='Viscous effects')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_1_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_1_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_1_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_1_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_1/(2*np.pi),f_s_Tds_0_1,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()

#plt.figure()
#plt.semilogx(temperature_Tds_0_1,geff_Tds_0_1)
#plt.title(r'$N_{\rm eff}={}$'.format(Neff_today_Tds_0_1))
#plt.xlabel(r'$T\;\;\;{\rm[ GeV ]}$')
#plt.ylabel(r'$g_*$')
#plt.legend()
#plt.grid()
#
#plt.figure()
#plt.semilogx(kc_horizon_crossing_Tds_0_1/(2*np.pi),geff_Tds_0_1)
#plt.title(r'$N_{\rm eff}={}$'.format(Neff_today_Tds_0_1))
#plt.axvline(1/(1e13*np.pi*2))
#plt.xlabel(r'$f\;\;\;{\rm[ Hz ]}$')
#plt.ylabel(r'$g_*$')
#plt.legend()
#plt.grid()


#%%

# Análisis para Tds = 0.25 GeV

n = 300

delta_i = int(len(kc_vec)/n)

range_i = [j*delta_i for j in range(n)]

maxs_no_int = []
maxs_position_no_int = []

for i in range_i:
    
    step = spectrum_100_no_int_Tds_0_25[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_no_int.extend(maxs_per_step)
    maxs_position_no_int.extend(maxs_position_per_step) 


maxs_int = []
maxs_position_int = []

for i in range_i:
    
    step = spectrum_100_int_Tds_0_25_tau1e13[i:i+int(len(kc_vec)/n)]
    maxs_per_step = [max(step)]
    maxs_position_per_step = [k+i for k, j in enumerate(step) if j == maxs_per_step]
    
    maxs_int.extend(maxs_per_step)
    maxs_position_int.extend(maxs_position_per_step)
    

plt.figure()
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_25,linewidth=2,label='Free evolution')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_25_tau1e13,linewidth=2,label='Viscous effects')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.legend()
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')

deltaomega_omega_Tds_0_25_tau1e13 = np.abs(np.asarray(maxs_int) - np.asarray(maxs_no_int))/((np.asarray(maxs_int) + np.asarray(maxs_no_int))/2)
kc_vec_deltaomega_omega_Tds_0_25_tau1e13 = kc_vec[maxs_position_int]

import matplotlib.ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_25_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_25_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_25/(2*np.pi),f_s_Tds_0_25,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()

#plt.figure()
#plt.semilogx(temperature_Tds_0_25,geff_Tds_0_25)
#plt.title(r'$N_{\rm eff}={}$'.format(Neff_today_Tds_0_25))
#plt.xlabel(r'$T\;\;\;{\rm[ GeV ]}$')
#plt.ylabel(r'$g_*$')
#plt.legend()
#plt.grid()
#
#plt.figure()
#plt.semilogx(kc_horizon_crossing_Tds_0_25/(2*np.pi),geff_Tds_0_25)
#plt.title(r'$N_{\rm eff}={}$'.format(Neff_today_Tds_0_25))
#plt.axvline(1/(1e13*np.pi*2))
#plt.xlabel(r'$f\;\;\;{\rm[ Hz ]}$')
#plt.ylabel(r'$g_*$')
#plt.legend()
#plt.grid()

#%%

# GRAFICO COMPARATIVO

from matplotlib import ticker as mtick

fig = plt.figure()
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_25_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_25_tau1e13,label=r'$\Delta \Omega/\Omega\;(T_{ds}=0.25\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_25/(2*np.pi),f_s_Tds_0_25,label=r'$\rho_s/\rho_{\rm c}\;(T_{ds}=0.25\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_07_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_07_tau1e13,label=r'$\Delta \Omega/\Omega\;(T_{ds}=0.07\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_07/(2*np.pi),f_s_Tds_0_07,label=r'$\rho_s/\rho_{\rm c}\;(T_{ds}=0.07\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_1_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_1_tau1e13,label=r'$\Delta \Omega/\Omega\;(T_{ds}=0.1\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_1/(2*np.pi),f_s_Tds_0_1,label=r'$\rho_s/\rho_{\rm c}\;(T_{ds}=0.1\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_002_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_002_tau1e13,label=r'$\Delta \Omega/\Omega\;(T_{ds}=0.002\,{\rm GeV})$',linewidth=2)
plt.semilogx(kc_horizon_crossing_Tds_0_002/(2*np.pi),f_s_Tds_0_002,label=r'$\rho_s/\rho_{\rm c}\;(T_{ds}=0.002\,{\rm GeV})$',linewidth=2)
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend(loc=2)
plt.grid()


#%%

#PARA IMPORTAR


from mpl_toolkits.axes_grid.inset_locator import inset_axes

fig = plt.figure(dpi=300#,facecolor='white')#figsize=(9, 4),facecolor='white')
                    )
ax = fig.add_subplot(111)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,label=r'${\rm Free}$ ${\rm evolution}$')#,color='black')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,label=r'${\rm Viscous}$ ${\rm effects}$')#,color='dimgray')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.axvline(1/(1e13*np.pi*2),ymax=0.5,color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.legend(loc=2,facecolor='white')
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')
#plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])

# this is an inset axes over the main axes
inset_axes = inset_axes(ax, 
                    width="50%", # width = 30% of parent_bbox
                    height="40%", # height : 1 inch
                    loc=1)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,color='C0')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,color='C1')
plt.xlim(1*10**(-17),3*10**(-5))
plt.ylim(2.9*10**(-16),1.1*10**(-15))
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('pGW_spectra_100_Tds_0_125_tau1e13')


import matplotlib.ticker as mtick

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_125_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_125_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)#,color='black')
plt.semilogx(kc_horizon_crossing_Tds_0_125/(2*np.pi),f_s_Tds_0_125,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)#,color='dimgray')
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('delta_omega_Tds_0_125_tau1e13')


#%%

#para importar en español

from mpl_toolkits.axes_grid.inset_locator import inset_axes

fig = plt.figure(dpi=300#,facecolor='white')#figsize=(9, 4),facecolor='white')
                    )
ax = fig.add_subplot(111)
plt.loglog(kc_vec/(2.*np.pi),spec,linewidth=2,label=r'${\rm Libre}$ ${\rm s/GdLR}$',color='silver')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,label=r'${\rm Libre}$')#,color='black')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,label=r'${\rm Efectos}$ ${\rm viscosos}$')#,color='dimgray')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.axvline(1/(1e13*np.pi*2),ymax=0.5,color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.annotate(r'${\rm Aniquilación}\; e^{\pm}$',
         xy=(10**(-11.8),  8*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(10**(-11.8), 3.7*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Transición}\; {\rm QCD}$',
         xy=(4*10**(-9),  4.7*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(4*10**(-9), 3.5*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Transición}\; {\rm EW}$',
         xy=(1.7*10**(-6),  3.4*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-6), 2*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.legend(loc=2,facecolor='white')
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')
#plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])

# this is an inset axes over the main axes
inset_axes = inset_axes(ax, 
                    width="50%", # width = 30% of parent_bbox
                    height="40%", # height : 1 inch
                    loc=1)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,color='C0')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,color='C1')
plt.xlim(1*10**(-17),3*10**(-5))
plt.ylim(2.9*10**(-16),1.1*10**(-15))
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('pGW_spectra_100_Tds_0_125_tau1e13_tesis')

#%%
import matplotlib.ticker as mtick

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_125_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_125_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)#,color='black')
plt.semilogx(kc_horizon_crossing_Tds_0_125/(2*np.pi),f_s_Tds_0_125,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)#,color='dimgray')
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.annotate(r'${\rm Aniq.}\; e^{\pm}$',
         xy=(2.7*10**(-12),  0.07),
         ha='center',
         #xycoords='data',
         xytext=(2.7*10**(-12), 0.09),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Desacople}\;s$',
         xy=(1.7*10**(-8),  0.021),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-8), 0.005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}$',
         xy=(4*10**(-9),  0.04),
         ha='center',
         #xycoords='data',
         xytext=(4*10**(-9), 0.021),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm EW}$',
         xy=(1.3*10**(-6),  0.018),
         ha='center',
         #xycoords='data',
         xytext=(1.3*10**(-6), 0.0005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Igualdad}\; {\rm Mat-Rad}$',
         xy=(9*10**(-17),  0.06),
         ha='center',
         #xycoords='data',
         xytext=(9*10**(-17), 0.081),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('delta_omega_Tds_0_125_tau1e13_tesis')

#%%

#para importar en inglés

from mpl_toolkits.axes_grid.inset_locator import inset_axes

fig = plt.figure(dpi=300#,facecolor='white')#figsize=(9, 4),facecolor='white')
                    )
ax = fig.add_subplot(111)
#plt.loglog(kc_vec/(2.*np.pi),spec,linewidth=2,label=r'${\rm Libre}$ ${\rm s/GdLR}$',color='silver')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,label=r'${\rm Free}$ ${\rm evolution}$')#,color='black')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,label=r'${\rm Viscous}$ ${\rm effects}$')#,color='dimgray')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.axvline(1/(1e13*np.pi*2),ymax=0.5,color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.annotate(r'${\rm e^{\pm}\; Annih.}\; $',
         xy=(10**(-11.8),  8*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(10**(-11.8), 3.7*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}\; {\rm Transition} $',
         xy=(4*10**(-9),  4.7*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(4*10**(-9), 3.5*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'$ {\rm EW}\;{\rm Transition}\;$',
         xy=(1.7*10**(-6),  3.4*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-6), 2*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.legend(loc=2,facecolor='white')
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')
#plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])

# this is an inset axes over the main axes
inset_axes = inset_axes(ax, 
                    width="50%", # width = 30% of parent_bbox
                    height="40%", # height : 1 inch
                    loc=1)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_125,linewidth=2,color='C0')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_125_tau1e13,linewidth=2,color='C1')
plt.xlim(1*10**(-17),3*10**(-5))
plt.ylim(2.9*10**(-16),1.1*10**(-15))
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('pGW_spectra_100_Tds_0_125_tau1e13_ingles')

#%%
import matplotlib.ticker as mtick

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_125_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_125_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)#,color='black')
plt.semilogx(kc_horizon_crossing_Tds_0_125/(2*np.pi),f_s_Tds_0_125,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)#,color='dimgray')
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.annotate(r'$e^{\pm}\;{\rm Annih.}\; $',
         xy=(2.7*10**(-12),  0.07),
         ha='center',
         #xycoords='data',
         xytext=(2.7*10**(-12), 0.09),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'$s\;{\rm Decoupling}$',
         xy=(1.7*10**(-8),  0.021),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-8), 0.005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}$',
         xy=(4*10**(-9),  0.04),
         ha='center',
         #xycoords='data',
         xytext=(4*10**(-9), 0.021),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm EW}$',
         xy=(1.3*10**(-6),  0.018),
         ha='center',
         #xycoords='data',
         xytext=(1.3*10**(-6), 0.0005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Equality}$',
         xy=(9*10**(-17),  0.06),
         ha='center',
         #xycoords='data',
         xytext=(9*10**(-17), 0.081),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('delta_omega_Tds_0_125_tau1e13_ingles')

#%%

#para importar en inglés 0.175

from mpl_toolkits.axes_grid.inset_locator import inset_axes

fig = plt.figure(dpi=300#,facecolor='white')#figsize=(9, 4),facecolor='white')
                    )
ax = fig.add_subplot(111)
#plt.loglog(kc_vec/(2.*np.pi),spec,linewidth=2,label=r'${\rm Libre}$ ${\rm s/GdLR}$',color='silver')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_175,linewidth=2,label=r'${\rm Free}$ ${\rm evolution}$')#,color='black')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_175_tau1e13,linewidth=2,label=r'${\rm Viscous}$ ${\rm effects}$')#,color='dimgray')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.axvline(1/(1e13*np.pi*2),ymax=0.5,color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.annotate(r'${\rm e^{\pm}\; Annih.}\; $',
         xy=(10**(-11.8),  8*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(10**(-11.8), 3.7*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}\; {\rm Transition} $',
         xy=(4*10**(-9),  4.7*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(4*10**(-9), 3.5*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'$ {\rm EW}\;{\rm Transition}\;$',
         xy=(1.7*10**(-6),  3.4*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-6), 2*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.legend(loc=2,facecolor='white')
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')
#plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])

# this is an inset axes over the main axes
inset_axes = inset_axes(ax, 
                    width="50%", # width = 30% of parent_bbox
                    height="40%", # height : 1 inch
                    loc=1)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_175,linewidth=2,color='C0')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_175_tau1e13,linewidth=2,color='C1')
plt.xlim(1*10**(-17),3*10**(-5))
plt.ylim(2.9*10**(-16),1.1*10**(-15))
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('pGW_spectra_100_Tds_0_175_tau1e13_ingles')

#%%
import matplotlib.ticker as mtick

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_175_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_175_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)#,color='black')
plt.semilogx(kc_horizon_crossing_Tds_0_175/(2*np.pi),f_s_Tds_0_175,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)#,color='dimgray')
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.annotate(r'$e^{\pm}\;{\rm Annih.}\; $',
         xy=(2.7*10**(-12),  0.041),
         ha='center',
         #xycoords='data',
         xytext=(2.7*10**(-12), 0.061),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'$s\;{\rm Decoupling}$',
         xy=(1.8*10**(-8),  0.025),
         ha='center',
         #xycoords='data',
         xytext=(1.8*10**(-8), 0.007),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}$',
         xy=(2.5*10**(-9),  0.04),
         ha='center',
         #xycoords='data',
         xytext=(2.5*10**(-9), 0.021),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm EW}$',
         xy=(1.3*10**(-6),  0.018),
         ha='center',
         #xycoords='data',
         xytext=(1.3*10**(-6), 0.0005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Equality}$',
         xy=(9*10**(-17),  0.035),
         ha='center',
         #xycoords='data',
         xytext=(9*10**(-17), 0.051),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('delta_omega_Tds_0_175_tau1e13_ingles')

#%%

#para importar en inglés 0.15

from mpl_toolkits.axes_grid.inset_locator import inset_axes

fig = plt.figure(dpi=300#,facecolor='white')#figsize=(9, 4),facecolor='white')
                    )
ax = fig.add_subplot(111)
#plt.loglog(kc_vec/(2.*np.pi),spec,linewidth=2,label=r'${\rm Libre}$ ${\rm s/GdLR}$',color='silver')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_15,linewidth=2,label=r'${\rm Free}$ ${\rm evolution}$')#,color='black')
#plt.loglog(kc_vec[maxs_position_no_int]/(2*np.pi),maxs_no_int)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_15_tau1e13,linewidth=2,label=r'${\rm Viscous}$ ${\rm effects}$')#,color='dimgray')
#plt.loglog(kc_vec[maxs_position_int]/(2*np.pi),maxs_int)
plt.axvline(1/(1e13*np.pi*2),ymax=0.5,color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.annotate(r'${\rm e^{\pm}\; Annih.}\; $',
         xy=(10**(-11.8),  8*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(10**(-11.8), 3.7*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}\; {\rm Transition} $',
         xy=(4*10**(-9),  4.7*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(4*10**(-9), 3.5*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'$ {\rm EW}\;{\rm Transition}\;$',
         xy=(1.7*10**(-6),  3.4*10**(-16)),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-6), 2*10**(-15)),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.legend(loc=2,facecolor='white')
plt.ylim(10**(-16),10**(-12))
plt.grid(True, which='major')
plt.xlabel(r'$f\;{\rm [ Hz ]}$')
plt.ylabel(r'$\Omega_{\rm GW}(f,\eta_0)$')
#plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])

# this is an inset axes over the main axes
inset_axes = inset_axes(ax, 
                    width="50%", # width = 30% of parent_bbox
                    height="40%", # height : 1 inch
                    loc=1)
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_no_int_Tds_0_15,linewidth=2,color='C0')
plt.loglog(kc_vec/(2.*np.pi),spectrum_100_int_Tds_0_15_tau1e13,linewidth=2,color='C1')
plt.xlim(1*10**(-17),3*10**(-5))
plt.ylim(2.9*10**(-16),1.1*10**(-15))
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_\tau$')
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('pGW_spectra_100_Tds_0_15_tau1e13_ingles')

#%%
import matplotlib.ticker as mtick

fig = plt.figure(dpi=300)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
plt.semilogx(kc_vec_deltaomega_omega_Tds_0_15_tau1e13/(2*np.pi),deltaomega_omega_Tds_0_15_tau1e13,label=r'$\Delta \Omega/\Omega$',linewidth=2)#,color='black')
plt.semilogx(kc_horizon_crossing_Tds_0_15/(2*np.pi),f_s_Tds_0_15,label=r'$\rho_s/\rho_{\rm c}$',linewidth=2)#,color='dimgray')
plt.axvline(1/(1e13*np.pi*2),color='black',linestyle='--',linewidth=1,label=r'$f_{\tau}$')
plt.annotate(r'$e^{\pm}\;{\rm Annih.}\; $',
         xy=(2.7*10**(-12),  0.055),
         ha='center',
         #xycoords='data',
         xytext=(2.7*10**(-12), 0.07),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'$s\;{\rm Decoupling}$',
         xy=(1.7*10**(-8),  0.021),
         ha='center',
         #xycoords='data',
         xytext=(1.7*10**(-8), 0.005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm QCD}$',
         xy=(3*10**(-9),  0.04),
         ha='center',
         #xycoords='data',
         xytext=(3*10**(-9), 0.021),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm EW}$',
         xy=(1.3*10**(-6),  0.018),
         ha='center',
         #xycoords='data',
         xytext=(1.3*10**(-6), 0.0005),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.annotate(r'${\rm Equality}$',
         xy=(9*10**(-17),  0.047),
         ha='center',
         #xycoords='data',
         xytext=(9*10**(-17), 0.061),
         fontsize=10,
         arrowprops=dict(arrowstyle="->"))
plt.xlim(10**(-20),10**(-2))
plt.xlabel(r'$f\;{\rm[ Hz ]}$')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('delta_omega_Tds_0_15_tau1e13_ingles')