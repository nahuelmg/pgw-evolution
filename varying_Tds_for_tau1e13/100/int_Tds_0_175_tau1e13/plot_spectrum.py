from matplotlib import pyplot as plt
import numpy as np
from cosmology import H0, mpl, a, a_eval, H_eval, H_inf, Neff_today
from params import h1_initial, h3_initial, spectrum_not_today, deg_s, rdof, reg_r0

#%%

kc_vec = np.loadtxt('kc_vec', unpack = True)

number_of_ranges = 20

h_prime_list = [np.loadtxt('h_prime_u0_{}'.format(i),unpack = True) for i in range(1,number_of_ranges + 1)]

h_prime = [sum(x) for x in zip(*h_prime_list)]

h_prime = np.array(h_prime,dtype=float)

if h1_initial:
    
    prim_spec = H_inf**2 / ( np.pi**2 * mpl**2 )
    
    if spectrum_not_today:
        
        prop = 1/(12*(a_eval*H_eval)**2)
        
    else:
        
        prop = 1/(12*(a[-1]*H0)**2)
    
    spectrum = prop * 2 * prim_spec * kc_vec**2. * h_prime**2.
        
    np.savetxt('spectrum_100_int_Tds_0_175_tau1e13',spectrum)

if rdof:
    
    from cosmology import g
    
    geff_inf = g[0]  

else:
    
    geff_inf = 100 # esto nunca debería activarse, lo pongo por las dudas para que no salte un error

sigma = 2 * deg_s * np.pi**2 / 30

mpl_hz = mpl * 1.519*1e24 # masa de Planck reducida en Hz

H_inf_hz = H_inf * 1.519*1e24 # H_inf en Hz



    
if h3_initial:
    
    if reg_r0:
        
        prim_spec = ( 1 / ( 23040 * np.pi ) ) * ( geff_inf**2 / sigma**2 ) * ( a[0] * H_inf_hz / kc_vec ) * ( kc_vec / ( a[0] * mpl_hz ) )**4 * np.heaviside( a[0] * H_inf_hz - kc_vec , 0 )
        
        if spectrum_not_today:
            
            prop = 1/(12*(a_eval*H_eval)**2)
            
        else:
            
            prop = 1/(12*(a[-1]*H0)**2)
        
        spectrum = prop * 2 * prim_spec * kc_vec**2. * h_prime**2.
            
        np.savetxt('spectrum_001_int_Tds_0_175_tau1e13',spectrum)
        
    else:
        
        prim_spec =  ( ( 1 / 3 * geff_inf**2 ) / ( sigma**2 * 18432 ) ) * ( kc_vec / ( a[0] * mpl_hz ) )**4 * np.heaviside( a[0] * H_inf_hz - kc_vec , 0 )
        
        if spectrum_not_today:
            
            prop = 1/(12*(a_eval*H_eval)**2)
            
        else:
            
            prop = 1/(12*(a[-1]*H0)**2)
        
        spectrum = prop * 2 * prim_spec * kc_vec**2. * h_prime**2.
            
        np.savetxt('spectrum_001_int_Tds_0_175_tau1e13',spectrum)
        
        
    

np.savetxt('kc_vec_int_Tds_0_175_tau1e13',kc_vec)

np.savetxt('Neff_today_int_Tds_0_175_tau1e13',[Neff_today, Neff_today])

#%%
#
##plt.figure()
###plt.figure(dpi=300)
#plt.loglog(kc_vec/(2.*np.pi),spectrum)
###plt.legend()
###plt.ylim(1e-41,10**(-35.5))
##plt.grid(True)
##plt.xlabel('f  [ Hz ]')
##plt.ylabel('$\Omega_{GW}(f,\eta_0)$')
###plt.title('pGW spectrum today $Neff={}$'.format(Neff_today))
###plt.savefig('pGW_spectrum')
