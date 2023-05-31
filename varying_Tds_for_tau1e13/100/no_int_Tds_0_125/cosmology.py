#%%
import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import quad as quad
from scipy.interpolate import UnivariateSpline as uspl
#plt.ion()
from params import h, T0, V_phi_1_4, T_ls, T_nu, ms, Tds, omega_k0, omega_lambda0, omega_gamma0, do_s_massive, deg_s, do_s_nomassive, do_cosmo_s, geff_today, rdof, kill_interaction_after_decoupling, T_eval


#%% Parámetros cosmológicos relevantes

if rdof:
    
    H0 = h * 3.241e-18 # en Hz
    
    T, g, gs = np.loadtxt('effective_relativistic_dof_saikawa+18.dat',usecols=(0,1,3),unpack=True)
    
    # T en GeV
    # OJO mirar bien la suavidad de g y gs
    
    # extiendo T, g y gs hasta hoy
    
    T_log_res = (np.log10(T[-1])-np.log10(T[0])) / len(T)
    
    number_T_ext = (np.log10(T[0])-T_log_res-np.log10(T0)) / T_log_res
    
    T_ext = np.logspace(np.log10(T0), np.log10(T[0])-T_log_res, number_T_ext)
    
    T = T[::-1]             # invierte el orden
    T_ext = T_ext[::-1]
    T = np.concatenate((T,T_ext), axis=None)
    #T = T[::-1]
    
    g = g[::-1]
    g = np.concatenate(( g, g[-1] * np.ones(len(T_ext))), axis=None)
    #g = g[::-1]
    
    gs = gs[::-1]
    gs = np.concatenate(( gs, gs[-1] * np.ones(len(T_ext))), axis=None)
    #gs = gs[::-1]
    
    
    #solo voy a tomar los valores hasta la temperatura de reheating indicada por
    #la transferencia total de energía del inflaton a la radiación
    #H_inf^2=V_\phi/(3*mpl**2)=(1/(3*mpl**2))*geff*np.pi**2/30 * T**4
    
    mpl= 2.4353e18 # en GeV tal que mpl^2 = 1/(8 \pi G)
    
    H_inf = np.sqrt( V_phi_1_4**4. / ( 3. * mpl**2. ) ) # en GeV
    
    T_rh = V_phi_1_4 * ( 30. / ( np.pi**2. * g[0] ) )**(0.25) # en GeV
    
    # aquí se define mi temperatura inicial aproximada porque uso el g sin
    # saber si todavía le voy a agregar un campo escalar con deg_s grados de
    # libertad, de todas formas solo cambiará de g[0] \simeq 105 a
    # g[0] \simeq 107 asi que voy a recortar todo desde esta T en adelante,
    # al final calculo la Treheating real que pasará a main_integration
    # y veremos que es casi igual
    
    mask = [ T < T_rh ]
    
    mask = tuple(mask)
    
    T = T[mask]
    g = g[mask]
    gs = gs[mask]
    
    #plt.figure()
    #plt.semilogx(T,g,'g.')
    #plt.semilogx(T,gs,'r.')
    
    #%%
    
    # Agrego un campo escalar no masivo con degeneración deg_s. Tiene una
    # temperatura de desacople Tds y una vez desacoplado termina la interacción
    # o solamente pasa que la temperatura evoluciona desacoplada y listo?
    
    # Acá debería agregar un campo escalar masivo pero con ms/H_inf << 1 que además
    # tenga una temperatura de desacople mas chica que su masa. Sería una
    # cold relic. Cuán chica sea la temperatura de desacople en este caso indicará
    # cuanta cantidad de ese campo queda hasta hoy. De todas formas podría ser
    # parte de la materia oscura
    
    def step(x,deltax,valor_intermedio): # acá estamos definiendo un escalón suave 
        smooth = 1/(np.exp(-x/deltax)+1)
        return smooth

#    def step(x,deltax,valor_intermedio): # acá estamos definiendo un escalón suave 
#        sharp = np.heaviside(x,valor_intermedio)
#        return sharp

    if do_cosmo_s:
        
        if do_s_massive: # REVISAAAAR BIEN 
            
            g_SM = g
            gs_SM = gs
            
            gs = gs_SM + deg_s*step(T-ms,0.05*ms,1)
            g = g_SM + deg_s*step(T-ms,0.05*ms,1)
            
            omega_rad0 = omega_gamma0 * g[-1] / 2 # cuenta la radiación no fotónica también
            omega_m0 = 1.0 - omega_rad0 - omega_k0 - omega_lambda0    # 0.3
                    
            a_final = 1             # hoy
            points_number_a = 1e4 + 1 # porque luego le voy a sacar 1
            
            gs0 = gs[-1]
            
            a = a_final * ( gs0 / gs )**(1/3) * T0 / T
            
            a_initial = a[0]
            
            #aeq = omega_rad0 / omega_m0 # con a0=1 y ojo que si no hay energia oscura no lo standard
            
            f_s = omega_gamma0*(gs[-1]/gs)**(4/3)/(omega_m0*a+g/2*omega_gamma0*(gs[-1]/gs)**(4/3))*deg_s*step(T-ms,0.05*ms,1)
        
        elif do_s_nomassive:
            
            g_SM = g
            
            gs_SM =gs
            
            gs_SM_T = uspl(T[::-1],gs_SM[::-1],k=4,s=0)
            
            gs = gs_SM+deg_s*(step(T-Tds,0.05*Tds,0)+(gs_SM/gs_SM_T(Tds))*step(-T+Tds,0.05*Tds,1))
            
            g = g_SM+deg_s*(step(T-Tds,0.05*Tds,0)+(gs_SM/gs_SM_T(Tds))**(4/3)*step(-T+Tds,0.05*Tds,1))
            
            gs_TT = uspl(T[::-1],gs[::-1],k=4,s=0)
            
            omega_rad0 = omega_gamma0 * g[-1] / 2 # cuenta la radiación no fotónica también
            omega_m0 = 1.0 - omega_rad0 - omega_k0 - omega_lambda0    # 0.3
            
            a_final = 1             # hoy
            points_number_a = 1e4 + 1 # porque luego le voy a sacar 1
            
            gs0 = gs[-1]
            
            a = a_final * ( gs0 / gs )**(1/3) * T0 / T
            
            a_initial = a[0]
            
            #aeq = omega_rad0 / omega_m0 # con a0=1 y ojo que si no hay energia oscura no lo standard
            
            Neff_today = (g[-1]/2-1)*(8/7)*(11/4)**(4/3)
            
            if kill_interaction_after_decoupling:
                
                f_s = ( deg_s/2 * omega_gamma0 / ( g * omega_gamma0 / 2 + omega_m0 * a * ( gs / gs[-1] )**(4/3) + omega_lambda0 * a**4 * ( gs / gs[-1] )**(4/3) ) * ( step(T-Tds,0.05*Tds,0) + step(-T+Tds,0.05*Tds,1) * ( gs_SM / gs_SM_T(Tds) )**(4/3)) ) * step(T-Tds,0.05*Tds,1)
            
            else:
                
                f_s = deg_s/2 * omega_gamma0 / ( g * omega_gamma0 / 2 + omega_m0 * a * ( gs / gs[-1] )**(4/3) + omega_lambda0 * a**4 * ( gs / gs[-1] )**(4/3) ) * ( step(T-Tds,0.15*Tds,0) + step(-T+Tds,0.15*Tds,1) * ( gs_SM / gs_SM_T(Tds) )**(4/3))
            
    else:
        
        omega_rad0 = omega_gamma0 * g[-1] / 2 # cuenta la radiación no fotónica también
        omega_m0 = 1.0 - omega_rad0 - omega_k0 - omega_lambda0    # 0.3
        
        a_final = 1             # hoy
        points_number_a = 1e4 + 1 # porque luego le voy a sacar 1
        
        gs0 = gs[-1]
        
        a = a_final * ( gs0 / gs )**(1/3) * T0 / T
        
        a_initial = a[0]
        
        #aeq = omega_rad0 / omega_m0 # con a0=1 y ojo que si no hay energia oscura no lo standard
        
        f_s = np.zeros(len(a))
    
    me = 0.5*1e-3 # en GeV masa del electrón
    
    gs_T = uspl(T[::-1],gs[::-1],k=4,s=0)
    
    f_nu = ( 21/8 * omega_gamma0 / ( g * omega_gamma0 / 2 + omega_m0 * a * ( gs / gs[-1] )**(4/3) + omega_lambda0 * a**4 * ( gs / gs[-1] )**(4/3) )  ) * ( step(T-T_nu,0.05*T_nu,0) + step(-T+T_nu,0.05*T_nu,1) * ( gs / gs_T(T_nu) )**(4/3))
 
        
    #%%
    
    T_a = uspl(a,T,k=4,s=0)
#    a_T = uspl(T,a,k=4,s=0)
    g_a = uspl(a,g,k=4,s=0)
    gs_a = uspl(a,gs,k=4,s=0)
    
    T_a_minus_a_for_ls = uspl(a,T_a(a)-T_ls,s=0)
    
    a_ls = T_a_minus_a_for_ls.roots()
    
    T_a_minus_a_for_nu = uspl(a,T_a(a)-T_nu,s=0)
    
    a_nu = T_a_minus_a_for_nu.roots()
    
    T_a_minus_a_for_eval = uspl(a,T_a(a)-T_eval,s=0)
    
    a_eval = T_a_minus_a_for_eval.roots()
    
    T_a_minus_a_for_ds = uspl(a,T_a(a)-Tds,s=0)
    
    a_ds = T_a_minus_a_for_ds.roots()
    
    f_nu_a = uspl(a,f_nu,k=4,s=0)
    
    f_s_a = uspl(a,f_s,k=4,s=0)
    
    f_s_T = uspl(T[::-1],f_s,k=4,s=0)

    f_nu_T = uspl(T[::-1],f_nu,k=4,s=0)
    
    #%%
    
    def H(x):
        return H0*np.sqrt(omega_m0*x**(-3)+g_a(x)/2*(gs_a(x)/gs0)**(-4/3)*omega_gamma0*x**(-4)+omega_lambda0)
    
    def Feta(x):
        return (H0*np.sqrt(omega_m0*x+g_a(x)/2*(gs_a(x)/gs0)**(-4/3)*omega_gamma0+omega_lambda0*x**4))**(-1) # = (x**2*H(x))**(-1)
    
    a = np.logspace(np.log10(a_initial),np.log10(a_final),points_number_a)
    
    Feta_a = uspl(a,Feta(a),k=4,s=0)
    
    eta = np.zeros(len(a))
    integral=0
    
    for a_i, a_j, n in zip(a[:-1],a[1:],range(len(a)-1)):
        integral += quad(Feta_a,a_i,a_j)[0] #OJO nunca mire los errores en quad
        eta[n+1] = integral
        
    a = a[1::]
    
    eta = eta[1::] # en segundos

    eta_a = uspl(a,eta,s=0)
    
    eta_ls = eta_a(a_ls) # en segundos
    
    eta_nu = eta_a(a_nu) # en segundos
    
    eta_eval = eta_a(a_eval) # en segundos
    
    eta_ds = eta_a(a_ds) # en segundos
    
    a_eta = uspl(eta,a,k=4,s=0)
    
    dadeta_eta = a_eta.derivative() # en Hz
    
    dadeta_a2_eta = uspl(eta,dadeta_eta(eta)/a_eta(eta)**2,k=4,s=0) # en Hz
    
    H_phys = dadeta_a2_eta(eta)
    
    kc_from_horizon_crossing = H_phys * a # en Hz
    
    H_eval = dadeta_a2_eta(eta_eval) # en Hz
    
    H_nu = dadeta_a2_eta(eta_nu) # en Hz
    
    H_ds = dadeta_a2_eta(eta_ds) # en Hz
        
    f_nu = f_nu_a(a)
    
    f_s = f_s_a(a)
    
    T_reheating = a[-1] * ( gs0 / gs[0] )**(1/3) * T0 / a[0] # en GeV

else:
    
    H0 = h * np.float64(3.241e-18) # en Hz
    
    omega_rad0 = omega_gamma0 * geff_today / 2 # cuenta la radiación no fotónica también
    omega_m0 = 1.0 - omega_rad0 - omega_k0 - omega_lambda0    # 0.3
    
    aeq = omega_rad0 / omega_m0 # con a0=1 y ojo que si no hay energia oscura no lo standard
    
    #%%  Cáluclo del factor de escala a(\eta), \eta = tiempo conforme
    
    a_final = 1             # hoy
    points_number_a = 1e4 + 1 #este mas un es porque el primero da 0 y lo voy a sacar despues
    a_initial = 1e-30
    
    def Feta(x):
        return (H0*np.sqrt(omega_m0*x+omega_rad0+omega_lambda0*x**4))**(-1) # = (x**2*H(x))**(-1)
    
    a = np.logspace(np.log10(a_initial),np.log10(a_final),points_number_a)
    
    Feta_a = uspl(a,Feta(a),k=4,s=0)
    
    eta = np.zeros(len(a))
    integral=0
    
    for a_i, a_j, n in zip(a[:-1],a[1:],range(len(a)-1)):
        integral += quad(Feta_a,a_i,a_j)[0] #OJO nunca mire los errores en quad
        eta[n+1] = integral
    
    a = a[1::]
    
    eta = eta[1::]  # en s
    
    eta_a = uspl(a,eta,s=0)  # en s
    
    a_eval = T0/T_eval
    
    eta_eval = eta_a(a_eval) # en s
    
    a_eta = uspl(eta,a,k=4,s=0)
    
    dadeta_eta = a_eta.derivative()  # en Hz
    
    dadeta_a2_eta = uspl(eta,dadeta_eta(eta)/a_eta(eta)**2,k=4,s=0) # en Hz
    
    H_phys = dadeta_a2_eta(eta)
    
    frequency_from_horizon_crossing = H_phys * a / ( 2 * np.pi ) # en Hz

    H_eval = dadeta_a2_eta(eta_eval) # en Hz
    
#%%

np.savetxt('f_s_Tds_0_125',f_s)
np.savetxt('temperature_Tds_0_125',T_a(a))
np.savetxt('geff_Tds_0_125',g_a(a))
np.savetxt('kc_horizon_crossing_Tds_0_125',kc_from_horizon_crossing)

#%%
#from matplotlib import pyplot as plt

#plt.ion()

#plt.figure()
#plt.semilogx(frequency_from_horizon_crossing,g_a(a)/200)
#plt.semilogx(T_a(a),g_a(a)/200)
#plt.semilogx(frequency_from_horizon_crossing,f_s)
#plt.figure(2)
#plt.semilogx(T_a(a),f_s)
#plt.semilogx(T,0.40523/(1 + 8*1e-10/T))
#plt.ylim(0,0.6)
#plt.xlim(8*1e-11,1.5*1e-1)
##plt.semilogx(a,f_nu_a(a))
