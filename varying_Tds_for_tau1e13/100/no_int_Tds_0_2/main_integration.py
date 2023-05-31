     #%%
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import UnivariateSpline as uspl
from scipy.optimize import root as root
from params import h1_initial, h2_initial, h3_initial, tau, step_max_num, u_max_num, do_s_interaction, do_neutrino_damping, do_viscosity, spectrum_not_today
from cosmology import a, eta, eta_nu, eta_eval, f_nu, f_s
from scipy.special import sici as sici

#%%

l = int(sys.argv[1])
i_initial = int(sys.argv[2])
i_final = int(sys.argv[3])

kc_min = np.float64(1e-18)  # en Hz
kc_max = np.float64(1e-3)   # en Hz
points_number_kc = np.float64(int(sys.argv[4])) # 5*1e4

kc_vec = np.logspace(np.log10(kc_min),np.log10(kc_max),points_number_kc)

#%%
if l == 1 :
    np.savetxt('kc_vec',kc_vec)

#%%

h_u0 = np.zeros(len(kc_vec))
h_prime_u0 = np.zeros(len(kc_vec))
spectrum = np.zeros(len(kc_vec))

#%%

# condiciones iniciales
y0 = np.asarray([np.float64(h1_initial), np.float64(h2_initial), np.float64(h3_initial)])

#%%

# función escalón suave o sharp

#def step(x,deltax,valor_intermedio): # acá estamos definiendo un escalón suave 
#    smooth = 1/(np.exp(-x/deltax)+1)
#    return smooth

def step(x,deltax,valor_intermedio): # acá estamos definiendo un escalón suave 
    sharp = np.heaviside(x,valor_intermedio)
    return sharp


# integracion

saving = [int(i_initial+i*(i_final-i_initial)/10) for i in range(1,10)] # para guardado parcial

for i in range(i_initial, i_final):
    
    if spectrum_not_today:
        
        mask = [ eta < eta_eval ]
        
        mask = tuple(mask)
        
        eta = eta[mask]
    
    b = 0.5 #constante que depende de la linealizacion del modelo hydro
    
    kc = kc_vec[i]

    u = eta * kc
    
    u_nu = eta_nu * kc
    
    u_end_num = np.float64(min(u[-1],u_max_num))
    
    u_ini = u[0]
    
    a_u = uspl(u,a,k=4,s=0)
    
    dadu_u = a_u.derivative()
    
    dadu_a_u = uspl(u,dadu_u(u)/a_u(u),k=4,s=0)
    
    if do_s_interaction and do_neutrino_damping:
        
        f_nu_u = uspl(u,f_nu,k=4,s=0)
        
        f_s_u = uspl(u,f_s,k=4,s=0)
        
        viscosity_nu = uspl(u,24*f_nu_u(u)*(dadu_a_u(u))**2.*((np.heaviside((u-u_nu)-1e-2,1))*(1./(8.*(u-u_nu)**4.))*(((u-u_nu)**2.+6.)*(u-u_nu)*np.cos((u-u_nu))+((u-u_nu)**2.-6.)*np.sin((u-u_nu))+((u-u_nu))**4.*sici((u-u_nu))[0])+(np.heaviside(-(u-u_nu)+1e-2,0))*((u-u_nu)/15-(u-u_nu)**3/630+(u-u_nu)**5/37800)),k=4,s=0)
        # en viscosity_nu no se puede usar el smoothstep porque es algo que pega bien así
        
        def step_function_for_nu(x):
            y = step(x-u_nu,0.05*u_nu,0)
            return y
        
        def func(u,y):
            h1, h2, h3 = y
            dydu = [ h2, - h1 - 2. * dadu_a_u(u) * h2 - step_function_for_nu(u) * viscosity_nu(u) * h2 + 6. * f_s_u(u) * (dadu_a_u(u))**2. * 0.53333333 *h3, - 1./(kc*tau) * h3 - b * h2 ]
            return dydu
        
        ui = u_ini
        uf = u_end_num
        
        u_span = np.asarray([ui, uf])
        
        sol = solve_ivp(func,u_span,y0,method='RK45',max_step=step_max_num)
        
    elif do_s_interaction:
        
        if kc < 1 / ( 1000 * tau ) :
            
            def func(u,y):
                h1, h2 = y
                dydu = [ h2, - h1 - 2. * dadu_a_u(u) * h2 ]
                return dydu
            
            y0_aca = y0[0:2]
            
            ui = u_ini
            uf = u_end_num
            
            u_span = np.asarray([ui, uf])
            
            sol = solve_ivp(func,u_span,y0_aca,method='RK45',max_step=step_max_num)

        else:
            
            f_s_u = uspl(u,f_s,k=4,s=0)
            
            def func(u,y):
                h1, h2, h3 = y
                dydu = [ h2, - h1 - 2. * dadu_a_u(u) * h2 + 6. * f_s_u(u) * (dadu_a_u(u))**2. * 0.53333333 * h3, - 1. / ( kc * tau ) * h3 - b * h2 ]
                return dydu
                
            ui = u_ini
            uf = u_end_num
            
            u_span = np.asarray([ui, uf])
            
            sol = solve_ivp(func,u_span,y0,method='RK45',max_step=step_max_num)
    
    elif do_neutrino_damping:
        
        f_nu_u = uspl(u,f_nu,k=4,s=0)
        
        viscosity_nu = uspl(u,24*f_nu_u(u)*(dadu_a_u(u))**2.*((np.heaviside((u-u_nu)-1e-2,1))*(1./(8.*(u-u_nu)**4.))*(((u-u_nu)**2.+6.)*(u-u_nu)*np.cos((u-u_nu))+((u-u_nu)**2.-6.)*np.sin((u-u_nu))+((u-u_nu))**4.*sici((u-u_nu))[0])+(np.heaviside(-(u-u_nu)+1e-2,0))*((u-u_nu)/15-(u-u_nu)**3/630+(u-u_nu)**5/37800)),k=4,s=0)
        # en viscosity_nu no se puede usar el smoothstep porque es algo que pega bien así
        
        def step_function_for_nu(x):
            y = step(x-u_nu,0.05*u_nu,0)
            return y
        
        def func(u,y):
            h1, h2 = y
            dydu = [ h2, - h1 - 2. * dadu_a_u(u) * h2 - step_function_for_nu(u) * viscosity_nu(u) * h2 ]
            return dydu
        
        ui = u_ini
        uf = u_end_num
        
        u_span = np.asarray([ui, uf])
        
        sol = solve_ivp(func,u_span,np.asarray([h1_initial, h2_initial]),method='RK45',max_step=step_max_num)
    
    elif do_viscosity:
        
        f_s_u = uspl(u,f_s,k=4,s=0)
        
        def func(u,y):
            h1, h2 = y
            dydu = [ h2, - h1 - 2. * dadu_a_u(u) * h2 - 6. * f_s_u(u) * (dadu_a_u(u))**2. * ( b * tau * kc ) * h2 ]
            return dydu
        
        ui = u_ini
        uf = u_end_num
        
        u_span = np.asarray([ui, uf])
        
        sol = solve_ivp(func,u_span,np.asarray([h1_initial, h2_initial]),method='RK45',max_step=step_max_num)
        
    else:
        
        def func(u,y):
            h1, h2 = y
            dydu = [ h2, - h1 - 2. * dadu_a_u(u) * h2]
            return dydu
        
        ui = u_ini
        uf = u_end_num
        
        u_span = np.asarray([ui, uf])
        
        sol = solve_ivp(func,u_span,np.asarray([h1_initial, h2_initial]),method='RK45',max_step=step_max_num)
    
    h1_0 = sol.y[0][-1]
    
    h2_0 = sol.y[1][-1]
        
    if u[-1] < u_max_num:
        
        h_prime_u0[i] = h2_0
        
    else:
        
        h_uf = h1_0
        h_prime_uf = h2_0
        
        def WKB_sol_uf(p):
            amp, delta = p
            y = np.zeros(2)
            y = [amp*np.sin(u_end_num+delta)/a_u(u_end_num) - h_uf,
                 amp*(np.cos(u_end_num+delta)/a_u(u_end_num)-np.sin(u_end_num+delta)*dadu_a_u(u_end_num)/a_u(u_end_num)) - h_prime_uf]
            return y
        
        amp_guess = 1.0
        delta_guess = 0.0
        guess = np.asarray([amp_guess, delta_guess])
        
        nsol = root(WKB_sol_uf, guess)
        
        amp = nsol.x[0]
        delta = nsol.x[1]
        
        def h_prime_analytical(u):
            y = amp*(np.cos(u+delta)/a_u(u)-np.sin(u+delta)*dadu_a_u(u)/a_u(u))
            return y
        
        def h_analytical(u):
            y = amp*np.sin(u+delta)/a_u(u)
            return y

        h_prime_u0[i] = h_prime_analytical(u[-1])
        
    if i in saving:
        
        np.savetxt('h_prime_u0_{}'.format(l), h_prime_u0)

#%%
np.savetxt('h_prime_u0_{}'.format(l), h_prime_u0)
