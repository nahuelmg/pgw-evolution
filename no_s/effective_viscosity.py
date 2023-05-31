import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as uspl
from params import tau
from cosmology import f_s, a, eta
#from scipy.integrate import quad as quad

plt.ion()


#%%

b = 0.5

valor3c2_8pig = 1.607e26 #en kg/m

a_eta = uspl(eta,a,k=4,s=0)

dadeta_eta = a_eta.derivative()  # en Hz

dadeta_a2_eta = uspl(eta,dadeta_eta(eta)/a_eta(eta)**2,k=4,s=0) # en Hz

dadeta_a2 = dadeta_a2_eta(eta)

rho_c = valor3c2_8pig * dadeta_a2**2 

viscosity_phys = 8 * b / 15 * f_s * rho_c * a * tau

z = 1/a-1

#%%
plt.figure(1)
plt.loglog(z, viscosity_phys, label='Effective viscosity')
plt.axhline(1e8,color='C1',label='Bound for viscosity today')
plt.xlabel('Redshift $z$')
plt.ylabel('Effective viscosity [ Pa s ]')
plt.legend()
plt.grid()