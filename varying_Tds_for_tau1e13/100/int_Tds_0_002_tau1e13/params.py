import numpy as np

# parámetros de cosmology.py

h = .7

T0 = 2.3486*1e-13  # en GeV temperatura fotones hoy

spectrum_not_today = False

T_eval = 5*1e-13 # en GeV representa el tiempo final en el cual evaluaremos el espectro
               # debe ser mayor estricto que T0

rdof = True # activa el cambio en los grados de libertad relativistas que
            # que no son s. A s hay que agregarlo a parte abajo.
            # OJOOO si False entonces do_cosmo_s y do_neutrino_damping deben
            # ser False para que tenga sentido. Si no mirar con mucho detalle.

geff_today = np.float64(3.3830836)  # este es solo para cuando rdof = False
                                    # 3.3830836 es el del modelo standard
                                    # con neutrinos no masivos. Esto fija la
                                    # cantidad de radiación.
                                    # Si se quiere comparar con los do_cosmo_s
                                    # hay que ver entonces cuanto da
                                    # g[-1] para ese caso que incluye a s.

omega_gamma0 = 2.4728*1e-5/h**2

omega_lambda0 = 0.7

omega_k0 = 0.

V_phi_1_4 = 1.5e16 # en GeV raíz cuarta del potencial de inflación

T_ls = 2.585e-10 # en GeV desacople de fotones

T_nu = 0.002 # en GeV desacople neutrinos

do_neutrino_damping = False

do_cosmo_s = True #activa el campos s para los grados de libertad relativistas

deg_s = 2 # degeneracion de s, sin degeneración y con part=antipart deg_s=1 (como fotones con una sola polarización)

do_s_massive = False # OJO CON EL MASIVO HAY ERRORES si do_s es true hace que s sea masivo

ms = 1e2 # en GeV masa de s para caso masivo T0 << ms << H_inf

do_s_nomassive = not do_s_massive # si do_s es true hace que s sea no masivo

Tds = 0.002 # en GeV temperatura de desacople de s para caso no masivo

do_s_interaction = True # activa la interacción pGW / tensores s,
                        # si es True do_cosmo_s debe serlo también
                        # si es falso, kill_interaction_after_decoupling no
                        # tiene relevancia

do_viscosity = False

kill_interaction_after_decoupling = False # si True hace que la interacción
                                          # de los tensores y las pGW
                                          # desaparezca luego de desacople.
                                          # si False, la interacción sigue,
                                          # lo que podría pasar es que siga
                                          # pero que tau vaya como 1/Ts y no
                                          # como 1/Tfotones. (Siempre asumimos
                                          # que la masa es mas baja que 3K)
                                          # También aplica para la viscosidad

if do_cosmo_s is False:
    
    do_s_interaction = False

if do_s_interaction or do_neutrino_damping:
    
    do_viscosity = False

#%% parámetros de main_integration.py

# condiciones iniciales y tau

h1_initial = 1  # condición inicial función de transferencia gravitones
h2_initial = 0 # condición inicial derivada de la función de transferencia gravitones
h3_initial = 0 # condición inicial tensores del fluido

reg_r0 = True # si True entonces regulariza con el cut off, si False entonces hace continuación analítica cut off-independent

tau = 1e13 # este es el tiempo de relajación hoy, por ejemplo tau\sim1/g**4 Ts (con la temperatura física hoy)

# parámetros de integración

step_max_num = 1e-2 # paso de integración máximo

u_max_num = 100 # u máximo hasta que el que integra luego pega con las WKB
