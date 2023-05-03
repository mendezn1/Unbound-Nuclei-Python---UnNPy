import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import breitwigner as bw
import functions_edecay as funcs
import lorentz_boost as lb
import os

path = os.getcwd()
#save_directory = '\data\\'
#filename_f = 'boosted_data_fragments.txt'
#filename_n = 'boosted_data_neutrons.txt'

amu = 931.49410242 # atomic mass in MeV/c^2 (NIST)
m_neutron = 1.00866491590 * amu  # neutron amu mass from in MeV/c^2 (AME2020) (0.000000000047 - sigma)
m_electron = 0.000548579909065 * amu # electron amu mass from in MeV/c^2 (NIST)
m_fragment = 51.963213646 * amu  - 20*m_electron # Ca-52 fragment mass in MeV/c^2 (AME2020)(0.000000720 - sigma)
m_parent = 52.968451000 * amu - 20*m_electron # Ca-53 fragment mass in MeV/c^2 (AME2020) (0.000047000 - sigma)
m_beam =  53.963029359 * amu - 21*m_electron # Sc-54 secondary beam mass in MeV/c^2 (AME2020) (0.000015000 - sigma)
m_ex_neutron = 8071.31806/1000  # neutron mass excess in MeV (AME2020) (0.00044 - sigma)
m_ex_fragment = -34266.272/1000  # Ca-52 mass excess in MeV (AME2020) (0.671 - sigma)
m_ex_parent = -29387.707/1000 # Ca-53 mass excess in MeV (AME2020) (43.780 - sigma)
m_ex_beam =  -34437.934/1000 # Sc-54 mass excess in MeV (AME2020) (13.973 - sigma)
s_n_parent = m_ex_fragment+m_ex_neutron-m_ex_parent # MeV | Neutron separation energy of Ca-53 (3.195 MeV (Calculated))

# Differenct experiment parameters
ke_beam = 180*52.9684 # Kinetic energy of the secondary Sc-54 beam
###gamma_beam   = 1+ke_beam/m_beam # Gamma factor for ### NOT USING JUST FOR CALCULATING PURPOSESES
gamma_parent = 1+ke_beam/m_parent
###gamma_fragment = 1+ke_beam/m_fragment ### NOT USING JUST FOR CALCULATING PURPOSESES

beta_parent = np.sqrt(1-(1/gamma_parent**2))

def kinematics(events, eZero, width, angMom, flagGamma, fragMass):
    # Get BW distribution
    bw_neutron_1  = bw.MakeBW(eZero,width,angMom,flagGamma,fragMass)
    bw_pdf_energy = bw_neutron_1[0]
    bw_pdf_prob   = bw_neutron_1[1]

    # Sample from the BW
    rng = np.random.default_rng(seed=12345)
    bw_energy = rng.choice(bw_pdf_energy, size=events, p=bw_pdf_prob/np.sum(bw_pdf_prob))

    # Total energy of parent
    e_parent = m_parent+s_n_parent+bw_energy

    # Solve for momentum in center-of-mass frame
    p_neutron = np.sqrt(((((e_parent**2)+(m_neutron**2)-(m_fragment**2))/(2*e_parent))**2)-(m_neutron**2))

    theta = np.random.uniform(0,180,events)#np.radians(np.repeat(90,events))#np.random.uniform(0,180,events)#np.radians(thta)
    phi   = np.random.uniform(0,360,events)#np.radians(ph)

    vec_f = []
    vec_n = []
    for j in range(0,len(p_neutron)):
        four_vec_f = np.array([np.sqrt(m_fragment**2+p_neutron[j]**2),       # Energy Total
                              p_neutron[j]*np.sin(theta[j])*np.cos(phi[j]),  # px
                              p_neutron[j]*np.sin(theta[j])*np.sin(phi[j]),  # py
                              p_neutron[j]*np.cos(theta[j])])                # pz
        four_vec_n = np.array([np.sqrt(m_neutron**2+p_neutron[j]**2),        # Energy Total
                              -p_neutron[j]*np.sin(theta[j])*np.cos(phi[j]), # px
                              -p_neutron[j]*np.sin(theta[j])*np.sin(phi[j]), # py
                              -p_neutron[j]*np.cos(theta[j])])               # pz
        vec_f.append(four_vec_f)
        vec_n.append(four_vec_n)

    four_vecs_f = vec_f
    four_vecs_n = vec_n

    boosted_vecs_f = []
    for k in four_vecs_f:
        boosted_f = lb.boost(gamma_parent,k,0,0,-1)
        boosted_vecs_f.append(boosted_f)
    boosted_new_vecs_f = boosted_vecs_f

    boosted_vecs_n = []
    for k in four_vecs_n:
        boosted_n = lb.boost(gamma_parent,k,0,0,-1)
        boosted_vecs_n.append(boosted_n)
    boosted_new_vecs_n = boosted_vecs_n

    data = []
    data.append(boosted_new_vecs_f)
    data.append(boosted_new_vecs_n)

    #f_4vec = np.array(boosted_new_vecs_f)
    #n_4vec = np.array(boosted_new_vecs_n)

    #np.savetxt(path+save_directory+filename_f, boosted_new_vecs_f, delimiter=',', fmt='%f')
    #np.savetxt(path+save_directory+filename_n, boosted_new_vecs_n, delimiter=',', fmt='%f')

    return data

#raw = kinematics(int(1e4),0.5,0.5,1,0,52)
#raw

#bw_neutron_1  = bw.MakeBW(0.5,0.5,4,0,52)
#bw_pdf_energy1 = bw_neutron_1[0]
#bw_pdf_prob1   = bw_neutron_1[1]
#bw_neutron_2  = bw.MakeBW(2.5,3.6,2,0,52)
#bw_pdf_energy2 = bw_neutron_2[0]
#bw_pdf_prob2   = bw_neutron_2[1]
#bw_neutron_3  = bw.MakeBW(1.9,1.8,2,0,52)
#bw_pdf_energy3 = bw_neutron_3[0]
#bw_pdf_prob3   = bw_neutron_3[1]
#bw_pdf_tot = bw_pdf_prob1+bw_pdf_prob2+bw_pdf_prob3
#np.sum(bw_pdf_prob3)
#
#plt.plot(bw_pdf_energy1, bw_pdf_prob1, label='9/2+', color='red')
#plt.plot(bw_pdf_energy3, bw_pdf_prob3, label='5/2+', color='blue')
#plt.plot(bw_pdf_energy2, bw_pdf_prob2, label='3/2+', color='green')
#plt.plot(bw_pdf_energy1, bw_pdf_tot,   label='Combo', color='black')
#plt.title('Breit-Wigner Lineshapes')
#plt.xlabel(r'$E_d$ (MeV)')
#plt.ylabel('Probability')
#plt.xlim(0,5)
#plt.ylim(0,0.002)
#plt.legend()
#plt.show()


###kinematics_data_1 = np.array(kinematics(int(1e4),0.5,0.5,1,0,52,0,0)[0])[:,0].flatten()
###kinematics_data_2 = np.array(kinematics(int(1e4),0.5,0.5,1,0,52,0,0)[0])[:,0].flatten()
###np.min(kinematics_data_1)
###np.min(kinematics_data_2)
###np.max(kinematics_data_1)
###np.max(kinematics_data_2)
###plt.hist(kinematics_data_1, bins=100, histtype='step')
###plt.hist(kinematics_data_2, bins=100, histtype='step')
###plt.show()

###### Get BW distribution
#####bw_neutron_1  = bw.MakeBW(0.5,0.5,1,0,52) #eZero,width,angMom,flagGamma,fragMass
#####bw_pdf_energy = bw_neutron_1[0]
#####bw_pdf_prob   = bw_neutron_1[1]
#####
###### Sample from the BW
#####rng = np.random.default_rng(seed=12345)
#####bw_energy = rng.choice(bw_pdf_energy, size=int(1e6), p=bw_pdf_prob/np.sum(bw_pdf_prob))
#####
#####p_neutron = []
#####for bw_e_val in bw_energy:
#####    func = lambda p_neutron : -(m_parent+s_n_parent+bw_e_val) + np.sqrt(m_fragment**2+p_neutron**2) + np.sqrt(m_neutron**2+p_neutron**2)
#####    guess = 30
#####    p_solved = fsolve(func, guess)
#####    p_neutron.append(float(p_solved))
#####
#####e_parent = m_parent+s_n_parent+bw_energy         #          0.044640
#####
#####p_solve = np.sqrt(((((e_parent**2)+(m_neutron**2)-(m_fragment**2))/(2*e_parent))**2)-(m_neutron**2))
#####
#####for i in p_solve:
#####    print(i)
#####
#####for i in np.sort(p_solve):
#####    print(i)
#####
#####p_solve = p_solve[~np.isnan(p_solve)]
#####
#####plt.hist(p_neutron, histtype='step',bins=100, range=[0,100], label='solve', linestyle = 'solid')
#####plt.hist(p_solve  , histtype='step',bins=100, range=[0,100], label='calc', linestyle = (0, (5, 5)))
#####plt.title('CoM Momentum')
#####plt.xlabel('MeV/c')
#####plt.ylabel('Counts')
#####plt.xlim(0,100)
#####plt.legend()
#####plt.show()
#####np.sum(p_neutron)-np.sum(p_solve)
