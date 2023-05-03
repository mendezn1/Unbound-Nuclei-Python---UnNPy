import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import kinematics as kin

# Parameters for Breit-Wigner
eZero     = 0.5 # MeV
width     = 0.5 # MeV
angMom    = 4   # L
flagGamma = 0   # traditionaly 0
fragMass  = 52  # fragment Mass

# constants definition
c = 299792458 # speed of light (m/s) (NIST)
amu = 931.494102 # atomic mass in MeV/c^2 (NIST)
m_neutron = 1.008664915 * amu  # neutron amu mass from in MeV/c^2 (AME2020)
m_electron = 0.000548579909065 * amu # electron amu mass from in MeV/c^2 (NIST)
m_fragment = 51.963213646 * amu  - 20*m_electron # Ca-52 fragment mass in MeV/c^2 (AME2020)(0.000000720 - sigma)
m_parent = 52.968451000 * amu - 20*m_electron # Ca-53 fragment mass in MeV/c^2 (AME2020) (0.000047000 - sigma)

ke_beam = 180*52.9684 # Kinetic energy of the parent Sc-53 beam
gamma_parent = 1+ke_beam/m_parent

# Define number of events for PDF sampling
events = int(3e3)
events_sci_notat = '{:.0e}'.format(events)

kinematics_data_1 = kin.kinematics(events,eZero,width,angMom,flagGamma,fragMass) # 9/2+ state
kinematics_data_2 = kin.kinematics(events,2.5  ,3.6  ,2     ,0        ,52) # 3/2+ state
kinematics_data_3 = kin.kinematics(events,1.9  ,1.8  ,2     ,0        ,52) # 5/2+ state

fragment_4vector_1 = np.array(kinematics_data_1[0])
fragment_4vector_2 = np.array(kinematics_data_2[0])
fragment_4vector_3 = np.array(kinematics_data_3[0])

neutron_4vector_1  = np.array(kinematics_data_1[1])
neutron_4vector_2  = np.array(kinematics_data_2[1])
neutron_4vector_3  = np.array(kinematics_data_3[1])

e_fragment_1 = fragment_4vector_1[:,0].flatten()
e_fragment_2 = fragment_4vector_2[:,0].flatten()
e_fragment_3 = fragment_4vector_3[:,0].flatten()

p_fragment_1 = fragment_4vector_1[:,1:].reshape(events, 3)
p_fragment_2 = fragment_4vector_2[:,1:].reshape(events, 3)
p_fragment_3 = fragment_4vector_3[:,1:].reshape(events, 3)

e_neutron_1 = neutron_4vector_1[:,0].flatten()
e_neutron_2 = neutron_4vector_2[:,0].flatten()
e_neutron_3 = neutron_4vector_3[:,0].flatten()

p_neutron_1 = neutron_4vector_1[:,1:].reshape(events, 3)
p_neutron_2 = neutron_4vector_2[:,1:].reshape(events, 3)
p_neutron_3 = neutron_4vector_3[:,1:].reshape(events, 3)

theta_1 = []
phi_1 = []
for i in p_neutron_1:
    x = np.rad2deg(np.arccos(i[2]/np.linalg.norm(i)))
    if   i[0]>0:
        phi = np.arctan(i[1]/i[0])
    elif i[0]<0  and i[1]>=0:
        phi = np.arctan(i[1]/i[0]) + np.pi
    elif i[0]<0  and i[1]<0:
        phi = np.arctan(i[1]/i[0]) - np.pi
    elif i[0]==0 and i[1]>0:
        phi = np.pi/2
    elif i[0]==0 and i[1]<0:
        phi = -np.pi/2
    elif i[0]==0 and i[1]==0:
        phi = 0
    theta_1.append(x)
    phi_1.append(np.rad2deg(phi))

theta_2 = []
phi_2 = []
for i in p_neutron_2:
    x = np.rad2deg(np.arccos(i[2]/np.linalg.norm(i)))
    if   i[0]>0:
        phi = np.arctan(i[1]/i[0])
    elif i[0]<0  and i[1]>=0:
        phi = np.arctan(i[1]/i[0]) + np.pi
    elif i[0]<0  and i[1]<0:
        phi = np.arctan(i[1]/i[0]) - np.pi
    elif i[0]==0 and i[1]>0:
        phi = np.pi/2
    elif i[0]==0 and i[1]<0:
        phi = -np.pi/2
    elif i[0]==0 and i[1]==0:
        phi = 0
    theta_2.append(x)
    phi_2.append(np.rad2deg(phi))

theta_3 = []
phi_3 = []
for i in p_neutron_3:
    x = np.rad2deg(np.arccos(i[2]/np.linalg.norm(i)))
    if   i[0]>0:
        phi = np.arctan(i[1]/i[0])
    elif i[0]<0  and i[1]>=0:
        phi = np.arctan(i[1]/i[0]) + np.pi
    elif i[0]<0  and i[1]<0:
        phi = np.arctan(i[1]/i[0]) - np.pi
    elif i[0]==0 and i[1]>0:
        phi = np.pi/2
    elif i[0]==0 and i[1]<0:
        phi = -np.pi/2
    elif i[0]==0 and i[1]==0:
        phi = 0
    theta_3.append(x)
    phi_3.append(np.rad2deg(phi))

good_1 = np.logical_and(theta_1<=np.rad2deg(np.arctan(1/8)), phi_1<=np.rad2deg(np.arctan(0.8/8)))
good_2 = np.logical_and(theta_2<=np.rad2deg(np.arctan(1/8)), phi_2<=np.rad2deg(np.arctan(0.8/8)))
good_3 = np.logical_and(theta_3<=np.rad2deg(np.arctan(1/8)), phi_3<=np.rad2deg(np.arctan(0.8/8)))

e_fragment_1_good = e_fragment_1[good_1]
e_neutron_1_good  = e_neutron_1 [good_1]
p_fragment_1_good = p_fragment_1[good_1]
p_neutron_1_good  = p_neutron_1 [good_1]
e_fragment_2_good = e_fragment_2[good_2]
e_neutron_2_good  = e_neutron_2 [good_2]
p_fragment_2_good = p_fragment_2[good_2]
p_neutron_2_good  = p_neutron_2 [good_2]
e_fragment_3_good = e_fragment_3[good_3]
e_neutron_3_good  = e_neutron_3 [good_3]
p_fragment_3_good = p_fragment_3[good_3]
p_neutron_3_good  = p_neutron_3 [good_3]

#edecay_1 = np.sqrt((m_fragment)**2 + (m_neutron)**2 + 2*(e_fragment_1*e_neutron_1 - np.einsum('ij, ij->i',p_fragment_1,p_neutron_1))) - m_fragment - m_neutron
#edecay_2 = np.sqrt((m_fragment)**2 + (m_neutron)**2 + 2*(e_fragment_2*e_neutron_2 - np.einsum('ij, ij->i',p_fragment_2,p_neutron_2))) - m_fragment - m_neutron
#edecay_3 = np.sqrt((m_fragment)**2 + (m_neutron)**2 + 2*(e_fragment_3*e_neutron_3 - np.einsum('ij, ij->i',p_fragment_3,p_neutron_3))) - m_fragment - m_neutron
#edecay_t = np.append(edecay_1,(edecay_2,edecay_3))

bin = 50
erange = [0,5]

edecay_1_good = np.sqrt((m_fragment)**2 + (m_neutron)**2 + 2*(e_fragment_1_good*e_neutron_1_good - np.einsum('ij, ij->i',p_fragment_1_good,p_neutron_1_good))) - m_fragment - m_neutron
edecay_2_good = np.sqrt((m_fragment)**2 + (m_neutron)**2 + 2*(e_fragment_2_good*e_neutron_2_good - np.einsum('ij, ij->i',p_fragment_2_good,p_neutron_2_good))) - m_fragment - m_neutron
edecay_3_good = np.sqrt((m_fragment)**2 + (m_neutron)**2 + 2*(e_fragment_3_good*e_neutron_3_good - np.einsum('ij, ij->i',p_fragment_3_good,p_neutron_3_good))) - m_fragment - m_neutron
edecay_t_good = np.concatenate((edecay_1_good,edecay_2_good,edecay_3_good))

#plt.hist(edecay_1_good, bins=bin, range=erange, histtype='step', color='red'  , label=r'$J^{\pi}$ = 9/2+')
#plt.hist(edecay_2_good, bins=bin, range=erange, histtype='step', color='green', label=r'$J^{\pi}$ = 3/2+')
#plt.hist(edecay_3_good, bins=bin, range=erange, histtype='step', color='blue' , label=r'$J^{\pi}$ = 5/2+')
#plt.hist(edecay_t_good, bins=bin, range=erange, histtype='step', color='black', label=r'Total')
#plt.title(f"Events = {events_sci_notat} | fragMass = {fragMass} | Bins = {(erange[1]-erange[0])/bin} MeV")
#plt.xlabel('MeV')
#plt.ylabel('Counts')
#plt.legend()
#plt.xlim(0,5)
#plt.show()

plt.hist(edecay_1_good, bins=bin, range=erange, histtype='step', color='red'  , label=r'$J^{\pi}$ = 9/2+')
plt.hist(edecay_2_good, bins=bin, range=erange, histtype='step', color='green', label=r'$J^{\pi}$ = 3/2+')
plt.hist(edecay_3_good, bins=bin, range=erange, histtype='step', color='blue' , label=r'$J^{\pi}$ = 5/2+')
plt.hist(edecay_t_good, bins=bin, range=erange, histtype='step', color='black', label=r'Total')
plt.title(r"$^{53}$ Ca Edecay | "+f"Events = {events}")
plt.xlabel('MeV')
plt.ylabel(f'Counts / {(erange[1]-erange[0])/bin} MeV')
plt.legend()
plt.xlim(0,5)
plt.show()

#plt.hist(edecay_1, bins=bin, range=erange, histtype='step', label=r'$J^{\pi}$ = 9/2+', color='red')
#plt.hist(edecay_2, bins=bin, range=erange, histtype='step', label=r'$J^{\pi}$ = 3/2+', color='blue')
#plt.hist(edecay_3, bins=bin, range=erange, histtype='step', label=r'$J^{\pi}$ = 5/2+', color='green')
#plt.hist(edecay_t, bins=bin, range=erange, histtype='step', label=r'Total', color='black')
#plt.title(f"Events = {events_sci_notat} | fragMass = {fragMass} | Bins = {(erange[1]-erange[0])/bin} MeV")
#plt.xlabel('MeV')
#plt.ylabel('Counts')
#plt.legend()
#plt.xlim(0,5)
#plt.show()
