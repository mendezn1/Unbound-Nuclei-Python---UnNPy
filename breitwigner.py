import numpy as np
from scipy import special
from scipy import integrate

a = 1e-10 # lower limit in MeV (can be adjusted)
b = 100 # Upper Limit in MeV (can be adjusted)
n = round((b-a)*1000) # histogram bins
hbarc = 197.3269631 # MeV*fm
nmass = 939.565346 #MeV
amu = 931.494028 # NIST
nucRad = 1.4 # fm
const = 1/3.14159265 # constant for Breit-Wigner
energy = np.linspace(a,b,n+1)
eN = np.linspace(a,b,n+1)[:-1]+(b-a)/(n*2)

def MakeBW(eZero, width, angMom, flagGamma, fragMass):
    #print(f"AsymBW: E0={eZero} Width={width} L={angMom} flagGamma={flagGamma} fragMass={fragMass}\n")

    redMass = (((nmass/amu)*fragMass)/(nmass/amu + fragMass))*amu # reduced mass
    k = np.sqrt(2 * redMass)/hbarc # k in term os square-root(Energy)
    radius = nucRad * (fragMass**(1/3)) + ((nmass/amu)**(1/3)) #Nuclear Radius A1^(1/3) + A2^(1/3)
    x0 = k*radius*np.sqrt(eZero)
    jl0 = x0 * special.spherical_jn(angMom,x0) # Regular spherical Bessel Function of kind/order angMom @ x0
    nl0 = x0 * special.spherical_yn(angMom,x0) # Irregular spherical Bessel Function of kind/order angMom @ x0
    penL0 = x0/(nl0*nl0 + jl0*jl0) # Penetrability as a function of E,W,L

    energy = np.linspace(a,b,n+1) # bounds for integration
    eN = np.linspace(a,b,n+1)[:-1]+(b-a)/(n*2) # get bin centers
    ###x = k * radius * np.sqrt(eN) # X = k*R with energy included
    x = k * radius * np.sqrt(energy) # X = k*R with energy included
    jl = x * special.spherical_jn(angMom,x) # Regular spherical Bessel Function  @ x
    nl = x * special.spherical_yn(angMom,x) # Irregular spherical Bessel Function  @ x
    penL = x / (nl*nl + jl*jl) # Penetrability as a function of E,W,L
    bound = 0 # Define the boundary condition
    shift = 0 # Define the shift function

    #zwk, check this for angMom==0
    if angMom == 0:
      shift = 0
      bound = 0
      penL = 1
      penL0 = 1

    if angMom == 1:
      shift = - 1 / ( 1 + x**2) # From lane and thomas, i.e. x*(F'*F+G'*G)/(F*F+G*G)
      bound = -1/ ( 1 + x0**2)  # Set boudary so shift = 0 at eigen value (i.e. Decay energy)

    if angMom == 2 : # Same as above for l = 2
      shift = -3*(x**2+6)/(x**4+9+3*x**2)
      bound = -3*(x0**2+6)/(x0**4+9+3*x0**2)

    redGamma = 0
    delta = 0
    gamma = 0
    obsGamma = 0

    if flagGamma == 1 : # For use with Gflag - reduced width
      redGamma = width
      delta = -( shift - bound ) * redGamma * redGamma # Full shift function
      gamma = 2.0 * redGamma * redGamma * penL # Width
    if flagGamma == 0 :# default --observed width CRH 5/12/08
      obsGamma = width
      delta = - ( shift - bound)* (obsGamma / (2 * penL0))
      gamma = (obsGamma / penL0)*penL

    ###bw = const * ((gamma/2)/((eZero + delta - eN)*(eZero + delta - eN) + (gamma*gamma)/4))
    bw = const * ((gamma/2)/((eZero + delta - energy)*(eZero + delta - energy) + (gamma*gamma)/4))

    probability = np.array([])
    for i in range(0,len(eN)):
        bw_range = [bw[i],bw[i+1]]
        energy_range = [energy[i],energy[i+1]]
        holder_prob = integrate.trapezoid(bw_range,energy_range)
        probability = np.append(probability, holder_prob)

    return np.array([eN,probability]) #returns [energy(MeV),Probability]

#data = MakeBW(0.5,0.5,4,0,52)
#np.sum(data[1])

#### Plotting for checks
###import matplotlib.pyplot as plt
###
###data1 = MakeBW(0.5,0.5,1,0,52)
###data2 = MakeBW(0.5,0.5,1,0,52)
###data3 = MakeBW(0.5,0.5,1,0,52)
###data4 = MakeBW(0.5,0.5,1,0,52)
###data5 = MakeBW(0.5,0.5,1,0,52)
###
###nrgy1 = data1[0]
###prob1 = data1[1]
###nrgy2 = data2[0]
###prob2 = data2[1]
###nrgy3 = data3[0]
###prob3 = data3[1]
###nrgy4 = data4[0]
###prob4 = data4[1]
###nrgy5 = data5[0]
###prob5 = data5[1]
###
####len(nrgy)
####len(prob)
####
####np.sum(prob1)
###
###plt.scatter(nrgy1, prob1, marker='.')
###plt.scatter(nrgy2, prob2, marker='.')
###plt.scatter(nrgy3, prob3, marker='.')
###plt.scatter(nrgy4, prob4, marker='.')
###plt.scatter(nrgy5, prob5, marker='.')
###plt.xlim(0.47,0.48)
###plt.ylim(0.001312,0.001322)
###plt.title('Breit-Wigner Probabilities')
###plt.xlabel('MeV')
###plt.ylabel('Probability')
###plt.show()
