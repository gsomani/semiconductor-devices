import scipy as sp
from scipy import constants

# physical and mathematical constants
m_0=constants.m_e
k=constants.value('Boltzmann constant in eV/K')
h=constants.h
eV=constants.eV
pi=constants.pi

# system constants
N_d=2e16 # per cubic centimeter
E_d=-0.02
m_h=0.5*m_0
m_e=0.05*m_0
g=2

def E_g(T):
    return 1.52-(5.405*10**-4)*(T**2)/(T+ 204)

def n_i(T):
    N_v=2*(2*pi*m_h*k*T*eV/(h*h))**1.5 
    N_c=2*(2*pi*m_e*k*T*eV/(h*h))**1.5    
    return 1e-6*sp.sqrt(N_v*N_c)*sp.exp(-E_g(T)/(2*k*T)) # 1e-6 to calculate n_i per cubic centimeter

def n_0(T):
    N_c=1e-6*2*(2*pi*m_e*k*T*eV/(h*h))**1.5 # 1e-6 to calculate N_c per cubic centimeter
    r=sp.exp(-E_d/(k*T))
    x=(sp.sqrt(1+4*g*r*N_d/N_c)-1)/(2*g*r)
    return N_c*x

