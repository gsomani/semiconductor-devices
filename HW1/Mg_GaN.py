import scipy as sp
from scipy import constants

# physical and mathematical constants
m_0=constants.m_e
k=constants.value('Boltzmann constant in eV/K')
h=constants.h
eV=constants.eV
pi=constants.pi

# system constants
m_h=0.6*m_0
g=4
E_a=0.2

def fermi_level(T,N_a):
    N_v=1e-6*2*(2*pi*m_h*k*T*eV/(h*h))**1.5  # 1e-6 to calculate N_v per cubic centimeter
    r=sp.exp(E_a/(k*T))
    x=(sp.sqrt(1+4*g*r*N_a/N_v)-1)/(2*g*r)
    return -k*T*sp.log(x)*1e3 # 1e3 to calculate fermi level in meV



