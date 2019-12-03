import matplotlib.pyplot as plt
import scipy as sp
import Mg_GaN as gan

x=sp.linspace(16,19,100)
N_a=10**x  # N_a is per cubic centimeter
T=300
fermi = lambda N : gan.fermi_level(T,N) 
E_f=sp.vectorize(fermi)(N_a)
plt.grid(True, which="both")
plt.semilogx(N_a,E_f) # E_f in meV
plt.title('Fermi level variation with Mg doping')
plt.xlabel('Mg doping ($cm^{-3}$)')
plt.ylabel('Fermi Level ($E_F-E_V$) (meV)')
plt.show()

T=sp.array([5,50,100,300,500,800])
N_a=10**18
fermi = lambda T : gan.fermi_level(T,N_a) 
E_f=sp.vectorize(fermi)(T)
for i in range(6):
    print("Fermi Level at "+str(T[i])+" K = {:.0f} meV".format(E_f[i])) 
plt.grid(True, which="both")
plt.title('Fermi level variation with Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Fermi Level ($E_F-E_V$) (meV)')
plt.plot(T,E_f,marker='o') # E_f in meV
plt.show()

