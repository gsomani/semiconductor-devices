import matplotlib.pyplot as plt
import scipy as sp
import GaAs

T=sp.linspace(1,1000,1000)
n_i=GaAs.n_i(T)
plt.grid(True, which="both")
plt.semilogy(T,n_i)
plt.title('Intrinsic carrier concentration variation with Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Intrinsic carrier concentration ($cm^{-3}$)')
plt.show()

T=sp.linspace(1,300,300)
n_0=GaAs.n_0(T)
n_c=0.1*GaAs.n_0(300)
for i in range(300):
    if n_0[i] > n_c :
           break;
T_f=((n_0[i]-n_c)*T[i-1]+(n_c-n_0[i-1])*T[i])/(n_0[i]-n_0[i-1])
print("Freeze out temperature = "+"{:.2f} K".format(T_f)) 
plt.title('Electron concentration variation with Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Electron concentration / $10^{16}$ $cm^{-3}$')
plt.grid(True, which="both")
plt.plot(T,n_0*1e-16)
plt.show()

plt.title('Electron concentration variation with Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Electron concentration ($cm^{-3}$)')
plt.grid(True, which="both")
plt.semilogy(1000/T,n_0)
plt.show()
