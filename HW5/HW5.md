
# NE 205: Semiconductor Devices and IC Technology

## Homework 5

### Name : Gaurav Somani
### SR No. : 16082


```sos
import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
import math
from operator import sub

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'
    
MU='\u03BC'    
approx= '\u2248'
exponent_2= '\u00b2'
epsilon_0 = 0.01*constants.epsilon_0 # multiply by 0.01 to covert to F/cm 
q=constants.e
k=constants.k

T=300

Vth=k*T/q
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)
epsilon_r_si=11.7
epsilon_r_sio2=3.8
ni=1e10
epsilon_si=epsilon_r_si*epsilon_0
epsilon_sio2=epsilon_r_sio2*epsilon_0
a0=epsilon_si*q
a_si=Vth*epsilon_si/q
```

\begin{equation}
\begin{aligned}
        \boxed {
            T = 300\ K \\[10pt]
            (\epsilon_r)_{Si} = 11.7 \\[10pt]
            (\epsilon_r)_{SiO_2} = 3.8 \\[10pt]
            \epsilon_{Si} = \epsilon_0(\epsilon_r)_{Si} \\[10pt]
            \epsilon_{SiO_2} = \epsilon_0(\epsilon_r)_{SiO_2} \\[10pt]
            n_i = 10^{10} cm^{-3}  
         }   
        \end{aligned}
          \label{A}
    \tag{A}
\end{equation}

### 1)

\begin{equation}
\begin{aligned}
        \boxed {
            \sigma = 50\ \mu m  \\[10pt]
            t_{ox}= 50\ nm
         }   
        \end{aligned}
        \label{pr1}
    \tag{1}
\end{equation}


\begin{equation}
\begin{aligned}
        \boxed {
            V_{FB} = - \frac{q\ \sigma}{\epsilon_{SiO_2}} d \\[10pt]
            \left\lvert (V_{FB})_{max} \right\rvert = \frac{q\ \sigma\ t_{ox}}{\epsilon_{SiO_2}}\\[10pt]
         }   
        \end{aligned}
        \label{eq1.1}
    \tag{1.1}
\end{equation}

where $d$ represents the distance of center of oxide charge from $metal-SiO_2$ interface.


```sos
sigma=1e11
t_ox=5e-6 # 50 nm
Vfb_max=q*sigma*t_ox/epsilon_sio2
Vfb_center=Vfb_max/2

Vfb_m = round(Vfb_max,2)
Vfb_c = round(Vfb_center,2)
```


```sos
%expand {{ }} 

Using ([A](#mjx-eqn-A)) and ([1](#mjx-eqn-pr1)), $\left\lvert (V_{FB})_{max} \right\rvert = -${{Vfb_m}} $V$

Using ([1.1](#mjx-eqn-eq1.1)) and $\left\lvert (V_{FB})_{max} \right\rvert$, $V_{FB}$ is calculated. 

#### a)  
$d = t_{ox}$
    $\implies V_{FB} = -\left\lvert (V_{FB})_{max} \right\rvert = -${{Vfb_m}} $V$

#### b) 

$d = \frac{t_{ox}}{2}$
    $\implies V_{FB} = -\left\lvert \frac{(V_{FB})_{max}}{2} \right\rvert = -${{Vfb_c}} $V$

#### c) 

$d = 0$
    $\implies V_{FB} = 0$
    
#### d) 

$d = \frac{t_{ox}}{2}$
    $\implies V_{FB} = -\left\lvert \frac{(V_{FB})_{max}}{2} \right\rvert = -${{Vfb_c}} $V$
```


Using ([A](#mjx-eqn-A)) and ([1](#mjx-eqn-pr1)), $\left\lvert (V_{FB})_{max} \right\rvert = -$0.24 $V$

Using ([1.1](#mjx-eqn-eq1.1)) and $\left\lvert (V_{FB})_{max} \right\rvert$, $V_{FB}$ is calculated. 

#### a)  
$d = t_{ox}$
    $\implies V_{FB} = -\left\lvert (V_{FB})_{max} \right\rvert = -$0.24 $V$

#### b) 

$d = \frac{t_{ox}}{2}$
    $\implies V_{FB} = -\left\lvert \frac{(V_{FB})_{max}}{2} \right\rvert = -$0.12 $V$

#### c) 

$d = 0$
    $\implies V_{FB} = 0$
    
#### d) 

$d = \frac{t_{ox}}{2}$
    $\implies V_{FB} = -\left\lvert \frac{(V_{FB})_{max}}{2} \right\rvert = -$0.12 $V$



### 2)

\begin{equation}
\begin{aligned}
        \boxed {
             W = 50\ \mu m  \\[10pt]
            L = 1.5\ \mu m \\[10pt]
            \mu = 600\ cm^2/Vs \\[10pt]    
            C_{min} = 2.9 \times 10^{-8} F/cm^2 \\[10pt]
            C_{max} = 1.72 \times 10^{-7} F/cm^2  
         }   
        \end{aligned}
        \label{pr2}
    \tag{2}
\end{equation}

\begin{equation}
\begin{aligned}
        \boxed {
            t_{ox} = \frac{\epsilon_{SiO_2}}{C_{ox}} \\[10pt]
            (W_D)_{max} = \epsilon_{Si}(\frac{1}{C_{min}} - \frac{1}{C_{max}}) \\[10pt]
            \frac{ ln\frac{N_A}{n_i}}{N_A}= \frac{q^2\ (W_D)_{max}^2}{4\ kT \epsilon_{Si}}\\[10pt]
         }   
        \end{aligned}
        \label{eq2.1}
    \tag{2.1}
\end{equation}

\begin{equation}
\begin{aligned}
        \boxed {
         (I_D)_{sat} = \frac{W}{L}\frac{\mu C_{max} (V{gs}-V_{th})^2}{2}
         }   
        \end{aligned}
        \label{eq2.2}
    \tag{2.2}
\end{equation}


```sos
c_min=2.9e-8
c_max=1.72e-7
t_ox=epsilon_sio2/c_max
wd_max = epsilon_si/c_min - epsilon_si/c_max

def N_A(wd_max):
    Na_min=math.exp(1)*ni
    W=2*(a_si*math.log(Na_min/ni)/Na_min)**0.5
    while(W>wd_max):
        Na_min=Na_min*2
        W=2*(a_si*math.log(Na_min/ni)/Na_min)**0.5
    
    Na_max=Na_min
    Na_min=Na_max/2
    while(Na_max/Na_min-1 > 1e-6):
        Na=0.5*(Na_min+Na_max)
        W=2*(a_si*math.log(Na/ni)/Na)**0.5
        if(W>wd_max):
            Na_min=Na
        else:
            Na_max=Na
    return Na     

Na=N_A(wd_max)

L=1.5e-4
W=5e-3
Vgt=1.5
mu=600
Id_sat=0.5*(mu*c_max*W/L)*Vgt*Vgt

t_oxide=round(t_ox*1e7,1)
Wd_maximum=round(wd_max*1e7)
Na_sub=round(Na*1e-16,2)
Id_s=round(Id_sat*1e3,2)
```


```sos
%expand {{ }} 

#### a) 

Using ([A](#mjx-eqn-A)),([2](#mjx-eqn-pr2)) and ([2.1](#mjx-eqn-eq2.1)),

$t_{ox} = $ {{t_oxide}} $nm$

#### b) 

Using ([A](#mjx-eqn-A)),([2](#mjx-eqn-pr2)) and ([2.1](#mjx-eqn-eq2.1)),

$(W_D)_{max} = $ {{Wd_maximum}} $nm$
        
$\implies N_A = $ {{Na_sub}} $\times 10^{16} cm^{-3}$
        
c) Using ([2](#mjx-eqn-pr2)) and ([2.2](#mjx-eqn-eq2.2)),

$(I_D)_{sat} = $ {{Id_s}} $mA$
```


#### a) 

Using ([A](#mjx-eqn-A)),([2](#mjx-eqn-pr2)) and ([2.1](#mjx-eqn-eq2.1)),

$t_{ox} = $ 19.6 $nm$

#### b) 

Using ([A](#mjx-eqn-A)),([2](#mjx-eqn-pr2)) and ([2.1](#mjx-eqn-eq2.1)),

$(W_D)_{max} = $ 297 $nm$
        
$\implies N_A = $ 1.05 $\times 10^{16} cm^{-3}$
        
c) Using ([2](#mjx-eqn-pr2)) and ([2.2](#mjx-eqn-eq2.2)),

$(I_D)_{sat} = $ 3.87 $mA$



### 3)

\begin{equation}
\begin{aligned}
        \boxed {
            \phi_{ms} = 1\ V  \\[10pt]
            t_{ox} = 50\ nm \\[10pt]
            N_A = 10^{13} cm^{-3}  
         }   
        \end{aligned}
        \label{pr3}
    \tag{3}
\end{equation}

For $\psi_s = 2 \psi_B = 2(E_i – E_F)_{bulk}$ ,

\begin{equation}
\begin{aligned}
        \boxed {
            V_{th} = \phi_{ms} + \frac{Q_s}{C_{ox}} + 2\psi_B  \\[10pt]
            \psi_B = \frac{kT}{q} ln\frac{N_A}{n_i}\\[10pt]
            Q_s \approx 2\sqrt{q\ \epsilon_{Si} N_A \psi_B} \\[10pt]
         }   
        \end{aligned}
        \label{eq3.1}
    \tag{3.1}
\end{equation}

For $n_{inv}=10^{16} cm^{-3}$ ,

\begin{equation}
\begin{aligned}
        \boxed {
            V_{th} = \phi_{ms} + \frac{Q_s}{C_{ox}} + \psi_s  \\[10pt]
            \psi_s = \frac{kT}{q}ln\frac{n_{inv}N_A }{n_i^2}\\[10pt]
            Q_s \approx \sqrt{2\ \epsilon_{Si}( N_A\ (q\psi_s - kT) + \frac{k T\ n_i^2}{N_A}\ e^{\frac{q\ \psi_s}{kT}})}\\[10pt]
         }   
        \end{aligned}
        \label{eq3.2}
    \tag{3.2}
\end{equation}


```sos
Na=1e13
n=1e16
t_ox=5e-6
phi_ms=1
C_ox=epsilon_sio2/t_ox

psi_b=Vth*math.log(Na/ni)
Qd_max=2*(a0*Na*psi_b)**0.5
Vth_inv = phi_ms + Qd_max/C_ox + 2*psi_b

n_inv=1e16
n_p0=(ni*ni)/Na
psi_ninv=Vth*(math.log(Na/ni)+math.log(n_inv/ni))
Qs_ninv=(2*a0*(Na*(psi_ninv-Vth) +  Vth*ni*ni*math.exp(psi_ninv/Vth)/Na))**0.5
Vth_ninv=phi_ms + Qs_ninv/C_ox + psi_ninv

Vt_conventionl = round(Vth_inv,2)
Vt_new = round(Vth_ninv,2)
```


```sos
%expand {{ }} 

For $\psi_s = 2 \psi_B = 2(E_i – E_F)_{bulk}$ , using ([A](#mjx-eqn-A)),([3](#mjx-eqn-pr3)) and ([3.1](#mjx-eqn-eq3.1)),
 $V_{T} = $ {{Vt_conventionl}} $V$.

For $n_{inv}=10^{16} cm^{-3}$ , using ([A](#mjx-eqn-A)),([3](#mjx-eqn-pr3)) and ([3.2](#mjx-eqn-eq3.2)), $V_{T} = $ {{Vt_new}} $V$.
```


For $\psi_s = 2 \psi_B = 2(E_i – E_F)_{bulk}$ , using ([A](#mjx-eqn-A)),([3](#mjx-eqn-pr3)) and ([3.1](#mjx-eqn-eq3.1)),
 $V_{T} = $ 1.37 $V$.

For $n_{inv}=10^{16} cm^{-3}$ , using ([A](#mjx-eqn-A)),([3](#mjx-eqn-pr3)) and ([3.2](#mjx-eqn-eq3.2)), $V_{T} = $ 1.67 $V$.



### 4)

\begin{equation}
\begin{aligned}
        \boxed {
            W = 20\ \mu m  \\[10pt]
            L = 2\ \mu m \\[10pt]
            V_{FB} = -1\ V  \\[10pt]
            V_{DS} = 1\ V  \\[10pt]
            \mu_e = 500\ cm^2/Vs \\[10pt]
            \mu_h = 100\ cm^2/Vs \\[10pt]
            N_A = 10^{16}\ cm^{-3} \\[10pt]
            t_{ox} = 50\ nm         
            }   
        \end{aligned}
        \label{pr4}
    \tag{4}
\end{equation}

\begin{equation}
\begin{aligned}
        \boxed {
            \sigma_{FB} = q\ \mu_h p = q\ \mu_h N_A \\[10pt]
            \sigma_{inv} = q\ \mu_e n = q\ \mu_e N_A
         }   
        \end{aligned}
        \label{eq4.1}
    \tag{4.1}
\end{equation}

\begin{equation}
\begin{aligned}
        \boxed {
           \psi_B = \frac{kT}{q} ln\frac{N_A}{n_i}\\[10pt]
            V_{th} = V_{FB} + \frac{2\sqrt{q\ \epsilon_{Si} N_A \psi_B}}{C_{ox}} + 2 \psi_B  \\[10pt]
            V_{GS} = V_{th} + 0.5\ V  \\[10pt]
            \implies V_{GS}- V_{FB} = \frac{2\sqrt{q\ \epsilon_{Si} N_A \psi_B}}{C_{ox}} + 2 \psi_B + 0.5\ V \\[10pt]
            (p_s)_{source} = N_A e^{-\frac{q\ (\psi_s)_{source}}{kT}}\\[10pt]
            (n_s)_{source} = \frac{n_i^2}{(p_s)_{source}} \\[10pt]
            (p_s)_{drain} = N_A e^{-\frac{q\ (\psi_s)_{drain}}{kT}}\\[10pt]
            (n_s)_{drain} = \frac{n_i^2}{(p_s)_{drain}} e^{-\frac{q\ V_{DS}}{kT}} \\[10pt] 
               \frac{(Q_s)_{drain}}{C_{ox}} + (\psi_s)_{drain} = V_{GS} - V_{FB} \\[10pt] 
            (Q_s)_{drain} \approx \sqrt{2\ \epsilon_{Si} N_A\ (q(\ (\psi_s)_{drain}) - kT)} \\[10pt]
                        \frac{(Q_s)_{source}}{C_{ox}} + (\psi_s)_{source} = V_{GS}- V_{FB} \\[10pt] 
               (Q_s)_{source} \approx \sqrt{2\ \epsilon_{Si}( N_A\ (q(\psi_s)_{drain} - kT) + \frac{k T\ n_i^2}{N_A}\ e^{\frac{q\ (\psi_s)_{source}}{kT}})} \\[10pt]
         }   
        \end{aligned}
        \label{eq4.2}
    \tag{4.2}
\end{equation}

\begin{equation}
\begin{aligned}
        \boxed {
            J = \sigma_{FB} \frac{V_{DS}}{L} = q\ \mu_h N_A \frac{V_{DS}}{L}
         }   
        \end{aligned}
        \label{eq4.3}
    \tag{4.3}
\end{equation}


```sos
W=2e-3
L=2e-4
Vfb=-1
mu_e=500
mu_h=100
Na=1e16
Vds=1
charge=q*Na
sigma_fb=charge*mu_h
sigma_inv=charge*mu_e
psi_b=Vth*math.log(Na/ni)
Qth=2*(a0*Na*psi_b)**0.5
Vth_fb = Qth/C_ox + 2*psi_b
Vgs_fb = Vth_fb + 0.5

def psi_source(Vg):
    psi_min=0
    psi_max=Vg
    psi=0.5*(psi_max+psi_min)
    while(psi_max-psi_min>1e-3):
        psi=0.5*(psi_max+psi_min)
        Qs = (2*a0*(Na*(psi-Vth) +  Vth*ni*ni*math.exp(psi/Vth)/Na))**0.5
        V = Qs/C_ox + psi
        if(V>Vg):
            psi_max=psi
        else:
            psi_min=psi
    return psi        
        
def psi_drain(Vg):
    psi_min=0
    psi_max=Vg
    psi=0.5*(psi_max+psi_min)
    while(psi_max-psi_min>1e-3):
        psi=0.5*(psi_max+psi_min)
        Qs = (2*a0*Na*(psi-Vth))**0.5
        V = Qs/C_ox + psi
        if(V>Vg):
            psi_max=psi
        else:
            psi_min=psi
    return psi

psi_d=psi_drain(Vgs_fb)
psi_s=psi_source(Vgs_fb)
p_s=Na*math.exp(-psi_s/Vth)
p_d=Na*math.exp(-psi_d/Vth)
n_s=ni*ni/p_s
n_d=ni*ni*math.exp(-Vds/Vth)/p_d

J=1e-3*sigma_fb*Vds/L # in kA per square centimeter

conductivity_fb=round(sigma_fb,2)
conductivity_inv=round(sigma_inv,2)
J_ch=round(J,2)
p_source=round(p_s*1e-2,2)
p_drain=round(p_d*1e2,2)
n_source=round(n_s*1e-17,2)
n_drain=round(n_d*1e-5,2)
```


```sos
%expand {{ }} 

#### a)

Using ([4](#mjx-eqn-pr4)) and ([4.1](#mjx-eqn-eq4.1)),

(Conductivity under flat band) $\sigma_{FB} = $ {{conductivity_fb}} $(\Omega\ cm)^{-1}$ ;

(Conductivity at inversion) $\sigma_{inv} = $ {{conductivity_inv}} $(\Omega\ cm)^{-1}$

#### b) 

Using ([A](#mjx-eqn-A)),([4](#mjx-eqn-pr4)) and ([4.2](#mjx-eqn-eq4.2)),

$(n_s)_{source} \approx $ {{n_source}} $\times 10^{17} cm^{-3}$; 
$(p_s)_{source} \approx $ {{p_source}} $\times 10^{2} cm^{-3}$;

$(n_s)_{drain} \approx $ {{n_drain}} $\times 10^{5} cm^{-3}$;
$(p_s)_{drain} \approx $ {{p_drain}} $\times 10^{-2} cm^{-3}$

where $n_s$ represents electron density at the $Si-SiO_2$ interface and $p_s$ represents hole density at the $Si-SiO_2$ interface. 
        
#### c)  
Using ([4](#mjx-eqn-pr4)) and ([4.3](#mjx-eqn-eq4.3)),

$J = $ {{J_ch}} $kA/cm^2$
```


#### a)

Using ([4](#mjx-eqn-pr4)) and ([4.1](#mjx-eqn-eq4.1)),

(Conductivity under flat band) $\sigma_{FB} = $ 0.16 $(\Omega\ cm)^{-1}$ ;

(Conductivity at inversion) $\sigma_{inv} = $ 0.8 $(\Omega\ cm)^{-1}$

#### b) 

Using ([A](#mjx-eqn-A)),([4](#mjx-eqn-pr4)) and ([4.2](#mjx-eqn-eq4.2)),

$(n_s)_{source} \approx $ 3.63 $\times 10^{17} cm^{-3}$; 
$(p_s)_{source} \approx $ 2.76 $\times 10^{2} cm^{-3}$;

$(n_s)_{drain} \approx $ 1.21 $\times 10^{5} cm^{-3}$;
$(p_s)_{drain} \approx $ 1.31 $\times 10^{-2} cm^{-3}$

where $n_s$ represents electron density at the $Si-SiO_2$ interface and $p_s$ represents hole density at the $Si-SiO_2$ interface. 
        
#### c)  
Using ([4](#mjx-eqn-pr4)) and ([4.3](#mjx-eqn-eq4.3)),

$J = $ 0.8 $kA/cm^2$


