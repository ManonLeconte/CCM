import numpy as np
import matplotlib.pyplot as plt

# i0 = 1
i0 = 0.1
alpha = 0.5
n = 1
F = 96500
R = 8.314
T = 298
f = F/(R*T)
Es = 0.7
A = 1e-4
DO=6e-6
DR=6e-6
delta=1e-6
cR=0.01
cO=0.01

mO=DO/delta
mR=DR/delta

E = np.linspace(-3,3,10000)

#potentiel de Nernst
E_N=Es+1/(n*f)*np.log(cO/cR)

#courants limites de diffusion
ial = F*A*mR*cR
icl = -F*A*mO*cO

#courant total
id = (ial*np.exp(n*f*(E-E_N)) + icl)/(1+np.exp(n*f*(E-E_N)))


plt.grid()
plt.title("Courbe i = f(E) pour le couple du fer")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(-1,2)
plt.ylim(-1,1)

#Pour la diffusion limitante
plt.plot(E,id,label="courant total - Limitation par le transport de matière")

plt.legend()
plt.show()
