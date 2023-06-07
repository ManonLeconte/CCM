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

#courants limites de diffusion
ial = F*A*mR*cR
icl = -F*A*mO*cO

E = np.linspace(-3,3,10000)
#Pour les courants anodique/cathodique/totaux en transfert de charge limitant
ia_charge = i0*np.exp(alpha*n*f*(E-Es))
ic_charge = -i0*np.exp(-(1-alpha)*n*f*(E-Es))
i_charge = ia_charge + ic_charge
#Pour la diffusion limitante
id = (ial*np.exp(n*f*(E-Es)) + icl)/(1+np.exp(n*f*(E-Es)))
#Pour un contrôle mixte
it = (i_charge*id)/(i_charge+id)

plt.grid()
plt.title("Courbe i = f(E) pour le couple du fer")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(0,1.3)
plt.ylim(-1,1)
#Pour les courants anodique/cathodique/totaux en transfert de charge limitant
#plt.plot(E,ia_charge,label="courant anodique - charge",color="b")
#plt.plot(E,ic_charge,label="courant cathodique - charge",color="r")
plt.plot(E,i_charge,label="courant total - TE limitant",color="g")
#Pour la diffusion limitante
plt.plot(E,id,label="courant total - TM limitant")
#Pour un contrôle mixte
plt.plot(E,it,label="courant total - mixte")
plt.legend()
plt.show()
