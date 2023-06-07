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

# surtension = np.linspace(-0.2,0.2,10000)
# ia = i0*np.exp(alpha*n*f*surtension)
# ic = -i0*np.exp((1-alpha)*n*f*surtension)
E = np.linspace(-3,3,10000)
ia = i0*np.exp(alpha*n*f*(E-Es))
ic = -i0*np.exp(-(1-alpha)*n*f*(E-Es))
i = ia + ic

# plt.xlim(-0.1,0.1)
plt.grid()
plt.title("Courbe i = f(E) pour le couple du fer")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(-0.5,1.5)
plt.ylim(-2,2)
# plt.plot(surtension,ia,label="ia",color="b")
# plt.plot(surtension,ic,label="ic",color="r")
# plt.plot(surtension,i,label="i",color="g")
plt.plot(E,ia,label="courant anodique",color="b")
plt.plot(E,ic,label="courant cathodique",color="r")
plt.plot(E,i,label="courant total",color="g")
plt.legend()
plt.show()
