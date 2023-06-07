import numpy as np
import matplotlib.pyplot as plt


###Définitions générales

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

def id(potentiel,courant1,courant2):#expression générale du courant
    return (courant1*np.exp(n*f*(potentiel-Es)) + courant2)/(1+np.exp(n*f*(potentiel-Es)))

def idv(courant1,courant2):#expression du courant de demi-vague
    return (courant1+courant2)/2

def Edv(courant,potentiel,value):#recherche du potentiel  de demi-vague
    index=(np.abs(courant-value)).argmin()
    return potentiel[index]

###Les deux partenaires du couple sont présents en quantités égales

#courants limites de diffusion
ial1 = F*A*mR*cR
icl1 = -F*A*mO*cO

#courant total
id1 = id(E,ial1,icl1)
idv1=idv(ial1,icl1)
Edv1=Edv(id1,E,idv1)

###Les deux partenaires du couple sont présents en quantités différentes

#courants limites de diffusion
ial2 = F*A*mR*cR
icl2 = -2*F*A*mO*cO

#courant total
id2 = id(E,ial2,icl2)
idv2=idv(ial2,icl2)
Edv2=Edv(id2,E,idv2)

###Un seul partenaire est présent

#courants limites de diffusion
ial3 = F*A*mR*cR
icl3 = 0

#courant total
id3 = id(E,ial3,icl3)
idv3=idv(ial3,icl3)
Edv3=Edv(id3,E,idv3)


###Tracé du graphe

fig, (ax1, ax2,ax3) = plt.subplots(1, 3)

ax1.grid()
#ax1.xlabel("Potentiel E (V)")
#ax1.ylabel("Intensité du courant i (A)")
ax1.set_xlim(0,1.5)
ax1.set_ylim(-1.5,1)
ax1.set_title('Les deux partenaires sont présents.')

ax1.axhline(y=0, xmin=-1, xmax=2, color="black") #axe des abscisses
ax1.scatter(Edv1,0) # Ajout d'un point
ax1.annotate("Potentiel de demi-vague",(Edv1,-0.1)) # Ajout d'une annotation
ax1.plot(E,id1)


ax2.grid()
#ax2.xlabel("Potentiel E (V)")
#ax2.ylabel("Intensité du courant i (A)")
ax2.set_xlim(0,1.5)
ax2.set_ylim(-1.5,1)
ax2.set_title('Les deux partenaires sont présents.')

ax2.axhline(y=0, xmin=-1, xmax=2, color="black") #axe des abscisses
ax2.scatter(Edv2,idv2) # Ajout d'un point
ym=idv2/2.5
ax2.axvline(Edv2,ymin=.6,ymax=.6+ym, linestyle=":") #axe des abscisses
ax2.annotate("Potentiel de demi-vague",(.75,-0.1)) # Ajout d'une annotation
ax2.plot(E,id2)


ax3.grid()
#ax2.xlabel("Potentiel E (V)")
#ax2.ylabel("Intensité du courant i (A)")
ax3.set_xlim(0,1.5)
ax3.set_ylim(-1.5,1)
ax3.set_title('Seul le réducteur est présent.')

ax3.axhline(y=0, xmin=-1, xmax=2, color="black") #axe des abscisses
ax3.scatter(Edv3,idv3) # Ajout d'un point
ym=idv3/2.5
ax3.axvline(Edv3,ymin=.6,ymax=.6+ym, linestyle=":") #axe des abscisses
ax3.annotate("Potentiel de demi-vague",(Edv3,-0.1)) # Ajout d'une annotation
ax3.plot(E,id3)


plt.show()
