import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import optimize

"""
Tracé des courbes intensité-potentiel pour la pile Daniell
"""

#%% Définition des constantes

i0 = 0.01
i0_eau = 0.00001
alpha = 0.5
n = 2 #nombre d'électrons échangés
F = 96500
R = 8.314
T = 298
A = 1e-4
DO=6e-6
DR=6e-6
delta=1e-6
cR=0.01
cO=0.01
E_Cu=0.48
E_Zn=0.2
eta_eau=-0.90 #surtension du couple H+/H2 sur une électrode de zinc [en V] On considère que c'est la même sur le cuivre.

f = F/(R*T)
E_cat=E_Cu+R*T/(n*F)*np.log(cO)
E_an=E_Zn+R*T/(n*F)*np.log(cR)
mO=DO/delta
mR=DR/delta


#%% Définition de l'abscisse et des ordonnées

E1 = np.linspace(eta_eau,E_an,10)
E2 = np.linspace(E_an,E_cat,100)
E3 = np.linspace(E_cat,0.75,100)
E=np.concatenate((E0,E1,E2,E3))

N=len(E)
abscisses=np.zeros(N)


#courants limites de diffusion
ial = F*A*mR*cR
icl = -F*A*mO*cO

#Pour les courants anodique/cathodique/totaux en transfert de charge limitant

def i_charge(i0,E,Es):
    #Contrôle de charges limitant
    return i0*np.exp(alpha*n*f*(E-Es))-i0*np.exp(-(1-alpha)*n*f*(E-Es))

def i_diff(i0,E,Es):
    #Pour la diffusion limitante
    return (ial*np.exp(n*f*(E-Es)) + icl)/(1+np.exp(n*f*(E-Es)))

def it(i0,E,Es):
    #Pour un contrôle mixte
    np.seterr(divide='ignore', invalid='ignore')
    return (i_charge(i0,E,Es)*i_diff(i0,E,Es))/(i_charge(i0,E,Es)+i_diff(i0,E,Es))

#Pour l'anode :
i_an1=0*E1
i_an2=i_charge(i0,E2,E_an)

#Pour la cathode :
i_cat1=i_diff(i0,E1,E_cat)
i_cat2=it(i0,E2,E_cat)
i_cat2[99]=0
i_cat3=i_charge(i0,E3,E_cat)

#Définition du potentiel mixte :

def EM(courant1,courant2,potentiel):#recherche du potentiel mixte
    courant=np.abs(courant1+courant2)
    index=courant.argmin()
    return index

Emixte=E2[EM(i_cat2,i_an2,E2)]
imixte=i_an2[EM(i_cat2,i_an2,E2)]

#%% Tracé du graphique

plt.grid()
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(0,0.5)
plt.ylim(-1,1)

#Axe des abscisses
plt.plot(E,abscisses,color='k')

#Pour l'anode :
plt.plot(E1,i_an1,'g--')
plt.plot(E2,i_an2,color="g")

#Pour la cathode :
plt.plot(E1,i_cat1,'b')
plt.plot(E2,i_cat2,'b')
plt.plot(E3,i_cat3,"b")

alpha=imixte/2
plt.scatter(Emixte,0) # Ajout d'un point
plt.axvline(Emixte,ymin=.5-alpha,ymax=.5+alpha, linestyle=":") #pointillés
plt.annotate("Potentiel mixte",(Emixte+0.01,0.05)) # Ajout d'une annotation

plt.show()
