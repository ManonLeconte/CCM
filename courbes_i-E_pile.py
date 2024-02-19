import numpy as np
import matplotlib.pyplot as plt

"""
Tracé des courbes intensité-potentiel pour la pile Daniell
Manon Leconte
"""

#%%Couleurs

bleuvert=(0,.5,.5)
violet=(.5,0,.5)
jaunefonce=(.5,.5,0)

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
E_Cu=0.34 #potentiel standard du couple Cu2+/Cu [en V]
E_Zn=-0.76 #potentiel standard du couple Zn2+/Zn [en V]
eta_eau=-0.90 #surtension du couple H+/H2 sur une électrode de zinc [en V] On considère que c'est la même sur le cuivre.

f = F/(R*T)
E_cat=E_Cu+R*T/(n*F)*np.log(cO)
E_an=E_Zn+R*T/(n*F)*np.log(cR)
mO=DO/delta
mR=DR/delta

#%% Définition de l'abscisse et des ordonnées

E0 = np.linspace(-1.5,eta_eau,100)
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
i_an0=i_charge(i0_eau,E0,eta_eau)+it(i0,E0,E_an)
i_an1=it(i0,E1,E_an)
i_an2=i_charge(i0,E2,E_an)

#Pour la cathode :
i_cat0=i_charge(i0_eau,E0,eta_eau)+it(i0,eta_eau,E_cat)
i_cat1=it(i0,E1,E_cat)
i_cat2=it(i0,E2,E_cat)
i_cat3=i_charge(i0,E3,E_cat)


#%% Tracé du graphique

plt.title("Courbes intensité-potentiel pour la pile Daniell")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(-1.5,0.75)
plt.ylim(-0.8,1)
#Anode
plt.text(E_an+.1,-0.1,'$E_{an,th}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.arrow(E_an,-.025,dx=0,dy=.05,ec='k')
plt.text(-0.68,0.5,'Zn $\longrightarrow$ Zn$^{2+}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.arrow(x=E_an,y=.05,dx=.05,dy=0,head_width=.03,head_length=.05,fc=jaunefonce,ec=jaunefonce)
plt.text(-.65,.05,'$\eta_a$', horizontalalignment = 'center', verticalalignment = 'center',color=jaunefonce)
#Cathode
plt.text(E_cat+.1,-0.1,'$E_{cat,th}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.arrow(E_cat,-.025,dx=0,dy=.05,ec='k')
plt.arrow(x=E_cat,y=.05,dx=-.1,dy=0,head_width=.03,head_length=.05,fc=jaunefonce,ec=jaunefonce)
plt.text(E_cat-.2,.05,'$\eta_c$', horizontalalignment = 'center', verticalalignment = 'center',color=jaunefonce)
plt.text(0.2,-0.3,'Cu $\leftarrow$ Cu$^{2+}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(-1.13,-0.7,'H$_2\leftarrow$ H$^{+}$', horizontalalignment = 'center', verticalalignment = 'center')
#Tension délivrée
plt.arrow(x=-.7,y=.25,dx=.75,dy=0,head_width=.03,head_length=.05,fc=jaunefonce,ec=jaunefonce)
plt.arrow(x=.1,y=.25,dx=-.75,dy=0,head_width=.03,head_length=.05,fc=jaunefonce,ec=jaunefonce)
plt.text(-.5,.3,'$E_{cat}-E_{an}$',color=jaunefonce)
#Courant
plt.arrow(-.025,icl,dx=0.05,dy=0,ec='k')
plt.arrow(-.025,-icl,dx=0.05,dy=0,ec='k')
plt.text(.15,-icl,'$i_a=i$', horizontalalignment = 'center', verticalalignment = 'center')


#Axes
plt.plot(E,abscisses,color='k')
plt.axvline(x=0,color='k')

#Pour l'anode :
plt.plot(E0,i_an0,color=violet,ls='--')
plt.plot(E1,i_an1,color=violet,ls='--')
plt.plot(E2,i_an2,color=violet,label="à l'anode de zinc")

#Pour la cathode :
plt.plot(E0,i_cat0,color=bleuvert,label="à la cathode de cuivre")
plt.plot(E1,i_cat1,color=bleuvert)
plt.plot(E2,i_cat2,color=bleuvert)
plt.plot(E3,i_cat3,color=bleuvert,ls='--')

plt.legend()
plt.show()
