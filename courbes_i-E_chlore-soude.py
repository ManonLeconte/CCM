import numpy as np
import matplotlib.pyplot as plt


"""
Tracé des courbes intensité-potentiel pour le procédé chlore-soude
Manon Leconte
"""

#Données issues des Techniques de l'ingénieur, J4804 v1


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
T = 70+273 #température de l'électrolyse
A = 1e-4
DO=6e-6
DR=6e-6
delta=1e-6
E_Cl=1.358 #potentiel standard du couple Cl2/Cl- [en V]
eta_Cl=0.010 #surtension du couple Cl2/Cl- sur l'anode [en V]
a_Cl=4.44 #activité des ions chlorures (c(NaCl) = 260 g/L)
a_Cl2=0.7 #activité du dichlore
E_eau=0.0 #potentiel standard du couple H+/H2 [en V]
eta_eau=-0.25 #surtension du couple H2O/H2 sur la cathode [en V]
a_H2=0.7
pH=8 #pH du compartiment cathodique

f = F/(R*T)

E_cat=E_eau+R*T/(n*F)*np.log(10**(-pH)/a_H2)
E_an=E_Cl+R*T/(n*F)*np.log(a_Cl2/a_Cl**2)

mO=DO/delta
mR=DR/delta

#%% Définition de l'abscisse et des ordonnées

E0 = np.linspace(-1.5,E_cat + eta_eau,100)
n0=len(E0)
E1 = np.linspace(E_cat + eta_eau,E_an+eta_Cl,100)
n1=len(E1)
E2 = np.linspace(E_an+eta_Cl,2.5,100)
n2=len(E2)


#courants limites de diffusion pour le couple du chlore
ial = F*A*mR*a_Cl
icl = -F*A*mO*a_Cl2

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
i_an0=np.zeros(n0)
i_an1=np.zeros(n1)
i_an2=it(i0,E2,E_an+eta_Cl)

#Pour la cathode :
i_cat0=i_charge(i0_eau,E0,E_cat + eta_eau)
i_cat1=np.zeros(n1)
i_cat2=np.zeros(n2)


#%% Tracé du graphique

fig = plt.figure()
ax = fig.add_subplot(111)

plt.title("Courbes intensité-potentiel pour le procédé chlore-soude")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(-1.5,2.5)
plt.ylim(-300,300)
#Anode
plt.text(E_an,-30,'$E_{an,th}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.arrow(E_an,-5,dx=0,dy=10,ec='k')
plt.arrow(x=E_an,y=20,dx=.2,dy=0,head_width=10,head_length=.1,fc=jaunefonce,ec=jaunefonce)
plt.text(1.4,40,'$\eta_a$', horizontalalignment = 'center', verticalalignment = 'center',color=jaunefonce)
plt.text(1.65,150,'Cl$^-\longrightarrow$ Cl$_2$', horizontalalignment = 'center', verticalalignment = 'center')
#Cathode
plt.text(E_cat,-30,'$E_{cat,th}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.arrow(E_cat,-5,dx=0,dy=10,ec='k')
plt.arrow(x=E_cat,y=20,dx=-.6,dy=0,head_width=10,head_length=.1,fc=jaunefonce,ec=jaunefonce)
plt.text(-.6,40,'$\eta_c$', horizontalalignment = 'center', verticalalignment = 'center',color=jaunefonce)
plt.text(-1,-150,'H$_2\leftarrow$ H$^{+}$', horizontalalignment = 'center', verticalalignment = 'center')
#Tension à imposer
plt.arrow(x=-.8,y=-60,dx=2.3,dy=0,head_width=10,head_length=.1,fc=jaunefonce,ec=jaunefonce)
plt.arrow(x=1.5,y=-60,dx=-2.35,dy=0,head_width=10,head_length=.1,fc=jaunefonce,ec=jaunefonce)
plt.text(.5,-90,'$E_{cat}-E_{an}$',color=jaunefonce)

#Axe des ordonnées
plt.axvline(x=0,color='k')

#Pour l'anode :
plt.plot(E0,i_an0,color=violet)
plt.plot(E1,i_an1,color=violet)
plt.plot(E2,i_an2,color=violet,label="à l'anode")

#Pour la cathode :
plt.plot(E0,i_cat0,color=bleuvert,label="à la cathode")
plt.plot(E1,i_cat1,color=bleuvert)
plt.plot(E2,i_cat2,color=bleuvert)


plt.legend()
plt.show()
