import numpy as np
import matplotlib.pyplot as plt


"""
Tracé des courbes intensité-potentiel pour le procédé chlore-soude
"""

#Données issues des Techniques de l'ingénieur, J4804 v1

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

plt.grid()
plt.title("Courbes intensité-potentiel pour le procédé chlore-soude")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(-1.5,2.5)
plt.ylim(-300,300)
plt.text(E_an,-30,'$E_{an}$ = 1,31 V', horizontalalignment = 'center', verticalalignment = 'center')
plt.axvline(x = E_an, ymin = 0, ymax  = 0.5,color ='gray',ls='--',lw=0.1)
plt.text(E_cat,-30,'$E_{cat}$ = -0,27 V', horizontalalignment = 'center', verticalalignment = 'center')
plt.axvline(x = E_cat, ymin = 0, ymax  = 0.5,color ='gray',ls='--',lw=0.1)
plt.text(1.65,150,'Cl$^-$ -> Cl$_2$', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(-1,-150,'H$_2\leftarrow$ H$^{+}$', horizontalalignment = 'center', verticalalignment = 'center')

#Pour l'anode :
plt.plot(E0,i_an0,'g')
plt.plot(E1,i_an1,'g')
plt.plot(E2,i_an2,color="g",label="à l'anode")

#Pour la cathode :
plt.plot(E0,i_cat0,'b',label="à la cathode de cuivre")
plt.plot(E1,i_cat1,'b')
plt.plot(E2,i_cat2,'b')

plt.legend()
plt.show()