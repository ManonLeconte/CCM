# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 19:40:31 2024

@author: Manon Leconte
"""

import numpy as np
import matplotlib.pyplot as plt


"""
Choix de l'anode pour le procédé chlore-soude
"""

#%%Couleurs

bleuvert=(0,.5,.5)
violet=(.5,0,.5)
rouge=(.5,0,0)
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
eta_Cl=0.010 #surtension du couple Cl2/Cl- sur l'anode en titane
eta_ClP=0.4 #surtension du couple Cl2/Cl- sur l'anode en platine
a_Cl=4.44 #activité des ions chlorures (c(NaCl) = 260 g/L)
a_Cl2=0.7 #activité du dichlore
E_O=1.23 #potentiel standard du couple H2O/O2 [en V]
eta_O=1 #surtension du couple H2O/O2 sur l'anode en platine
eta_OT=1.5 #surtension du couple H2O/O2 sur l'anode en titane

f = F/(R*T)

E_an=E_Cl+R*T/(n*F)*np.log(a_Cl2/a_Cl**2)

mO=DO/delta
mR=DR/delta

#%% Définition de l'abscisse et des ordonnées

E1 = np.linspace(1,E_an+eta_Cl,100)
n1=len(E1)
E2 = np.linspace(E_an+eta_Cl,E_an+eta_ClP,100)
n2=len(E2)
E3 = np.linspace(E_an+eta_ClP,E_an+eta_O,100)
n3=len(E3)
E4 = np.linspace(E_an+eta_O,E_an+eta_OT,100)
n4=len(E4)
E5 = np.linspace(E_an+eta_OT,3.5,100)
n5=len(E5)
E=np.concatenate((E1,E2,E3,E4,E5))

N=len(E)
abscisses=np.zeros(N)


#courants limites de diffusion pour le couple du chlore
ial = F*A*mR*a_Cl
icl = -F*A*mO*a_Cl2

#Pour les courants anodique/cathodique/totaux en transfert de charge limitant

def i_charge(i0,E,Es):
    #Contrôle de charges limitant
    return i0*np.exp(alpha*n*f*(E-Es))

def i_diff(i0,E,Es):
    #Pour la diffusion limitante
    return (ial*np.exp(n*f*(E-Es)) + icl)/(1+np.exp(n*f*(E-Es)))

def it(i0,E,Es):
    #Pour un contrôle mixte
    np.seterr(divide='ignore', invalid='ignore')
    return (i_charge(i0,E,Es)*i_diff(i0,E,Es))/(i_charge(i0,E,Es)+i_diff(i0,E,Es))


#Pour l'anode :
i_an1=np.zeros(n1)
i_ClT1=it(i0,E2,E_an+eta_Cl)
i_ClT2=it(i0,E3,E_an+eta_Cl)
i_ClT3=it(i0,E4,E_an+eta_Cl)
i_ClT4=it(i0,E5,E_an+eta_Cl)
i_ClP1=it(i0,E2,E_an+eta_ClP)
i_ClP2=it(i0,E3,E_an+eta_ClP)
i_ClP3=it(i0,E4,E_an+eta_ClP)
i_ClP4=it(i0,E5,E_an+eta_ClP)
i_OP1=i_charge(i0_eau,E2,E_an+eta_O)
i_OP2=i_charge(i0_eau,E3,E_an+eta_O)
i_OP3=i_charge(i0_eau,E4,E_an+eta_O)
i_OP4=i_charge(i0_eau,E5,E_an+eta_O)
i_OT1=i_charge(i0_eau,E2,E_an+eta_OT)
i_OT2=i_charge(i0_eau,E3,E_an+eta_OT)
i_OT3=i_charge(i0_eau,E4,E_an+eta_OT)
i_OT4=i_charge(i0_eau,E5,E_an+eta_OT)


#%% Tracé du graphique

fig = plt.figure()
ax = fig.add_subplot(111)

plt.title("Choix de l'anode pour le procédé chlore-soude")
plt.xlabel("Potentiel E (V)")
plt.ylabel("Intensité du courant i (A)")
plt.xlim(1,3.5)
plt.ylim(-300,300)
#Anode
plt.text(E_an,-30,'$E_{4,th}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.text(E_O,30,'$E_{3,th}$', horizontalalignment = 'center', verticalalignment = 'center')
plt.arrow(E_an,-5,dx=0,dy=10,ec='k')
plt.arrow(E_O,-5,dx=0,dy=10,ec='k')
#plt.arrow(x=E_an,y=20,dx=.2,dy=0,head_width=10,head_length=.1,fc=jaunefonce,ec=jaunefonce)
#plt.text(1.4,40,'$\eta_{Cl/Ti}$', horizontalalignment = 'center', verticalalignment = 'center',color=jaunefonce)
#plt.text(1.65,150,'Cl$^-\longrightarrow$ Cl$_2$', horizontalalignment = 'center', verticalalignment = 'center')

#Axe des abscisses
plt.plot(E,abscisses,color='k')

#Pour l'anode :
plt.plot(E2,i_ClT1,color=violet,label='Cl$^-\longrightarrow$ Cl$_2$ sur titane')
plt.plot(E3,i_ClT2,color=violet)
plt.plot(E4,i_ClT3,color=violet)
plt.plot(E5,i_ClT4,color=violet)
plt.plot(E2,i_ClP1,color=violet,ls='--',label='Cl$^-\longrightarrow$ Cl$_2$ sur platine')
plt.plot(E3,i_ClP2,color=violet,ls='--')
plt.plot(E4,i_ClP3,color=violet,ls='--')
plt.plot(E5,i_ClP4,color=violet,ls='--')
plt.plot(E2,i_OT1,color=rouge,label='H$_2$O $\longrightarrow$ O$_2$ sur titane')
plt.plot(E3,i_OT2,color=rouge)
plt.plot(E4,i_OT3,color=rouge)
plt.plot(E5,i_OT4,color=rouge)
plt.plot(E2,i_OP1,color=rouge,ls='--',label='H$_2$O $\longrightarrow$ O$_2$ sur platine')
plt.plot(E3,i_OP2,color=rouge,ls='--')
plt.plot(E4,i_OP3,color=rouge,ls='--')
plt.plot(E5,i_OP4,color=rouge,ls='--')

plt.legend()
plt.show()