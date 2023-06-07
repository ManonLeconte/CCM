# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import *

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#%% DONNEES

E1=0.77 #potentiel standard du couple Fe3+/Fe2+ vs ESH [en V]

E2=1.51 #potentiel standard du couple MnO4-/Mn2+ vs ESH [en V]

R=8.314 #constante des gaz parfaits [en J/K/mol]

T=298 #température [en K]

F=96500 #constante de Faraday [en C/mol]

n=100 #nombre de valeurs dans les tableaux

#%% FONCTION

def pot1(x):
        return E1+R*T/F*np.log(x/(1-x))
def pot2(x):
        return E2+R*T/(F*5)*np.log(5*(x-1))

X1 = np.linspace(0,0.99,n) #début, fin, nombre de points
X2 = np.linspace(1.01,2,n) #début, fin, nombre de points

E = np.concatenate((pot1(X1),pot2(X2)))
X = np.concatenate((X1,X2))

#%% GRAPHE

plt.xlim(0,2)
plt.title("Courbe de titrage potentiométrique à courant nul")
plt.xlabel("x")
plt.ylabel("Potentiel E [V]")

plt.plot(X,E)

plt.show()

#%% OPTIONNEL

#ax.grid(True) #Grille
#ax.legend(loc=2 ,prop={'size':20}) #Legende
