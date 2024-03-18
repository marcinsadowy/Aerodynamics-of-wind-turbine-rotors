
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%% plot airfoil
airfoil = pd.read_csv('SG6043.csv',sep='  ')

plt.figure(dpi=1200)
plt.axis('equal')
plt.grid(True)
plt.plot([0,1],[0,0], label="Chord line")
plt.plot(airfoil['x'],airfoil['y'])
plt.plot(0,0,marker='o', label="Leading edge")
plt.plot(1,0,marker='o', label="Trailing edge")
plt.legend()
plt.savefig('airfoil.eps',format='eps')

#%% calculate design parameters
Radius = 0.25
chordLength = 0.03
TSR = 5
Blades = 2
V_0 = 3
V_ang = TSR * V_0 / Radius
alpha = np.deg2rad(9.75)
C_l = 1.38
C_d = 0.056

#Cp target is around 0.4
#Cl/Cd target is 40 (38), which should occur around alpha = 8 deg
#Cl = 1.45, Cd = 0.038


design = pd.DataFrame()
design['r'] = [0.025, 0.05, 0.1, 0.15, 0.2, 0.25]
design['x'] = TSR * design['r'] / Radius
design['a'] = [0.3036, 0.3207,0.3294,0.3315,0.3323,0.3327] #manual input, solved in maple
design['a_prime'] = (1 - 3 * design['a']) / (4 * design['a'] - 1)
design['flowAngle'] = np.arctan((1 - design['a']) / ((1 + design['a_prime']) * design['x']))
design['f'] = (Blades / 2) * ((Radius - design['r']) / (design['r'] * np.sin(design['flowAngle'])))
design['F'] = (2 / np.pi) * np.arccos(np.exp(-1 * design['f']))
design['V_rel'] = np.sqrt(((1 - design['a']) ** 2) * (V_0 ** 2) + ((1 + design['a_prime']) ** 2) * (V_ang ** 2) * (design['r'] ** 2))
design['twistAngle'] = design['flowAngle'] - alpha
design['c'] = (8 * np.pi * design['r'] * design['a'] * design['F'] * (np.sin(design['flowAngle']) ** 2)) / (Blades * (1 - design['a']) * (C_l * np.cos(design['flowAngle']) + C_d * np.sin(design['flowAngle'])))

#print(np.rad2deg(design['twistAngle']))
print(design['c'])
#print((2 / 2) * ((0.25 - design['r']) / (design['r'] * np.sin(design['flowAngle']))))