
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%% plot airfoil
airfoil = pd.read_csv('SG6043.csv',sep='  ') / 25
airfoil['x'] = (airfoil['x'] * -1) + 0.04


plt.figure(dpi=1200)
plt.axis('equal')
plt.grid(True)
plt.plot([0,0.04],[0,0], label="Chord line")
plt.plot(airfoil['x'],airfoil['y'])
plt.plot(0.04,0,marker='o', label="Leading edge")
plt.plot(0,0,marker='o', label="Trailing edge")
plt.legend()
plt.savefig('airfoil.eps',format='eps')

#%%Reynold number calculation
rho = 1.2
c = 0.04
mu = 1.8*10**-5
#We use the relative velocity at 2/3 the length of the blade
U = np.sqrt(3**2+(7*3/(0.25)*(0.25*2/3))**2)
Re = rho * U * c / mu

#%% calculate design parameters
Radius = 0.25
chordLength = 0.03
TSR = 7
Blades = 2
V_0 = 3
V_ang = TSR * V_0 / Radius
alpha = np.deg2rad(8)
C_l = 1.45
C_d = 0.038

#Cp target is around 0.4
#Cl/Cd target is 40 (38), which should occur around alpha = 8 deg
#Cl = 1.45, Cd = 0.038


design = pd.DataFrame()
design['r'] = range(0,26)
design['r'] /= 100
design['x'] = TSR * design['r'] / Radius
design['a'] = [0.25, 0.2786203602, 0.2971555786, 0.3087706159, 0.3160606572, 0.3207429983, 0.3238490023, 0.3259804685, 0.3274908755, 0.3285926477, 0.3294171973, 0.3300483083, 0.3309322729, 0.3312478390, 0.3315057875, 0.3317191854, 0.3318976320, 0.3320482969, 0.3321766162, 0.3322867682, 0.3323820051, 0.3324648867, 0.3325374497, 0.3326013298, 0.3326578525, 0.3326578525] #manual input, solved in maple
design['a_prime'] = (1 - 3 * design['a']) / (4 * design['a'] - 1)
design['flowAngle'] = np.arctan((1 - design['a']) / ((1 + design['a_prime']) * design['x']))
design['f'] = (Blades / 2) * ((Radius - design['r']) / (design['r'] * np.sin(design['flowAngle'])))
design['F'] = (2 / np.pi) * np.arccos(np.exp(-1 * design['f']))
design['V_rel'] = np.sqrt(((1 - design['a']) ** 2) * (V_0 ** 2) + ((1 + design['a_prime']) ** 2) * (V_ang ** 2) * (design['r'] ** 2))
design['twistAngle'] = design['flowAngle'] - alpha
design['twistAngle_deg'] = np.rad2deg(design['twistAngle'])
design['c'] = (8 * np.pi * design['r'] * design['a'] * design['F'] * (np.sin(design['flowAngle']) ** 2)) / (Blades * (1 - design['a']) * (C_l * np.cos(design['flowAngle']) + C_d * np.sin(design['flowAngle'])))
design['width'] = np.cos(design['twistAngle'])*design['c']
#print(np.rad2deg(design['twistAngle']))
print(design['c'])
#print((2 / 2) * ((0.25 - design['r']) / (design['r'] * np.sin(design['flowAngle']))))

#%% plot cord length profile
plt.figure(dpi=213, figsize = (10.24, 2))
plt.axis('equal')
plt.xlim(0,0.26)
plt.plot(design['r'],(design['c']*-1 ))
plt.plot([0,0.25],[0,0])
plt.savefig('cord profile.pdf')
#%% plot airfoil profile slices

airfoil_r5 = pd.DataFrame()
airfoil_r5['x'] = (airfoil['x']*np.cos(-1*design['twistAngle'][5]) - airfoil['y'] * np.sin(design['twistAngle'][5])) * design['c'][5]/0.04
airfoil_r5['y'] = (airfoil['x']*np.sin(-1*design['twistAngle'][5]) + airfoil['y'] * np.cos(design['twistAngle'][5])) * design['c'][5]/0.04

airfoil_r10 = pd.DataFrame()
airfoil_r10['x'] = (airfoil['x']*np.cos(-1*design['twistAngle'][10]) - airfoil['y'] * np.sin(design['twistAngle'][10])) * design['c'][10]/0.04
airfoil_r10['y'] = (airfoil['x']*np.sin(-1*design['twistAngle'][10]) + airfoil['y'] * np.cos(design['twistAngle'][10])) * design['c'][10]/0.04

airfoil_r15 = pd.DataFrame()
airfoil_r15['x'] = (airfoil['x']*np.cos(-1*design['twistAngle'][15]) - airfoil['y'] * np.sin(design['twistAngle'][15])) * design['c'][15]/0.04
airfoil_r15['y'] = (airfoil['x']*np.sin(-1*design['twistAngle'][15]) + airfoil['y'] * np.cos(design['twistAngle'][15])) * design['c'][15]/0.04

airfoil_r20 = pd.DataFrame()
airfoil_r20['x'] = (airfoil['x']*np.cos(-1*design['twistAngle'][20]) - airfoil['y'] * np.sin(design['twistAngle'][20])) * design['c'][20]/0.04
airfoil_r20['y'] = (airfoil['x']*np.sin(-1*design['twistAngle'][20]) + airfoil['y'] * np.cos(design['twistAngle'][20])) * design['c'][20]/0.04

plt.figure(dpi=213, figsize = (1.5748, 2))
plt.axis('equal')
plt.xlim(0,0.04)
plt.plot(airfoil_r5['x'],airfoil_r5['y'],label='r5',linewidth=0.5)
plt.plot(airfoil_r10['x'],airfoil_r10['y'],label='r10',linewidth=0.5)
plt.plot(airfoil_r15['x'],airfoil_r15['y'],label='r15',linewidth=0.5)
plt.plot(airfoil_r20['x'],airfoil_r20['y'],label='r20',linewidth=0.5)
plt.grid(True)
plt.legend(fontsize=5)

#%% plot individual slices
plt.figure(dpi=213, figsize = (1.5748, 2))
plt.axis('equal')
plt.xlim(0,0.04)
plt.plot(airfoil_r5['x'],airfoil_r5['y'],label='r5',linewidth=0.5)
plt.grid(True)
plt.legend(fontsize=5)
plt.savefig('airfoil_r5.pdf')

plt.figure(dpi=213, figsize = (1.5748, 2))
plt.axis('equal')
plt.xlim(0,0.04)
plt.plot(airfoil_r10['x'],airfoil_r10['y'],label='r10',linewidth=0.5)
plt.grid(True)
plt.legend(fontsize=5)
plt.savefig('airfoil_r10.pdf')

plt.figure(dpi=213, figsize = (1.5748, 2))
plt.axis('equal')
plt.xlim(0,0.04)
plt.plot(airfoil_r15['x'],airfoil_r15['y'],label='r15',linewidth=0.5)
plt.grid(True)
plt.legend(fontsize=5)
plt.savefig('airfoil_r15.pdf')

plt.figure(dpi=213, figsize = (1.5748, 2))
plt.axis('equal')
plt.xlim(0,0.04)
plt.plot(airfoil_r20['x'],airfoil_r20['y'],label='r20',linewidth=0.5)
plt.grid(True)
plt.legend(fontsize=5)
plt.savefig('airfoil_r20.pdf')

print(design['twistAngle_deg'][5])
print(design['twistAngle_deg'][10])
print(design['twistAngle_deg'][15])
print(design['twistAngle_deg'][20])