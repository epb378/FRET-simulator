#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 14:28:35 2020

@author: littleneddyb
"""
from mpl_toolkits.mplot3d import Axes3D #(3D scatter https://matplotlib.org/3.1.1/gallery/mplot3d/scatter3d.html)
import numpy as np
import matplotlib.pyplot as plt
import math
import random
#set up grid of Xmax*Ymax*Zmax points
Xmax, Ymax, Zmax = 10, 10, 10
x = np.linspace(0,1,Xmax+2)
y = np.linspace(0,1,Ymax+2)
z = np.linspace(0,1,Zmax+2)
X,Y,Z = np.meshgrid(x,y,z)
positions = np.vstack([Z.ravel(),Y.ravel(),X.ravel()])
Ros=np.array([[4,2,1],[0.01,4,2],[0.01,0.01,4]])
PLQEs=np.array([0.2,0.9,0.95])
is_perimeter = (positions[0]==np.min(y)) | (positions[0] ==np.max(y)) | (positions[1]==np.min(x)) | (positions[1] ==np.max(x)) | (positions[2]==np.min(z)) | (positions[2] ==np.max(z))

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(*positions[::1], c=is_perimeter)
plt.show()
#populate grid with Nb blue molecules, Ng greens, and Nr reds
Nb = 20
blues=np.ones([Nb,1])
gtb=1
Ng = Nb*gtb
greens=np.ones([Ng,1])
greens=greens*2
rtg=0.1
Nr = math.floor(Ng*rtg)
reds=np.ones([Nr,1])
reds=reds*3 
Ntot=Nb+Ng+Nr
molecules=np.concatenate((reds,greens,blues))
#print(molecules)
np.random.shuffle(molecules)
rest=np.zeros([Ntot,4])
molecules=np.concatenate((molecules,rest), axis=1)
k=0
for i in range (0,Ntot):
    k=k+1
    molecules[i,1]=random.randint(1,Xmax)
    molecules[i,2]=random.randint(1,Ymax)
    molecules[i,3]=random.randint(1,Zmax)
    if i>0:
        for j in range(0,i):
            if molecules[i,1]==molecules[j,1] and molecules[i,2]==molecules[j,2] and molecules[i,3]==molecules[j,3]:
                i = i-1

#checking if my array is actually unique
#check=np.zeros([Ntot,1])
#check=10000*molecules[:,1] + 100*molecules[:,2] +molecules[:,3]
#uniquearrays , uniqueindices = np.unique(check, return_index=True)
distances=np.zeros([Ntot,Ntot])
Pfret=np.zeros([Ntot,Ntot])
Prad=np.zeros([Ntot,Ntot])

for i in range (0,Ntot):
    for j in range(0,Ntot):
        distances[i,j]=np.sqrt((molecules[i,1]-molecules[j,1])**2 + (molecules[i,2]-molecules[j,2])**2 + (molecules[i,3]-molecules[j,3])**2)
        
findblues=np.where(molecules[:,0]==1)
bluesarray=findblues[0]
excitint=random.randint(0,Nb-1)
excitindex=bluesarray[excitint]
print("molecule is number ")
print( excitindex)
molecules[excitindex,4]=1
zeroNtot=np.arange(Ntot)
zeroNtot.shape=(Ntot,1)
for i in range (0,100):
#    print(i)
    somethinghappened=0
    closest=1
    print("molecule is number ")
    print( excitindex)
    while somethinghappened == 0:
        a=distances[:,excitindex] #find distances to closest molecules
        a.shape=(Ntot,1)
        indexeddistances= np.concatenate((zeroNtot,a), axis=1)
        indexeddistances= indexeddistances[indexeddistances[:,1].argsort()] #find index of closest molecule
        if molecules[int(indexeddistances[closest,0]),0] > molecules[excitindex,0]:
            dice=np.random.random()
            Pfret=1/(1+(indexeddistances[closest,1]/Ros[int(molecules[excitindex,0]-1),int(molecules[int(indexeddistances[closest,0]),0]-1)])**6) #fret probability
            Pem=Pfret + (1-Pfret)*PLQEs[int(molecules[excitindex,0]-1)]
            if dice<=Pfret:
                molecules[excitindex,4]=0
                molecules[int(indexeddistances[closest,0]),4]=1
                excitindex=int(indexeddistances[closest,0])
#                print("FRET")
                somethinghappened=1
            elif Pfret< dice <= Pem:
                print("emitted")
                print(molecules[excitindex,0])
                print("light")    
                somethinghappened=1
            else:
                print("non-radiative recombination, bozo")
                somethinghappened=1
        else:
            closest= closest+1
    
        
            
    #find closest molecule lower or equal in energy to excitindex
    #if lower or equal energy roll a dice to see if fret or emit
    #if fret continue
    #if emit stop