print("Last update: 30/06/2021")

import numpy as np
import pandas as pd
import scipy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import argparse

p = argparse.ArgumentParser()
p.add_argument('-i', '--input', type=str, help='Input tensor data file', required=True)
values = p.parse_args()

df = pd.read_table(values.input,skiprows = 1, usecols = [1,3,4,5,6,7,8], delim_whitespace=True, names=['Rg','XX','YY','ZZ','XY','XZ','YZ'])

rx = np.array([])
ry = np.array([])
rz = np.array([])

A = np.array([])
B = np.array([])
C = np.array([])
D = np.array([])
E = np.array([])
F = np.array([])
G = np.array([])
H = np.array([])
I = np.array([])

glob = np.array([])
asph = np.array([])
acyl = np.array([])
anis = np.array([])



####Complete list of eigenvalues####
if input("Output detailed analysis of frames?[y/n]   ") =='y':
 with open('tensor_analysis.txt','w') as ten:
  for i in range(len(df.iloc[0:,])):
   rg1 = df['Rg'].iloc[i,]
   a = df['XX'].iloc[i,]
   b = df['YY'].iloc[i,]
   c = df['ZZ'].iloc[i,]
   d = df['XY'].iloc[i,]
   e = df['XZ'].iloc[i,]
   f = df['YZ'].iloc[i,]
   S = np.array([[a,d,e],
                 [d,b,f],
                 [e,f,c]])
   L = S @ S.T
   eivals, eivecs = la.eig(L)
   eivals = eivals.real
 
   A1, B1, C1 = eivecs[:, 0]
   D1, E1, F1 = eivecs[:, 1]
   G1, H1, I1 = eivecs[:, 2]
 
   A = np.append(A, A1)
   B = np.append(B, B1)
   C = np.append(C, C1)
   D = np.append(D, D1)
   E = np.append(E, E1)
   F = np.append(F, F1)
   G = np.append(G, G1)
   H = np.append(H, H1)
   I = np.append(I, I1)
 
   sort_perm = eivals.argsort()
   eivals.sort()
   rx1, ry1, rz1 = np.sqrt(eivals)

   rx = np.append(rx, rx1)
   ry = np.append(ry, ry1)
   rz = np.append(rz, rz1)
 
 

####Descriptors####
   glob1 = rx1/rz1
   asph1 = rz1 - 0.5*(rx1 + ry1)
   acyl1 = ry1-rx1
   anis1 = (((rz1 - 0.5*(rx1 + ry1))**2)+(3/4) * ((ry1-rx1)**2))/(rg1)**4
#  anis1 = ((0.5**0.5) * ((rx1-ry1)**2 + (ry1-rz1)**2 + (rz1-rx1)**2)**0.5)/(rx1**2 + ry1**2 + rz1**2)**0.5

   glob = np.append(glob, glob1)
   asph = np.append(asph, asph1)
   acyl = np.append(acyl, acyl1)
   anis = np.append(anis, anis1)
   print("Diagonalised eigenvalues for frame ",i+1," : " ,rx1,"   ", ry1,"   ", rz1,file=ten)
   print("Globularity: ", glob1,file=ten)
   print("Asphericity: ", asph1,file=ten)
   print("Acylindricity: ", acyl1,file=ten)
   print("Anisotropy (k^2): ", anis1,file=ten)

  avgglob = np.mean(glob)
  avgasph = np.mean(asph)
  avgacyl = np.mean(acyl)
  avganis = np.mean(anis)

  stdglob = np.std(glob)
  stdasph = np.std(asph)
  stdacyl = np.std(acyl)
  stdanis = np.std(anis)

  maxglob = np.max(glob)
  maxasph = np.max(asph)
  maxacyl = np.max(acyl)
  maxanis = np.max(anis)

  minglob = np.min(glob)
  minasph = np.min(asph)
  minacyl = np.min(acyl)
  minanis = np.min(anis)

####Average radii corresponding to eigenvals####
  avgrx = np.mean(rx)
  avgry = np.mean(ry)
  avgrz = np.mean(rz)

####Descriptive Statistics of eigenvals####
  stdvrx = np.std(rx)
  stdvry = np.std(ry)
  stdvrz = np.std(rz)

  maxrx = np.max(rx)
  maxry = np.max(ry)
  maxrz = np.max(rz)

  minrx = np.min(rx)
  minry = np.min(ry)
  minrz = np.min(rz)
  
  avgrg = np.mean(df['Rg'])
  avgx = np.mean(df['XX'])
  avgy = np.mean(df['YY'])
  avgz = np.mean(df['ZZ'])

  stdrg = np.std(df['Rg'])
  stdx = np.std(df['XX'])
  stdy = np.std(df['YY'])
  stdz = np.std(df['ZZ'])

  maxrg = np.max(df['Rg'])
  maxx = np.max(df['XX'])
  maxy = np.max(df['YY'])
  maxz = np.max(df['ZZ'])

  minrg = np.min(df['Rg'])
  minx = np.min(df['XX'])
  miny = np.min(df['YY'])
  minz = np.min(df['ZZ'])
  print("#############################################################################################################",file=ten)
  print("\t Mean\t\t\t SD\t\t\t Max\t\t Min",file=ten)
  print("Rg:\t",avgrg,"\t",stdrg,"\t",maxrg,"\t",minrg,file=ten)
  print("#############################################################################################################",file=ten)
  print("\t Mean\t\t\t SD\t\t\t Max\t\t\t Min",file=ten)
  print("λx:\t",avgrx,"\t",stdvrx,"\t",maxrx,"\t",minrx,file=ten)
  print("λy:\t",avgry,"\t",stdvry,"\t",maxry,"\t",minry,file=ten)
  print("λz:\t",avgrz,"\t",stdvrz,"\t",maxrz,"\t",minrz,file=ten)
  ####Descriptors####
  print("#############################################################################################################",file=ten)
  print("\t\t Mean\t\t\t SD\t\t\t Max\t\t\t Min",file=ten)
  print("Globularity:\t",avgglob,"\t",stdglob,"\t",maxglob,"\t",minglob,file=ten)
  print("Asphericity:\t",avgasph,"\t",stdasph,"\t",maxasph,"\t",minasph,file=ten)
  print("Acylindricity:\t",avgacyl,"\t",stdacyl,"\t",maxacyl,"\t",minacyl,file=ten)
  print("Anisotropy:\t",avganis,"\t",stdanis,"\t",maxanis,"\t",minanis,file=ten)
  print("#############################################################################################################",file=ten)


else:
 for i in range(len(df.iloc[0:,])):
  rg1 = df['Rg'].iloc[i,]
  a = df['XX'].iloc[i,]
  b = df['YY'].iloc[i,]
  c = df['ZZ'].iloc[i,]
  d = df['XY'].iloc[i,]
  e = df['XZ'].iloc[i,]
  f = df['YZ'].iloc[i,]
  S = np.array([[a,d,e],
                [d,b,f],
                [e,f,c]])
  L = S @ S.T
  eivals, eivecs = la.eig(L)
  eivals = eivals.real

  A1, B1, C1 = eivecs[:, 0]
  D1, E1, F1 = eivecs[:, 1]
  G1, H1, I1 = eivecs[:, 2]

  A = np.append(A, A1)
  B = np.append(B, B1)
  C = np.append(C, C1)
  D = np.append(D, D1)
  E = np.append(E, E1)
  F = np.append(F, F1)
  G = np.append(G, G1)
  H = np.append(H, H1)
  I = np.append(I, I1)

  sort_perm = eivals.argsort()
  eivals.sort()
  rx1, ry1, rz1 = np.sqrt(eivals)

  rx = np.append(rx, rx1)
  ry = np.append(ry, ry1)
  rz = np.append(rz, rz1)



####Descriptors####
  glob1 = rx1/rz1
  asph1 = rz1 - 0.5*(rx1 + ry1)
  acyl1 = ry1-rx1
  anis1 = (((rz1 - 0.5*(rx1 + ry1))**2)+(3/4) * ((ry1-rx1)**2))/(rg1)**4
#  anis1 = ((0.5**0.5) * ((rx1-ry1)**2 + (ry1-rz1)**2 + (rz1-rx1)**2)**0.5)/(rx1**2 + ry1**2 + rz1**2)**0.5

  glob = np.append(glob, glob1)
  asph = np.append(asph, asph1)
  acyl = np.append(acyl, acyl1)
  anis = np.append(anis, anis1)


####Descriptors' Stats####
avgglob = np.mean(glob)
avgasph = np.mean(asph)
avgacyl = np.mean(acyl)
avganis = np.mean(anis)

stdglob = np.std(glob)
stdasph = np.std(asph)
stdacyl = np.std(acyl)
stdanis = np.std(anis)

maxglob = np.max(glob)
maxasph = np.max(asph)
maxacyl = np.max(acyl)
maxanis = np.max(anis)

minglob = np.min(glob)
minasph = np.min(asph)
minacyl = np.min(acyl)
minanis = np.min(anis)

####Average radii corresponding to eigenvals####
avgrx = np.mean(rx)
avgry = np.mean(ry)
avgrz = np.mean(rz)

####Descriptive Statistics of eigenvals####
stdvrx = np.std(rx)
stdvry = np.std(ry)
stdvrz = np.std(rz)

maxrx = np.max(rx)
maxry = np.max(ry)
maxrz = np.max(rz)

minrx = np.min(rx)
minry = np.min(ry)
minrz = np.min(rz)

###Construct eigenvector matrix from avg eivecs####
aA = np.mean(A)
aB = np.mean(B)
aC = np.mean(C)
aD = np.mean(D)
aE = np.mean(E)
aF = np.mean(F)
aG = np.mean(G)
aH = np.mean(H)
aI = np.mean(I)

eigenvecs = np.array([[aA,aB,aC],[aD,aE,aF],[aG,aH,aI]])
 
###Plot of average eigenvalues####
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
 
#3rd eigenvector as the cross product between the 1st and 2nd eigenvectors
eigenvecs[:, 2] = np.cross(eigenvecs[:, 0], eigenvecs[:, 1])

#transpose 3x3 rotation matrix
rot_matrix = np.array([eigenvecs[0, :],
                       eigenvecs[1, :],
                       eigenvecs[2, :]])

center = [0,0,0]

#spherical angles
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

#cartesian coordinates corresponding to spherical angles
x = avgrx * np.outer(np.cos(u), np.sin(v))
y = avgry * np.outer(np.sin(u), np.sin(v))
z = avgrz * np.outer(np.ones_like(u), np.cos(v))

for i in range(len(x)):
 for j in range(len(x)):
  [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rot_matrix) + center

ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b')


####Descriptive Statistics####
avgrg = np.mean(df['Rg'])
avgx = np.mean(df['XX'])
avgy = np.mean(df['YY'])
avgz = np.mean(df['ZZ'])

stdrg = np.std(df['Rg'])
stdx = np.std(df['XX'])
stdy = np.std(df['YY'])
stdz = np.std(df['ZZ'])

maxrg = np.max(df['Rg'])
maxx = np.max(df['XX'])
maxy = np.max(df['YY'])
maxz = np.max(df['ZZ'])
 
minrg = np.min(df['Rg'])
minx = np.min(df['XX'])
miny = np.min(df['YY'])
minz = np.min(df['ZZ'])

print("#############################################################################################################")
print("\t Mean\t\t\t SD\t\t\t Max\t\t Min")
print("Rg:\t",avgrg,"\t",stdrg,"\t",maxrg,"\t",minrg)
#print("XX:  ",avgx,"   ",stdx,"   ",maxx,"   ",minx)
#print("YY:  ",avgy,"   ",stdy,"   ",maxy,"   ",miny)
#print("ZZ:  ",avgz,"   ",stdz,"   ",maxz,"   ",minz)
print("#############################################################################################################")
print("\t Mean\t\t\t SD\t\t\t Max\t\t\t Min")
print("λx:\t",avgrx,"\t",stdvrx,"\t",maxrx,"\t",minrx)
print("λy:\t",avgry,"\t",stdvry,"\t",maxry,"\t",minry)
print("λz:\t",avgrz,"\t",stdvrz,"\t",maxrz,"\t",minrz)
####Descriptors####
print("#############################################################################################################")
print("\t\t Mean\t\t\t SD\t\t\t Max\t\t\t Min")
print("Globularity:\t",avgglob,"\t",stdglob,"\t",maxglob,"\t",minglob)
print("Asphericity:\t",avgasph,"\t",stdasph,"\t",maxasph,"\t",minasph)
print("Acylindricity:\t",avgacyl,"\t",stdacyl,"\t",maxacyl,"\t",minacyl)
print("Anisotropy:\t",avganis,"\t",stdanis,"\t",maxanis,"\t",minanis)
print("#############################################################################################################")
if input("Show ellipsoid?[y/n]   ") =='y':
 plt.show()
