from curses import ncurses_version
from multiprocessing.connection import wait
import numpy as np 
from numpy import linalg as lg
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import math as mt

#BEGIN  Prueba Eq. points

#Parmetros
N = 10000
beta_R = 0.10 #Beta_x
delta_R = 0.5 #Delta_x
Nu = 0.25 #Nu
alpha = 0.5 #Alpha
beta_C = 0.2 #Beta_y
delta_C = 0.15 #Delta_y

#Parametros normalizados
tau = 1/delta_R
b_R = beta_R*tau
d_R = delta_R*tau
u = Nu*tau
a = alpha*tau
b_C = beta_C*tau
d_C = delta_C*tau

A_i = 0.1
A_f = 1
B_i = 0.1
B_f = 2

div = 100

b_R = np.linspace(A_i, A_f, div)*tau
b_C = np.linspace(B_i, B_f, div)*tau

#Puntos de equilibrio 1
x_1 = 0
y_1 = 0
z_1 = 0

#Puntos de equilibrio 2
x_2 = 1 - 1/b_R
y_2 = 0
z_2 = 0


# Puntos de equilibrio 3
#x_3 = d_C/a*(u + d_C)/(b_C - d_C)
#y_3 = -b_R/a*((1- b_R)/b_R - d_C/a*(u + d_C)/(b_C - d_C))
#z_3 = -b_R/a*d_C/(b_C - d_C)*((1- b_R)/b_R - d_C/a*(u + d_C)/(b_C - d_C))

#Jacobiano
def J(x,y,z,i,j):
    return [[-d_R + b_R[i] - 2*b_R[i]*x - a*y , -a*x , 0],
            [-a*y , -d_C - a*x , u + b_C[j]],
            [0 , a*x , -d_C - u]]

eigvals = np.zeros((div**2,3,3),dtype = 'complex_')

n = 0

for i in range(len(b_R)):

    for j in range(len(b_C)):

        x_3 = d_C/a*(u + d_C)/(b_C[j] - d_C)
        y_3 = -b_R[i]/a*((1- b_R[i])/b_R[i] - d_C/a*(u + d_C)/(b_C[j] - d_C))
        z_3 = -b_R[i]/a*d_C/(b_C[j] - d_C)*((1- b_R[i])/b_R[i] - d_C/a*(u + d_C)/(b_C[j] - d_C))

        J_1 = J(x_1,y_1,z_1,i,j)
        l1 = np.linalg.eigvals(J_1) #lg.eigenvals()
        eigvals[n][0][0] = l1[0]
        eigvals[n][1][0] = l1[1]
        eigvals[n][2][0] = l1[2]

        if (b_R[i] > 1):
            J_2 = J(x_2[i],y_2,z_2,i,j)
            l2 = np.linalg.eigvals(J_2) #lg.eigenvals()
            eigvals[n][0][1] = l2[0]
            eigvals[n][1][1] = l2[1]
            eigvals[n][2][1] = l2[2]
        
        else:
            eigvals[n][0][1] = np.nan
            eigvals[n][1][1] = np.nan
            eigvals[n][2][1] = np.nan

        if (d_C/a*(u + d_C)/(b_C[j] - d_C) > (1 - b_R[i])/b_R[i]):
            J_3 = J(x_3,y_3,z_3,i,j)
            l3 = np.linalg.eigvals(J_3) #lg.eigenvals()
            eigvals[n][0][2] = l3[0]
            eigvals[n][1][2] = l3[1]
            eigvals[n][2][2] = l3[2]
        else:
            eigvals[n][0][2] = 0.
            eigvals[n][1][2] = 0.
            eigvals[n][2][2] = 0.

        n = n + 1
 

black = (np.isnan(eigvals)).astype(int)
green = 2*(eigvals == 0.).astype(int)
p = 3*((np.isreal(eigvals)) & (eigvals > 0)).astype(int)
n = 4*((np.isreal(eigvals)) & (eigvals < 0)).astype(int)
pc = 5*((np.iscomplex(eigvals)) & (eigvals > 0)).astype(int)
nc = 6*((np.iscomplex(eigvals)) & (eigvals < 0)).astype(int)

teigvals = black + green + p + n + pc + nc

#np.savetxt("eig1.txt", np.array(list(itertools.zip_longest(*[teigvals[:,0,0], teigvals[:,1,0], teigvals[:,2,0]]))))
#np.savetxt("eig2.txt", np.array(list(itertools.zip_longest(*[teigvals[:,0,1], teigvals[:,1,1], teigvals[:,2,1]]))))
#np.savetxt("eig3.txt", np.array(list(itertools.zip_longest(*[teigvals[:,0,2], teigvals[:,1,2], teigvals[:,2,2]]))))

#END espai parÃ metres

#BEGIN mapeo
if (teigvals > 4):
    print ('YAS')
else:
    print('NO')
    
#END mapeo

