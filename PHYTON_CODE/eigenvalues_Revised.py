from curses import ncurses_version
from multiprocessing.connection import wait
import numpy as np 
from numpy import linalg as lg
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import math as mt
import matplotlib.colors as mcolors

#Parametros fijos normalizados
d_R = 1
a = 10
b_C = 1
d_C = 0.5

div = 1000

b_i = 0.001
b_f = 3
u_i =-3
u_f = 2



b_R = np.linspace(b_i, b_f, div)
u = np.logspace(u_i, u_f, div)
#Puntos de equilibrio 1
# x_1 = 0
# y_1 = 0
# z_1 = 0

#Puntos de equilibrio 2
# x_2 = 1 - 1/b_R
# y_2 = 0
# z_2 = 0


# Puntos de equilibrio 3
#x_3 = d_C/a*(u + d_C)/(b_C - d_C)
#y_3 = -b_R/a*((1- b_R)/b_R - d_C/a*(u + d_C)/(b_C - d_C))
#z_3 = -b_R/a*d_C/(b_C - d_C)*((1- b_R)/b_R - d_C/a*(u + d_C)/(b_C - d_C))

#Jacobiano
def J(x,y,z,i,j):
    return [[-d_R + b_R[i] - 2*b_R[i]*x - a*y , -a*x , 0],
            [-a*y , -d_C - a*x , u[j] + b_C],
            [a*y , a*x , -d_C - u[j]]]

eigvals = np.zeros((len(u),len(b_R)))


for i in range(len(b_R)):

    for j in range(len(u)):

        x_3 = d_C/a*(u[j] + d_C)/(b_C - d_C)
        y_3 = -b_R[i]/a*((1- b_R[i])/b_R[i] - d_C/a*(u[j] + d_C)/(b_C - d_C))
        z_3 = -b_R[i]/a*d_C/(b_C - d_C)*((1- b_R[i])/b_R[i] - d_C/a*(u[j] + d_C)/(b_C - d_C))

        if(b_R[i] < 1):
            eigvals[i][j] = 1 

        if(b_R[i] > 1):

            if(d_C/a*(u[j] + d_C)/(b_C - d_C) > (b_R[i] - 1)/b_R[i]):
                eigvals[i][j] = 2
        

            if(d_C/a*(u[j] + d_C)/(b_C - d_C) < (b_R[i] - 1)/b_R[i]):
                J_3 = J(x_3,y_3,z_3,i,j)
                l3 = np.linalg.eigvals(J_3) #lg.eigenvals()
                #print(l3)
                nr = ((np.isreal(l3[1])) & (l3[1] < 0))
                pr = ((np.isreal(l3[1])) & (l3[1] > 0))
                nc = ((np.iscomplex(l3[1])) & (l3[1] < 0))
                pc = ((np.iscomplex(l3[1])) & (l3[1] > 0))
                if (nr == True):
                    eigvals[i][j] = 3 #Yellow 
                if (pr == True):
                    eigvals[i][j] = 4 #Orange 
                if (nc == True):
                    eigvals[i][j] = 4 #Orange 
                if (pc == True):
                    eigvals[i][j] = 6 #Red 

 

np.savetxt("espai_parametres.txt", eigvals)

cmap1, norm1 = mcolors.from_levels_and_colors([1,2,3,4,5], ['0', 'xkcd:lightish green', 'xkcd:light yellow','xkcd:light orange'])
cmap2, norm2 = mcolors.from_levels_and_colors([5,6,7], ['xkcd:light orange', 'xkcd:cherry red'])

plt.pcolor(eigvals, cmap=cmap1, norm=norm1)
plt.pcolor(eigvals, cmap=cmap2, norm=norm2)
plt.ylim(250,1000)
plt.yticks([333,666,1000], (1,2,3))
plt.xticks([200, 400, 600, 800], (-2, -1, 0, 1))
plt.xlabel('u')
plt.ylabel("b{}".format('\u0052') )
plt.show()
#END espai paràmetres


