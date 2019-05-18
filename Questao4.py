# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:51:38 2019

@author: pablo gullith
"""
G = 6.674e-11
L = 10
M = 10000
sigma = M/L**2

import matplotlib.pyplot as plt
import numpy as np

#Guassiana dupla a partir do algoritmo de mark newman

def gaussxw(Ni,Nj):

   
    a = np.linspace(3,4*Ni-1,Ni)/(4*Ni+2)
    b = np.linspace(3,4*Nj-1,Nj)/(4*Nj+2)
    x = np.cos(np.pi*a+1/(8*Ni*Ni*np.tan(a)))
    y = np.cos(np.pi*b+1/(8*Nj*Nj*np.tan(b)))

    tol = 1e-15
    teta = 100.0
    while (teta > tol):
        P0x = np.ones(Ni,float)
        P1x = np.copy(x)
        for k in range(1,Ni):
            P0x,P1x = P1x,((2*k+1)*x*P1x-k*P0x)/(k+1)
        dPx = (Ni+1)*(P0x-x*P1x)/(1-x**2)
        dx = P1x/dPx
        x = x - dx
        teta = np.max(abs(dx))
    teta = 100.0
    while (teta > tol):
        P0y = np.ones(Nj,float)
        P1y = np.copy(y)
        for k in range(1,Nj):
            P0y,P1y = P1y,((2*k+1)*y*P1y-k*P0y)/(k+1)
        dPy = (Nj+1)*(P0y-y*P1y)/(1-y**2)
        dy = P1y/dPy
        y = y - dy
        teta = np.max(abs(dy))   

    wx = 2*(Ni+1)*(Ni+1)/(Ni*Ni*(1-x*x)*dPx**2)
    wy = 2*(Nj+1)*(Nj+1)/(Nj*Nj*(1-y*y)*dPy**2)

    return x,y,wx,wy

#Definições usuais

def gaussxw_a(Ni,Nj,xi,xf,yi,yf):
    x,y,wx,wy = gaussxw(Ni,Nj)
    return 0.5*(xf-xi)*x+0.5*(xf+xi),0.5*(xf-xi)*wx , 0.5*(yf-yi)*y+0.5*(yf+yi),0.5*(yf-yi)*wy

def f(x,y,z):
    return (x**2+y**2+z**2)**(-3/2)

def Gaussian(Ni,Nj,xi,xf,yi,yf,z):
    x,wx,y,wy = gaussxw_a(Ni,Nj,xi,xf,yi,yf);
    s = 0
    for i in range(Ni):
        for j in range(Nj):
            s = s + wx[i]*wy[j]*f(x[i],y[j],z)   
    return s    

xi = yi = -L/2
xf = yf = L/2

Ni = 100
Nj = 100

z = np.linspace(0,10,1000)
R = Gaussian(Ni,Nj,xi,xf,yi,yf,z)
Fz = G*sigma*z*R

plt.plot(z,Fz)
plt.xlabel('z')
plt.ylabel('Módulo')
plt.title('Força gravitacional')
plt.show()
plt.style.use('seaborn')

Ni = 1500
Nj = 1500

def S(i,N): 
    if (i == 0 or i == N):
        return(1/2)
    else:
        return 1

def Simp2(xi,xf,yi,yf,Ni,Nj,z):
    hx = (xf - xi)/Ni
    hy = (yf - yi)/Nj
    s = 0.0
    for i in range(Ni+1):
        for j in range(Nj+1):
            s = s + S(i,Ni)*S(j,Nj)*f(xi + i*hx, yi + j*hy, z)
    return (hx*hy*s)

R = Simp2(xi,xf,yi,yf,Ni,Nj,z)
Fz = G*sigma*z*R
plt.clf()
plt.plot(z,Fz)
plt.xlabel('z')
plt.ylabel('Módulo')
plt.title('Força gravitacional')
plt.show()
plt.style.use('seaborn')
#Comentários
# Na quadratura gaussiana algumas funções não são bem tratadas
# Para contornar isso, nós podemos usar outro método de integração ou 
# Fazer o aumento de N, assim chegamos a um maior número de pontos perto da
# origem, mas seriam necessários inúmeros deles. Com o maior número de pontos
# Maior a lentidão do programa rodar.
#Usamos 1500 porque com números menores, os grafícos não são satisfatórios