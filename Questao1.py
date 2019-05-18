# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 08:18:06 2019

@author: pablo gullith
"""
import numpy as np
import matplotlib.pyplot as plt

v = 0.001 #valor convertido de Cm^3 para M^3.
ro = 6.002e28
teta = 428
kb = 1.38e-23
#Gaussiana de mark newman
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 100.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

def f(x):
    return(x**4 * np.exp(x)/(np.exp(x)-1)**2)
    
def Int(N,a,b):
    x,w = gaussxwab(N,a,b)
    return np.sum(w*f(x))

def cv(T):
    b = teta/T
    N = 50
    return(9*v*ro*kb*Int(N,0,b)*b**(-3))

       
T = np.linspace(5,500,1000)
C = [cv(Ti) for Ti in T]
plt.plot(T,C)
plt.title('Gráfico do calor específico')
plt.xlabel('Temperatura (K)')
plt.ylabel('C_v')
plt.show()
plt.style.use('seaborn')
