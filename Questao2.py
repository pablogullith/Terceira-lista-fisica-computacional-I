# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 19:22:34 2019

@author: pablo gullith
"""
# Este é o valor correto da constante de Boltzmann:5.670400(40)e-08
# Temos o valor calculado nesse code, que é: 5.670366816083166e-08
# Aqui usamos o método de simpson com 1000 fatias, 
#Comparamos isso com o valor analítico real e possuimos uma aproximação de 8 
#algarismos significativos.

import numpy as n

def f(x):
    return (x**3/((1 - x)**5*(n.exp(x/(1-x)) - 1)))

def J(i,N): 
    if (i == 0 or i == N):
        return (1/3)
    elif (i % 2 == 0):
        return (2/3)
    else:
        return (4/3)

def Simp(a,b,N):
    h = (b-a)/N
    s = 0
    for i in range(N+1):
        s = s + J(i,N)*f(a + i*h)
    return (h*s)

a = 0 + 1e-15
b = 1 - 1e-15
N = 1000

Integral = Simp(a,b,N)

h = 6.62607004e-34
c = 299792458
Boltzmann = 1.38064852e-23 

Sigma = (2*n.pi*Boltzmann**4/(c**2*h**3))*Integral

print('Integral:',Integral)   
print("Constante de Boltzmann:",Sigma)

