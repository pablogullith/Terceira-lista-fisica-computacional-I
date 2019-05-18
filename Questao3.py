# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:55:16 2019

@author: pablo gullith
"""
#Bibliotecas 
from math import factorial
import matplotlib.pyplot as plt
import numpy as np
#Definição usual de Her
def Her(n,x): 
    if (n == 0):
        return 1
    elif (n == 1):
        return 2*x
    else:
        return 2*x*Her(n-1,x) - 2*(n-1)*Her(n-2,x)
#Definição de Sigma
def sigma(n,x):
    fatorial = factorial(n)
    delta = 1/np.sqrt(2**n * fatorial * np.sqrt(np.pi))
    return(delta*np.exp(- x**2/2)*Her(n,x))
#Primeiro Graph
x = np.linspace(-4,4,1000)
sigma0 = sigma(0,x)
sigma1 = sigma(1,x)
sigma2 = sigma(2,x)
sigma3 = sigma(3,x)
plt.plot(x,sigma0,label='sigma_0(x)') 
plt.plot(x,sigma1,label='sigma_1(x)')
plt.plot(x,sigma2,label='sigma_2(x)')  
plt.plot(x,sigma3,label='sigma_3(x)') 
plt.xlabel('x')
plt.ylabel('Sigma_n(x)')
plt.title('Funções de onda')
plt.legend()
plt.show()
plt.style.use('seaborn')
#Segundo Graph
plt.clf()
x = np.linspace(-10,10,1000)
sigma30 = sigma(30,x)
plt.plot(x,sigma30)
plt.xlabel('x')
plt.ylabel('Sigma_30(x)')
plt.title('Função de onda para n=30')
plt.show()
plt.style.use('seaborn')

#Integral Usando o Método da quadratura gaussiana feita pelo Mark Newman

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

def f(J):
    x = J/(1-J**2)
    n = 5
    return x**2 * sigma(n,x)**2 * (1 + J**2)/(1 - J**2)**2

def Quadrature(N,a,b):
    x,w = gaussxwab(N,a,b)
    return np.sum(w*f(x))


a = -1 + 10**(-10)
b = 1 - 10**(-10)
N = 100

R = Quadrature(N,a,b)

print("\nIncerteza Quantica calculada pelo método da quadratura:",np.sqrt(R))

