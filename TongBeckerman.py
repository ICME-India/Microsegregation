# -*- coding: utf-8 -*-
"""
Created on Tue May 10 12:51:20 2022

@author: hariharan
"""

import numpy as np
from scipy.integrate import odeint
from scipy.special import expi

class TongBeckerman():
    
    def __init__(self,co,ko,D,cl_tip,kv_tip,Pe,l,t):
        
        self.co=co
        self.ko=ko
        self.D=D
        self.cl_tip=cl_tip
        self.kv_tip=kv_tip
        self.cs_tip=kv_tip*cl_tip
        self.Fourier_solid=D*t/(l/2)**2
        self.Iv=np.array([self.ivantsov(i) for i in Pe])
        #self.beta_prime=np.divide((1-self.Iv),(2*self.Iv))
        self.beta_prime=(self.cs_tip-self.co)/(2*(self.co-self.cl_tip))
    
    def ivantsov(self,x):
    	if(x>20): 
    		z=1
    	z = -x*np.exp(x)*expi(-x)
    	return z
        
    def model(self,Clstar, fs):
        
        if fs<0.0002:
            k=self.kv_tip
        else:
            k=self.ko
        
        factor=np.exp((fs-1)/(2*self.beta_prime*fs))
        F1=2*self.beta_prime*fs*(1-factor)
        F2a=(1+(6*self.Fourier_solid))*(self.co-(Clstar*k))
        F2b1=factor/fs
        F2b2=2*self.beta_prime*(1+(6*self.Fourier_solid))*(1-factor)
        F2=F2a+((Clstar-self.co)*(F2b1-F2b2))
        dCdf = F2/F1
        return dCdf
    
    
    def solve(self):
        
        # Time Span of Interest
        fs = np.arange(0.0001,1.0, 0.0001)  # (t0, tf)
        
        print(self.beta_prime)
        # Solving ODE
        sol = odeint(self.model,self.cl_tip,fs)
        
        return sol