# -*- coding: utf-8 -*-
"""
Created on Tue May 10 12:51:20 2022

@author: hariharan
"""

import json

import pandas as pd
import numpy as np
import os.path
import sys
import TongBeckerman
import matplotlib.pyplot as plt

#Read json file for reading parameters

with open(sys.argv[1], 'r') as myfile:
    data=myfile.read()

# parse file
obj = json.loads(data)


co=np.array(obj['Alloy_composition'])
ko=np.array(obj['Eqm_partition_coefficient'])
mo=np.array(obj['Eqm_liquidus_slope'])
tl=obj['Liquidus_temperature']
delta_T=obj['Solidification_range']
solutes=obj["Solutes"]
R=obj["Tip_radius"]
G=obj["Thermal_gradient"]
Tdot=obj["Cooling_rate"]
Ttip=obj['Tip_temperature']
cl_tip=np.array(obj['Tip_composition'])
kv_tip=np.array(obj['kinetic_partition_coeff'])
D=np.array(obj["Solid_diffusivity"])
DL=np.array(obj["Liquid_diffusivity"])
fname=obj['File_name']

l=obj["PDAS"]
#l=np.sqrt(3*delta_T*R/G)
t=delta_T/Tdot
Pe=np.divide((R*Tdot/G),(2*DL))
file_exists = os.path.exists(fname)

if file_exists:
    nname=fname[:-4]+'_1.csv'
    print("{} File already present-Creating another file {}".format(fname, nname))
    fname=nname
    
    

calc= TongBeckerman.TongBeckerman(co,ko,D,cl_tip,kv_tip,Pe,l,t)

dataset=calc.solve()

fs = np.arange(0.0001,1.0, 0.0001)

T=[tl+np.sum(mo*(i-co)) for i in dataset]
cl=['clstar_'+i for i in solutes]
df=pd.DataFrame(np.row_stack(dataset),columns=cl)
df['fs']=fs
df['T']=T
df.to_csv(fname,index=True)
plt.plot(df["fs"],df["clstar_Mo"])
#plt.ylim(1250,1650)
plt.show()