# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 12:43:44 2023

@author: hariharan
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tc_python import *


def set_multiple_conditions(solutes,conc,comp_choice):
    if comp_choice=='mass':
        c="W("
    else:
        c="X("
    temp=["s-c"]
    for i in range(len(solutes)):
        temp.append(c+solutes[i]+")="+str(conc[i]))
    return (' '.join(temp))


def get_liquidus(calc_obj):
    l=(calc_obj
     .remove_condition(ThermodynamicQuantity.temperature())
     .set_phase_to_fixed("LIQUID", 1.0)
     .calculate()
     .get_value_of(ThermodynamicQuantity.temperature()))
    calc_obj.set_phase_to_entered("LIQUID",1.0).set_condition("T",1500).calculate()
    return l


database="TCNI10"
solutes=['Al','Co','Cr','Mo','Ti']
base_element="Ni"
comp=[0.01499,0.09398,0.1996,0.09850,0.02723]
comp_choice="mass" #mass for mass fraction, mole for mole fraction
dtemp=1.0 #slightly accurate near eutectic region
cols=['cl_'+i for i in solutes]
fname="haynes282_truncated.csv"



with TCPython() as start:
    # create and configure a single equilibrium calculation
    condition=set_multiple_conditions(solutes,comp, comp_choice)
    eq_calculation = (
        start
            .set_cache_folder(os.path.basename(__file__) + "_cache")
            .select_database_and_elements(database,[base_element]+ solutes)
            .without_default_phases()
            .select_phase("LIQUID")
            .select_phase("FCC_A1")
            .get_system()
            .with_single_equilibrium_calculation()
            .run_poly_command(condition)
            .set_condition("T",1000)
    )

    #liquidus=get_liquidus(eq_calculation)-273.15
    liquidus=1621.51-273.15
    dummyobj=eq_calculation.set_condition("T", liquidus+273.15).calculate()
    fs=1-dummyobj.get_value_of(MoleFractionOfAPhase("LIQUID"))
    #fs=0
    f=[fs]
    i=1
    comp_new=comp
    composition=[comp_new]
    ph=["LIQUID"]
    temp=[liquidus,liquidus-dtemp]
    while fs<0.99:
        cond=set_multiple_conditions(solutes,composition[i-1], comp_choice)
        loop_result=(eq_calculation.run_poly_command(cond).set_condition("T",temp[i]+273.15)
                     .calculate())
        comp_new=[]
        for s in solutes:
            comp_new.append(
                loop_result
                .get_value_of(CompositionOfPhaseAsWeightFraction("LIQUID", s)))
        temp.append(temp[i]-dtemp)
        fs=((1-loop_result.get_value_of(MoleFractionOfAPhase("LIQUID")))*(1-f[i-1]))+f[i-1]
        f.append(fs)
        composition.append(comp_new)
        pha=loop_result.get_stable_phases()
        StrPha = ' '.join([str(elem) for elem in pha])
        ph.append(StrPha)
        i=i+1
        print(i,fs)
    temp.pop(0)
    df=pd.DataFrame(np.row_stack(composition),columns=cols)
    df['fs']=f
    df['T']=temp
    df.to_csv(fname,index=False)
        
        
        