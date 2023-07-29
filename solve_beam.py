#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 18:58:16 2023

@author: thomasvandiver
"""

"""
The purpose of this code is two determine the Shear and Moment forces associated
with a two-span continuous beam.

We first create a function for determining the beam, build a DataFrame for the
results and then plot the Data
"""
import numpy as np 
import sympy as sp
import pandas as pd 

def solve_beam(l1, l2, q1, q2):
    l=l1+l2 #total length 
    Mx=sp.symbols('Mx') # create symbol Mx
    
    # Calculate Mx 
    Mx=sp.solveset(Mx*l1/3+q1*l1**3/24+Mx*l2/3+q2*l2**3/24, Mx).args[0]
    
    #solve equilibrium equqtions 
    Va, Vb1, Vb2, Vc=sp.symbols('Va Vb1 Vb2 Vc')
    Va, Vb1=sp.linsolve([Va+Vb1-q1*l1, Vb1*l1+Mx-(q1*l1**2)/2], (Va, Vb1)).args[0]
    Vc, Vb2=sp.linsolve([Vb2+Vc-q2*l2, Vb2*l2+Mx-(q2*l2**2)/2], (Vc, Vb2)).args[0]
    Vb=Vb1+Vb2
    
    x1=np.arange(0,l1+0.1,0.1) # create axis x1
    x2=np.arange(0,l2+0.1,0.1) #create axis x2
    
    beam1=pd.DataFrame({"x":x1})
    beam2=pd.DataFrame({"x":x2})
    
    beam1["M"]=Va*beam1.x-(q1*beam1.x**2)/2 # calculate M and store it
    beam2["M"]=Mx-(q2*beam2.x**2)/2+Vb2*beam2.x # Calculate M and store it
    
    beam1["V"]=Va-q1*beam1.x # calculate V and store it
    beam2["V"]=Vb2-q2*beam2.x # calculate Vand store it
    
    beam2.x=beam2.x+l1 #re-assign x for the second span
    
    beam=pd.concat([beam1, beam2]) # concatenate the two dataframes 
    
    return(beam) # return the result

"""
Now that we have a function for the Moment and Shear Calculations we can 
place the results in an organized DataFrame
"""

header=pd.MultiIndex.from_tuples([("combo 1", "M"), ("combo 1", "V"), ("combo 2", "M"), ("combo 2", "V")])

combos=pd.DataFrame(columns=header)
combos["x"]=solve_beam(4, 5, 3.2, 4.5)["x"]

combos["combo 1"]=solve_beam(4, 5, 3.2, 4.5)
combos["combo 2"]=solve_beam(4, 5, 4.5, 3.2)
combos=-combos.set_index("x")

combos=combos.astype("float")
combos.head()

"""
finally the result is plotted using matplotlib
"""
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(211)
ax.invert_yaxis()

combos.loc[:,pd.IndexSlice[:,"M"]].plot(ax=ax)

ax = plt.subplot(212)
ax.invert_yaxis()
combos.loc[:,pd.IndexSlice[:,"V"]].plot(ax=ax)

