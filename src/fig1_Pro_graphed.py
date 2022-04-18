#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 17:13:23 2022

name:   All_modeled.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_1_nutrient_additions/src

@author: dkm

goal: import and all treatments of HOOH on Pro iin unammended sea water 

working on: importing and visualising data from Morris et al 2011 fig 3a

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import *


############################

#  Data Import from csv   

############################



df_all = pd.read_csv("../data/avgs.csv")
df_all['avg_exp'] = df_all['avg_exp'].fillna(value = 0.0) #filling Nans with 0.0 in 'avg' column 
df_all = df_all.rename({'Time(days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'

####################################

# Slicing Data

####################################


treatments = [50,200,400,800,10000]  #HOOH nM []   #1000 or 10000 nM HOOH for last treatment....kdam no longer saturation if later case
treatments_smol = treatments[1:3]
dc = dict()

####slicing df into treatments and saving into dictionary (dc)######
for i in treatments:
    df_i = df_all[df_all["treatment"].isin([i])]
    name = ('df_' + str(i))
    dc.update({name : df_i})  #update dictionary of dfs with each loop itteration. 
    #print(dc)

########################

#graphing

#########################
plt.rcParams['font.size'] = 14

fig1,ax1 = plt.subplots()

colors = ('green', 'c', 'orange', 'r', 'k') #make into a set in loop? for c in colors, color = count(c)?????
markers = ('s','v','o','*','d')
#for df_i in dc:
for t in treatments: 
    count = treatments.index(t)
    #print(count)
    df = dc['df_'+str(t)]
    times = df['times']
    data = df['avg_exp'] #data was loggeed in original graph; transformed in excel before read in
    ax1.plot(times, data, linestyle = 'None', marker= markers[count], markersize= 10, label = (str(t) +' nM HOOH'), color = colors[count])  #color = colors(i))
    ax1.plot(times,data,linestyle='-', linewidth=0.25, color='black', marker = 'None')
#plt.show()
ax1.set(xlabel= 'Time (days)', ylabel='Biomass ( cells  ml$^{-1}$)')
ax1.set_title = ('Prochlorococcus U18301 w/ HOOH') #not showing up for some reason? 
ax1.legend(loc = 'lower left', prop={"size":10})
ax1.semilogy()
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)


############################

#model

#############################

'''
#dPdt = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
#dSdt = (Supply of nutriet) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)
#Qn = (9.4e-15*(1/14.0)*1e+9)  #Nitrogen Quota for Pro 
'''
step = 0.01
ndays = 2.1
times = np.linspace(0,ndays,int(ndays/step))
Qn = (9.4e-15*(1/(14.0))*1e+9)   #Nitrogen Quota for Pro from Bertillison? 

#S_base =   30.0             #3.0 ~160nM from BCC paper calculations 
     
P = (5e5) # units are cells per mL
#S = (S_base + 0)    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022

#k1= 0.002   #ALPHA  
#k2 = 0.22  #VMAX
#nrdelta = 0.06      #nutrient replete delta  #Kdam for this 
#nddelta = 0.19       #nutrient deplete delta kddam for this
#mumax = k2    #set max/min for this known from lit? 

#parameter values (k1 = alpha, k2 = Vmax, kdam = initial HOOH damage, kddam = HOOH damage when N runs out)
k1s = np.r_[[0.02, 0.02, 0.02, 0.02, 0.02]]*1e+100     
k2s = [0.32, 0.32, 0.32, 0.32, 0.32]
kdams = [0.07,0.2, 0.6, 1.3, 3.1]
kddams = [0.03, 0.5, 2, 2, 1.4]

params = list(zip(k1s,k2s,kdams,kddams))

fig2,ax2 = plt.subplots()
#df = dc['df_50']
#data = df['avg_exp']
#S_base =  (data.max()- data.iloc[0])*Qn 

for t in treatments: 
    count = treatments.index(t)
    k1 = params[count][0]
    k2 = params[count][1]
    kdam = params[count][2]
    kddam = params[count][3]
    ks = (k2/k1)   #set max value for this that is know from lit? (between 0.01 and 0.015 for N metabolism in )
    SsEuler = np.array([])
    PsEuler = np.array([])
    P = 1.3e5
    #print(kdam,kddam)
    S_base =   162.0             #3.0 ~160nM from BCC paper calculations 
    df = dc['df_50']
    data = df['avg_exp']
    S =  (data.max()- data.iloc[0])*Qn    #using QN and 50 treatment as base N for test.     
    # units are cells per mL
    #S = S_base    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022
    for t in times:
            PsEuler = np.append(PsEuler,P)
            SsEuler = np.append(SsEuler,S)
            #if (S>2e-3):
            #    delta = kdam
            #else:
            #    delta = kddam
            delta = kdam
            #print(S)
            #print(P)
            dPdt = k2 * P * S /((ks) + S) - delta*P
            dSdt = -P*(k2*(S)/((ks)+S))*Qn
            if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                    S = S + dSdt*step  #S = 4e-47
            else:
                    S = S + dSdt*step
            P = P + dPdt*step
    ax1.plot(times,(PsEuler), linestyle = 'dashed', color = colors[count]) 
    ax2.plot(times,(SsEuler), linestyle = 'dashed', color = colors[count])


    
#ax2_legend = [zip(treatments,colors)]  #ipz doesn't give a printable object for legend to show
ax2.set(xlabel= 'Time (days)', ylabel='Nitrogen (nM)')
ax2.set_title = ('Nutrient dynamics') #not showing up for some reason? 
#ax2.legend([(list(a2_legend))] ,loc = 'lower left')
ax2.semilogy()

fig3,ax3 = plt.subplots()

for y, c, m, in zip(kdams, colors, markers):
    ax3.plot(kdams, treatments,linestyle = 'none', color = c, marker = m)
#ax3.plot(kdams,treatments,linestyle = 'none', marker=markers[kdam[]],color=colors[count])
ax3.set(xlabel= 'HOOH treatment (nM)', ylabel='kdam (day$^{-1}$)')
ax3.semilogy()
plt.show()



fig1.savefig('../figures/biomass')
fig2.savefig('../figures/nutrients')
