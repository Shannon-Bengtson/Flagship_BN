# Load packages
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import itertools
import os
import json
from datetime import datetime
import pysmile
import pysmile_license
import sys
import json
sys.path.append('/src/python_classes')
from colormap import rgb2hex
from collections import ChainMap
import re
from scipy.stats import lognorm
import math

#
def create_BN(shape,scale,mu,sigma,L_means,x0_means,k_means,L_sigmas,x0_sigmas,k_sigmas):
    
    # Define binning of continuous hazard input (tsunami and shaking)
    shaking_continuous_for_frag_bounds = [0.01,1]
    shaking_continuous_for_frag_bins = [0.3333,0.6667,1]
    shaking_continuous_for_LS_bounds = [0.01,1]
    shaking_continuous_for_LS_bins = [0.06,0.1,0.2,0.35,0.45,1]
    tsunami_bounds = [0.15,6]
    tsunami_bins = [0.75,1.25,2.5,6]
    LSN_bounds = [0,100]
    LSN_bins = [10,20,30,100]
    
    shaking_continuous_for_frag_labels = [shaking_continuous_for_frag_bounds[0]]+shaking_continuous_for_frag_bins
    shaking_continuous_for_LS_labels = [shaking_continuous_for_LS_bounds[0]]+shaking_continuous_for_LS_bins
    tsunami_labels = [tsunami_bounds[0]]+tsunami_bins
    
    # Initial network
    net = pysmile.Network()

    # Add nodes to network, plus connections
    for node in ['ShakingAndLSNFragility','ShakingAndLSNFragilityA','ShakingAndLSNFragilityB','ShakingAndLSNFragilityC','ShakingAndLSNFragilityD','ShakingForFrag','LiquefactionModulator','CombinedFragility','Tsunami','TsunamiFragility','LandslideFragility','Landslide','ShakingForLS']:
        net.add_node(pysmile.NodeType.CPT,node)
    
    shakenode = net.add_node(pysmile.NodeType.EQUATION,"ShakingContinuous")
    net.set_node_equation(shakenode,f"ShakingContinuous=Gamma({shape},{scale})")
    
    shakenode = net.add_node(pysmile.NodeType.EQUATION,"ShakingContinuousForFrag")
    net.set_node_equation(shakenode,f"ShakingContinuousForFrag=ShakingContinuous")
    
    shakenode = net.add_node(pysmile.NodeType.EQUATION,"ShakingContinuousForLS")
    net.set_node_equation(shakenode,f"ShakingContinuousForLS=ShakingContinuous")
        
    tsunode = net.add_node(pysmile.NodeType.EQUATION,"TsunamiContinuous")
    net.set_node_equation(tsunode,f"TsunamiContinuous=Normal({mu},{sigma})")
    
    runoutmodels = net.add_node(pysmile.NodeType.EQUATION,"RunoutModels")
    net.set_node_equation(runoutmodels,f"RunoutModels=0.8")
    
    runoutmodifier = net.add_node(pysmile.NodeType.EQUATION,"RunoutModelModifier")
    net.set_node_equation(runoutmodifier,f"RunoutModelModifier=Choose(Landslide,0,Binomial(1,RunoutModels))")
    
    LSNSigmas = net.add_node(pysmile.NodeType.EQUATION,"LSNSigmas")
    net.set_node_equation(LSNSigmas,f"LSNSigmas={L_sigmas}/(1+Exp(-{k_sigmas}*(ShakingContinuous-{x0_sigmas})))")
    
    LSNMeans = net.add_node(pysmile.NodeType.EQUATION,"LSNMeans")
    net.set_node_equation(LSNMeans,f"LSNMeans={L_means}/(1+Exp(-{k_means}*(ShakingContinuous-{x0_means})))")
    
    LSNContinuous = net.add_node(pysmile.NodeType.EQUATION,"LSNContinuous")
    net.set_node_equation(LSNContinuous,f"LSNContinuous=Normal(LSNMeans,LSNSigmas)")
        
    # Add outcomes to each node
    for outcome in ['DS0','DS1','DS2','DS3','DS4']:
        net.add_outcome('CombinedFragility',outcome)
        net.add_outcome('ShakingAndLSNFragility',outcome)
        net.add_outcome('TsunamiFragility',outcome)
        net.add_outcome('ShakingAndLSNFragilityA',outcome)
        net.add_outcome('ShakingAndLSNFragilityB',outcome)
        net.add_outcome('ShakingAndLSNFragilityC',outcome)
        net.add_outcome('ShakingAndLSNFragilityD',outcome)
                 
    for lower,upper in zip(shaking_continuous_for_LS_labels[:-1],shaking_continuous_for_LS_labels[1:]):
        net.add_outcome('ShakingForLS',f"from_{lower}_to_{upper}".replace(".","p"))
        
    for lower,upper in zip(shaking_continuous_for_frag_labels[:-1],shaking_continuous_for_frag_labels[1:]):
        net.add_outcome('ShakingForFrag',f"from_{lower}_to_{upper}".replace(".","p"))
    
    for lower,upper in zip(tsunami_labels[:-1],tsunami_labels[1:]):
        net.add_outcome('Tsunami',f"from_{lower}_to_{upper}".replace(".","p"))
        
    for outcome in ['A','B','C','D']:
        net.add_outcome('LiquefactionModulator',outcome)
        
#     for outcome in ['low','medium','high']:
#         net.add_outcome('ShakingForFrag',outcome)
        
#     for outcome in ['zero','low','medium','high']:
#         net.add_outcome('LSN',outcome)
        
#     for outcome in ['zero','low','medium','high']:
#         net.add_outcome('Tsunami',outcome)
        
    for outcome in ['no','yes']:
        net.add_outcome('Landslide',outcome)
        
    for outcome in ['DS0','DS4']:
        net.add_outcome('LandslideFragility',outcome)
        
    # Add arcs
    net.add_arc('ShakingAndLSNFragility','CombinedFragility')
    net.add_arc('TsunamiFragility','CombinedFragility')
    net.add_arc('LandslideFragility','CombinedFragility')
    net.add_arc('ShakingAndLSNFragilityA','ShakingAndLSNFragility')
    net.add_arc('ShakingAndLSNFragilityB','ShakingAndLSNFragility')
    net.add_arc('ShakingAndLSNFragilityC','ShakingAndLSNFragility')
    net.add_arc('ShakingAndLSNFragilityD','ShakingAndLSNFragility')
    net.add_arc('LiquefactionModulator','ShakingAndLSNFragility')
    net.add_arc('ShakingForFrag','ShakingAndLSNFragilityA')
    net.add_arc('ShakingForFrag','ShakingAndLSNFragilityB')
    net.add_arc('ShakingForFrag','ShakingAndLSNFragilityC')
    net.add_arc('ShakingForFrag','ShakingAndLSNFragilityD')
    net.add_arc('Tsunami','TsunamiFragility')
#     net.add_arc('LSN','LiquefactionModulator')
    
#     net.add_arc('Landslide','RunoutModelModifier')
    net.add_arc('RunoutModelModifier','LandslideFragility')
#     net.add_arc('RunoutModels','RunoutModelModifier')
    
#     net.add_arc('ShakingContinuous','ShakingContinuousForFrag')  
#     net.add_arc('ShakingContinuous','ShakingContinuousForLS')  
    
    net.add_arc('ShakingContinuousForFrag','ShakingForFrag')    
    net.add_arc('ShakingContinuousForLS','ShakingForLS')
    net.add_arc('LSNContinuous','LiquefactionModulator')
    
    net.add_arc('ShakingForLS','Landslide')
    net.add_arc('TsunamiContinuous','Tsunami')

    # Delete the two default states
    for node in net.get_all_nodes():
        if net.get_node_id(node) not in ['ShakingContinuous','TsunamiContinuous','RunoutModelModifier','RunoutModels',
                                        'ShakingContinuousForFrag','ShakingContinuousForLS','LSNMeans','LSNSigmas','LSNContinuous']:
            net.delete_outcome(node,'State0')
            net.delete_outcome(node,'State1')
    
#     # this will need to be changed
#     for node in ['Landslide']:
#         net.set_node_definition(node,[1/net.get_outcome_count(node)]*net.get_outcome_count(node))

#     from IPython import embed
#     embed()

#     return(net)
    
    net.set_node_equation_bounds("ShakingContinuousForLS",shaking_continuous_for_frag_bounds[0],shaking_continuous_for_frag_bounds[1])
    iv = [pysmile.DiscretizationInterval("",x) for x in shaking_continuous_for_frag_bins]
    net.set_node_equation_discretization('ShakingContinuousForFrag',iv)
    
    net.set_node_equation_bounds("ShakingContinuousForLS",shaking_continuous_for_LS_bounds[0],shaking_continuous_for_LS_bounds[1])
    iv = [pysmile.DiscretizationInterval("",x) for x in shaking_continuous_for_LS_bins]
    net.set_node_equation_discretization('ShakingContinuousForLS',iv)
    
    net.set_node_equation_bounds("TsunamiContinuous",tsunami_bounds[0],tsunami_bounds[1])
    iv = [pysmile.DiscretizationInterval("",x) for x in tsunami_bins]
    net.set_node_equation_discretization('TsunamiContinuous',iv)
    
    iv = [pysmile.DiscretizationInterval("",x) for x in [0.5,1]]
    net.set_node_equation_discretization('RunoutModelModifier',iv)
    
    net.set_node_equation_bounds("LSNContinuous",LSN_bounds[0],LSN_bounds[1])
    iv = [pysmile.DiscretizationInterval("",x) for x in LSN_bins]
    net.set_node_equation_discretization('LSNContinuous',iv)
    
    
    net.set_node_definition('ShakingForLS',[
        1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,1
                                           ])
    
#     from IPython import embed
#     embed()
    
    net.set_node_definition('ShakingForFrag',[
        1,0,0,
        0,1,0,
        0,0,1,
                                           ])
     
        
    net.set_node_definition('Tsunami',[
        1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1
                                           ])
    
    net.set_node_definition('LiquefactionModulator',[
        1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1
                                           ])

        
    return(net)

def create_combining_CPT():

    # Setup dataframe defining how to combined individual fragilities, taking the maximum damage state
    damage_states_array = np.array((0,1,2,3,4))
    landslide_states_array = np.array((0,4))

    fragility_array = np.array(('ShakingForFrag','Tsunami','Landslide'))
    rows = list(itertools.product(list(itertools.product(damage_states_array,damage_states_array)),landslide_states_array))

    df_fragilities = pd.DataFrame.from_dict(dict(ChainMap(*[{i:{
        'ShakingAndLSN':x[0][0],
        'Tsunami':x[0][1],
        'Landslide':x[1]}} for i,x in zip(np.arange(0,len(rows),1),rows)]))).T

    for index,row in df_fragilities.iterrows():
        comb_frag = np.max(row)
        df_fragilities.loc[index,'Combined'] = int(comb_frag)
        

    # Rename the columns to include fragility
    [df_fragilities.rename(columns={x:x+'Fragility'},inplace=True) for x in df_fragilities if x!='Combined']

    # Replace the integers with strings 
    [df_fragilities.replace(x,ds,inplace=True) for x,ds in zip([0,1,2,3,4],['DS0','DS1','DS2','DS3','DS4'])]

    df_fragilities.reset_index(drop=True,inplace=True)

    # Stack the dataframe to setup for creating CPTs
    df_fragilities_melt = df_fragilities.melt(['ShakingAndLSNFragility','TsunamiFragility','LandslideFragility'])
    df_fragilities_melt = df_fragilities_melt.sort_values(['ShakingAndLSNFragility','TsunamiFragility','LandslideFragility'])

    # Create CPT for fragilities and add to network
    cpt_list = []
    for index,row in df_fragilities_melt.iterrows():
        if row.value=='DS0':
            lst = [1,0,0,0,0]
        elif row.value=='DS1':
            lst = [0,1,0,0,0]
        elif row.value=='DS2':
            lst = [0,0,1,0,0]
        elif row.value=='DS3':
            lst = [0,0,0,1,0]
        elif row.value=='DS4':
            lst = [0,0,0,0,1]

        cpt_list = cpt_list+lst
    
    return(cpt_list)
    
    
def create_PGA_CTP():
    
    pga_general = np.arange(0,5,1)
    PGA1 = pga_general
    PGA2 = pga_general
    PGA3 = pga_general
    PGA4 = pga_general
    
    LSA = np.arange(0,4,1)

    rows = list(itertools.product(list(itertools.product(list(itertools.product(list(itertools.product(PGA1,PGA2)),PGA3)),PGA4)),LSA))

    df_pga_cpt = pd.DataFrame.from_dict(dict(ChainMap(*[{i:{
        'ShakingAndLSNFragilityA':x[0][0][0][0],
        'ShakingAndLSNFragilityB':x[0][0][0][1],
        'ShakingAndLSNFragilityC':x[0][0][1],
        'ShakingAndLSNFragilityD':x[0][1],
        'LiquefactionModulator':x[1]}} for i,x in zip(np.arange(0,len(rows),1),rows)]))).T

    df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==0,'ShakingAndLSNFragility'] = df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==0,'ShakingAndLSNFragilityA']
    df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==1,'ShakingAndLSNFragility'] = df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==1,'ShakingAndLSNFragilityB']
    df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==2,'ShakingAndLSNFragility'] = df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==2,'ShakingAndLSNFragilityC']
    df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==3,'ShakingAndLSNFragility'] = df_pga_cpt.loc[df_pga_cpt.LiquefactionModulator==3,'ShakingAndLSNFragilityD']

    [df_pga_cpt.replace(x,ds,inplace=True) for x,ds in zip([0,1,2,3,4],['DS0','DS1','DS2','DS3','DS4'])]
    [df_pga_cpt['LiquefactionModulator'].replace(x,ds,inplace=True) for x,ds in zip(['DS0','DS1','DS2','DS3'],['A','B','C','D'])]
    
    df_pga_cpt_melt = df_pga_cpt.melt(['ShakingAndLSNFragilityA','ShakingAndLSNFragilityB','ShakingAndLSNFragilityC','ShakingAndLSNFragilityD','LiquefactionModulator'])
    df_pga_cpt_melt = df_pga_cpt_melt.sort_values(['ShakingAndLSNFragilityA','ShakingAndLSNFragilityB','ShakingAndLSNFragilityC','ShakingAndLSNFragilityD','LiquefactionModulator'])
    
    cpt_list = []
    for index,row in df_pga_cpt_melt.iterrows():
        if row.value=='DS0':
            lst = [1,0,0,0,0]
        elif row.value=='DS1':
            lst = [0,1,0,0,0]
        elif row.value=='DS2':
            lst = [0,0,1,0,0]
        elif row.value=='DS3':
            lst = [0,0,0,1,0]
        elif row.value=='DS4':
            lst = [0,0,0,0,1]
        else:
            print(df_pga_cpt_melt)

        cpt_list = cpt_list+lst
    
    return(cpt_list)

def normalise(df_ds):
    df_ds["DS3"] = df_ds.DS3-df_ds.DS4
    df_ds["DS2"] = df_ds.DS2-df_ds.DS3-df_ds.DS4
    df_ds["DS1"] = df_ds.DS1-df_ds.DS2-df_ds.DS3-df_ds.DS4

    # Now add zero
    ds_dict = {}

    for index,row in df_ds.iterrows():
        row['DS0'] = (1-np.sum(row))
        ds_dict.update({
            index:row
        })

    df_ds = pd.DataFrame.from_dict(ds_dict,orient='index')
    
    return(df_ds)

def create_bins(factor,mu_list,sigma_list,ds_list,IM,bins):

    df_ds = pd.DataFrame.from_dict({ds:[lognorm(s=sigma,scale=math.exp(mu*factor)).cdf(im) for im in IM] for mu,sigma,ds in zip(mu_list,sigma_list,ds_list)})
#     df_ds = pd.DataFrame.from_dict({ds:[0.5*(1+math.erf((im/factor-mu)/(sigma*(0.5**.5)))) for im in IM] for mu,sigma,ds in zip(mu_list,sigma_list,ds_list)})
    df_ds.index = IM
    
    df_ds = normalise(df_ds)    
    df_ds = df_ds[['DS0','DS1','DS2','DS3','DS4']]
    
    P_dict = {x:min(list(df_ds.index), key=lambda y:abs(y-x)) for key,x in bins.items()} 
    cpt = sum([list(df_ds[df_ds.index==P_dict[x]].T[P_dict[x]]) for key,x in bins.items()],[])

    return(cpt)




    


    
