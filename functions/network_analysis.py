#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 18:36:11 2018

@author: LAPT0084
"""
 
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point



def calculate_homophily(schools_net, schools_data_gdf):

    
    HI = []
    
    for den in schools_data_gdf.Denomination.unique():    
        nodes = schools_data_gdf.query('Denomination == @den').BRIN
        
        si = np.mean([d for n,d in schools_net.subgraph(nodes).degree(weight='weight')])
        all_degrees = schools_net.degree(weight='weight')
        di = np.mean([all_degrees[n] for n in nodes if n in schools_net.nodes()]) - si 
        
        N = schools_net.number_of_nodes()
        Ni = schools_net.subgraph(nodes).number_of_nodes()
        
        wi = 1.*Ni/N
        
        Hi = si/(si+di)
        IHi = (Hi - wi)/(1-wi)
        
        HI.append([den, wi, Hi, IHi])
    
    IH_df = pd.DataFrame(HI, columns=['Denomination',  'wi', 'Hi',  'IH'])
    
    #IH_df = IH_df[np.array(IH_df.IH) > 0.1]
    
    
    fig = plt.figure(figsize=[10,5])
    ax=fig.add_subplot(111)
    IH_df.plot(kind='bar', x='Denomination', y = 'IH', ax=ax, color='DarkOrange')
    ax.set_ylabel('Coleman Index', fontsize=14)
    ax.set_xlabel('Denomination', fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    
    fig = plt.figure(figsize=[20,10])
    ax=fig.add_subplot(111)
    IH_df.plot(kind='bar', x='Denomination', y = 'wi', ax=ax)
    ax.set_ylabel('Inbreeding Homophily Index', fontsize=20)
    ax.set_xlabel('Denomination', fontsize=20)
    #ax.set_xticklabels(fontsize=20)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return IH_df
    

def distance_compare(net, data, nodes, weight='inv_weight', sample='ALL', prepairs=0):
    
    

    distances_full = []
    
    weights = nx.get_edge_attributes(net, weight)
    
    #return weights

    dap = lambda n : distances_full.append(n)
    if sample=='ALL': 
    
        for i, node1 in enumerate(nodes): 
            for node2 in  nodes[i:]:
                try: 
                    network_dist = nx.shortest_path(net, node1, node2, weight=weight)
                    geo_dist = 0.001*np.array(data.query("BRIN == @node1").geometry)[0].distance(np.array(data.query("BRIN == @node2").geometry)[0])
                    dap([network_dist, geo_dist])
                    
                except: 
                    foo = 0
        
    elif isinstance(sample,(int, float)) or sample == 'pairs':
        
        
        
        if prepairs==0:
            print( 'selecting pairs')
            pairmatch = 1
            pairs=[]
            for i in range(sample):
    
                pair = np.random.choice(nodes, 2 , replace=False)
                while pairmatch == 1:
                    if pair in pairs or pair[::-1] in pairs:
                        pairmatch=1
                    else:
                        pairmatch=0
                pairs.append(pair)

        
        else:
            
            print( 'pairs provided')
            pairs=prepairs
        

        
        for pair in pairs:
            try:
                node1 = pair[0] 
                node2 = pair[1]
                network_dist = nx.shortest_path(net, node1, node2, weight=weight)
                geo_dist = 0.001*np.array(data.query("BRIN == @node1").geometry)[0].distance(np.array(data.query("BRIN == @node2").geometry)[0])
                dap([network_dist, geo_dist])
                    
            except: 
                foo = 0
            
    
    unweighted_distance = []
    weighted_distance = []
    
    for d in np.asarray(distances_full, dtype=object).T[0]: 
        unweighted_distance.append(len(d))
        if len(d) > 1:
            d_w = []
            for i in range(len(d)-1):
                try:
                    d_w.append(weights[(d[i+1],d[i])])

                except:
                    d_w.append(weights[(d[i],d[i+1])])
            weighted_distance.append(sum(d_w))
        else: 
            weighted_distance.append(0)
    
    #plt.plot(np.array(distances_full).T[1], weighted_distance, 'o')
    
    #plt.plot(np.array(distances_full).T[1], unweighted_distance, 'o')
    
    return  np.asarray(distances_full, dtype=object).T[1], weighted_distance, unweighted_distance

        

def compare_denom(net, data, den, sample=30):
    sdata = data
    sdata['pc2'] = [int(str(p)[:2]) for p in sdata['pc6']]
    
    Rnodes = np.random.choice(np.array(sdata.query('Denomination == @den').BRIN), sample)
    pc2s = sdata.query('BRIN in @Rnodes').groupby('pc2').count().reset_index()
    pc2s_all = sdata.query('Denomination == @den').groupby('pc2').count().reset_index()
    
    comp_schools = []
    for p in pc2s_all['pc2']:
        samp_pc2 = pc2s_all.query('pc2 == @p').BRIN.values[0]
        comp_schools.extend(np.random.choice(np.array(sdata.query('Denomination != @den and pc2 == @p').BRIN), samp_pc2, replace=False))
    Anodes = []
    for p in pc2s['pc2']:
        samp_pc2 = pc2s.query('pc2 == @p').BRIN.values[0]
        Anodes.extend(np.random.choice(np.array(sdata.query('Denomination != @den and pc2 == @p').BRIN), samp_pc2, replace=False))
        
    
    
    
    gd_den, wd_den, ud_den = distance_compare(net.subgraph([n for n in net.nodes() if n not in comp_schools]), sdata, Rnodes, weight='inv_weight')
    gd_res, wd_res, ud_res = distance_compare(net.subgraph(sdata.query('Denomination != @den').BRIN), sdata, Anodes, weight='inv_weight')
#    fig = plt.figure()
#    plt.plot(gd_res, wd_res,'o', color='DodgerBlue', alpha=0.7)
#    plt.plot(gd_den, wd_den,'o', color='red', alpha=0.7)
#    
#    plt.title(den)
#    plt.ylabel('netowork distance')
#    plt.xlabel('geographic distance (km)')
#    
#    distance_ratios = np.array(wd_res)/(np.array(gd_res) + 1.e-27)
#    distance_ratios_den = np.array(wd_den)/(np.array(gd_den) + 1.e-27)
#    
#    fig = plt.figure()
#    plt.hist(distance_ratios_den[30:], bins=np.arange(0,3.e-2, 6.e-4), histtype='step', color='r', linewidth=2.)
#    plt.hist(distance_ratios, bins=np.arange(0,3.e-2, 6.e-4), histtype='step', color='DodgerBlue', linewidth=2.)
#    
#    plt.title(den)
#    plt.ylabel('pairs of schools')
#    plt.xlabel('Network - Geographic distance ratio')
#    
    return gd_den, wd_den, ud_den, gd_res, wd_res, ud_res

def compare_denom_prepairs(net, data, den, sample=30):
    sdata = data
    sdata['pc2'] = [int(str(p)[:2]) for p in sdata['pc6']]

    
    nodes = np.array(data.query('Denomination == @den').BRIN)
    
    pairmatch = 1
    pairs=[]
    for i in range(sample):

        pair = np.random.choice(nodes, 2 , replace=False)
        while pairmatch == 1:
            
            if pair in pairs or pair[::-1] in pairs:
                pairmatch=1
            else:
                pairmatch=0
        pairs.append(pair)
    
    print( 'pairs generated')
    
    comp_pairs = []
    pairmatch = 1
    for pair in pairs:

        pc2s = sdata.query('BRIN in @pair')['pc2']
        comp_pair = [np.random.choice(np.array(sdata.query('Denomination != @d and pc2 == @p').BRIN),1)[0] for p in pc2s if (d := den)]

        while pairmatch == 1:
            if comp_pair in comp_pairs or comp_pair[::-1] in comp_pairs:
                pairmatch=1
            else:
                pairmatch=0
        comp_pairs.append(comp_pair)
    
    print( 'con_pairs generated')
    
    pc2s_all = sdata.query('Denomination == @den').groupby('pc2').count().reset_index()
    
    comp_schools = []
    for p in pc2s_all['pc2']:
        samp_pc2 = pc2s_all.query('pc2 == @p').BRIN.values[0]
        comp_schools.extend(np.random.choice(np.array(sdata.query('Denomination != @den and pc2 == @p').BRIN), samp_pc2, replace=False))
    
    
    
    gd_den, wd_den, ud_den = distance_compare(net.subgraph([n for n in net.nodes() if n not in comp_schools]), sdata, nodes, weight='inv_weight', sample='pairs', prepairs=pairs)
    gd_res, wd_res, ud_res = distance_compare(net.subgraph(sdata.query('Denomination != @den').BRIN), sdata, nodes, weight='inv_weight', sample='pairs', prepairs=comp_pairs)
    fig = plt.figure()
    plt.plot(gd_res, wd_res,'o', color='DodgerBlue', alpha=0.7)
    plt.plot(gd_den, wd_den,'o', color='red', alpha=0.7)
    
    plt.title(den)
    plt.ylabel('netowork distance')
    plt.xlabel('geographic distance (km)')
    
    distance_ratios = np.array(wd_res)/(np.array(gd_res) + 1.e-27)
    distance_ratios_den = np.array(wd_den)/(np.array(gd_den) + 1.e-27)
    
    fig = plt.figure()
    plt.hist(distance_ratios_den[30:], bins=np.arange(0,3.e-2, 6.e-4), histtype='step', color='r', linewidth=2.)
    plt.hist(distance_ratios, bins=np.arange(0,3.e-2, 6.e-4), histtype='step', color='DodgerBlue', linewidth=2.)
    
    plt.title(den)
    plt.ylabel('pairs of schools')
    plt.xlabel('Network - Geographic distance ratio')
    
    return gd_den, wd_den, ud_den, gd_res, wd_res, ud_res

    
        