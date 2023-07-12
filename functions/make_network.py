#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 09:10:58 2018

@author: LAPT0084
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 14:49:57 2018

@author: LAPT0084
"""


import networkx as nx 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import os

import create_secpri as cs
import data_loader as dl


basepath = os.path.dirname(__file__)
schools_path = basepath +  "/../data/"

def build_school_network(): 
    secpri = cs.create_secpri_multi()
    schools_net = make_network(secpri)
    return schools_net



def build_network_alinks(data_path = schools_path + r'/actual_links.csv'):

    alinks = pd.read_csv(data_path)
        
    data_net = nx.Graph()
    
    
    edges = np.array(alinks[['BRINVEST1', 'BRINVEST2']])
    weights = [float(w) for w in alinks.aantal]
    
    
    for n in set(list(alinks.BRINVEST1) + list(alinks.BRINVEST2)):
        data_net.add_node(n)
        
    for i, edge in enumerate(edges):
        c_ij = min(500., weights[i])
        data_net.add_edge(*edge, weight=c_ij , inv_weight=1./c_ij)
        
        
    return data_net





def make_network(school_frame, HH_dist_frame = 'default'):
    
    schools_net = nx.Graph()
    
    secondary_urns = school_frame.SecondaryURN.unique()
    #secondary_names = secpri['INSTELLINGSNAAM VO'].unique()
    schools_net.add_nodes_from(secondary_urns, typ='Sec')
    
    primary_urns = school_frame.Primary_School_URN.unique()
    #primary_names = school_frame['INSTELLINGSNAAM PO'].unique()
    schools_net.add_nodes_from(primary_urns, typ='Pri')
    
    edges = np.array(school_frame[['SecondaryURN', 'Primary_School_URN']])
    weights = np.array([float(w) for w in school_frame.pupil_count])
    for i, edge in enumerate(edges):
        if isinstance(HH_dist_frame, str):
            c_ij = edge_weight(weight=weights[i], ave_age_diff=3.)
        else:
            c_ij = edge_weight(weight=weights[i], ave_age_diff=3. , prop_c = np.array(HH_dist_frame.loc[edge[1]]))
        schools_net.add_edge(*edge, weight = c_ij , inv_weight=c_ij, pup_count=weights[i], inv_pup_count=weights[i])
        
    
    return schools_net

def edge_weight(weight=1., ave_age_diff=1., prop_c = np.array([0.2, 0.3, 0.35, 0.1, 0.05])):
    
    rho = lambda n: sum([k * (n-k) for k in np.arange(1, n)])/((n-1.)+0.0000000000000000001)
    nobs = lambda n: prop_c[n-1] * (n - 1) / n
    weight_new = weight * ave_age_diff * sum([rho(n) * nobs(n) for n in  np.arange(1, int(len(prop_c) + 1))])
    
    
    return weight_new


def geo_dist(data, node1, node2):
    
    return np.array(data.query("BRIN == @node1").geometry)[0].distance(np.array(data.query("BRIN == @node2").geometry)[0])
    


def make_random_network(net, school_frame, school_data):
    
    ran_net = nx.Graph()
    
    secondary_urns = [s for s in school_frame.SecondaryURN.unique() if s in map(str, school_data.BRIN)]

    ran_net.add_nodes_from(secondary_urns, typ='Sec')
    
    primary_urns = [p for p in school_frame.Primary_School_URN.unique() if p in map(str, school_data.BRIN)]
    #primary_names = school_frame['INSTELLINGSNAAM PO'].unique()
    ran_net.add_nodes_from(primary_urns, typ='Pri')
    
    for sec in secondary_urns[:1]:            
        
        link_weights = np.array([[pri , net.degree(pri, weight='weight')/geo_dist(school_data, sec, pri)] for pri in primary_urns])
        #link_weights = link_weights / sum(link_weights)
        
        lw_n = link_weights[:,0]
        lw_w = np.array(map(float, link_weights[:,1]))
        
        lw_w_edge = net.degree(sec, weight='weight') * lw_w / sum(lw_w)
        
        for i, p in enumerate(lw_n):
            ran_net.add_edge(p, sec, weight=lw_w_edge[i])
        
        
        
        
    return ran_net
        

def make_random_adjmat(net, school_frame, school_data, p=1., rule = 'power', sec_pri_only=True):
    
    school_data_match = school_data.query('BRIN in @net.nodes()')
    
    
    
    print( len(school_data_match))
    
    
    nodelist = np.array(school_data_match.BRIN)
    
    xs = np.array(school_data_match.x_coord)
    ys = np.array(school_data_match.y_coord)
    
    xmat = np.array([list(xs)]*len(xs)).T
    ymat = np.array([list(ys)]*len(ys)).T
    
    dist_mat = ((xmat - xs)**2 + (ymat - ys)**2)**0.5
    if rule == 'power': 
        inv_dist_mat = 1./dist_mat**p
        inv_dist_mat[np.where(inv_dist_mat == inv_dist_mat[0,0])] = 0.
    
    elif rule == 'gauss': 
        inv_dist_mat = np.exp(-1. * (dist_mat/1000.) ** p)
        inv_dist_mat = inv_dist_mat * -(np.identity(len(nodelist)) - 1.)
        #return dist_mat, inv_dist_mat
        
    else:
        print( 'enter valid rule: "power" or "gauss"')
        return 0, 0, 0
        
    
    
    degs = [net.degree(k, weight='weight') for k in school_data_match.BRIN]
    degmat = np.array([degs]*len(degs))
    deg_scaled_mat = inv_dist_mat * degmat
    
    if sec_pri_only == True:
    
        kinds = np.array(school_data_match.kind)
        kinds_mat = np.array([list(kinds)]*len(kinds)).T
        sec_pri_link = kinds_mat != kinds
        
        deg_scaled_mat_sp = sec_pri_link * deg_scaled_mat 
        
         
        
        contact_mat = degmat * deg_scaled_mat_sp/sum(deg_scaled_mat_sp)
    
    else: 
        contact_mat = deg_scaled_mat
    
    
    
        
    return contact_mat, nodelist, dist_mat


def compare_nets(adjmat, net, school_frame, school_data, rule='gauss', ps=np.arange(0.5, 0.85, 0.05)):
    aj_mats = []
    for p in ps:
        contact_mat = make_random_adjmat(net, school_frame, school_data, p, 'gauss')
        dist_dist_c = sum(contact_mat[0] * contact_mat[2]) / (sum(contact_mat[0])+1.e-15)
        h, b = np.histogram(dist_dist_c, bins=range(0,20000, 500))
        plt.plot(b[:-1],h, alpha=0.9)
        aj_mats.append(contact_mat)
        
    dist_dist_d = sum(adjmat * contact_mat[2]) / (sum(adjmat)+1.e-15)
    
    hd, bd = np.histogram(dist_dist_d, bins=range(0,20000, 500))
    
    plt.plot(bd[:-1],hd,'kx-', alpha=0.8)
    
    plt.legend(labels=list(ps) + ['data'])
    
    return aj_mats
    
    
    
