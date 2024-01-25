#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 12:15:40 2018

@author: LAPT0084
"""

import numpy as np
import networkx as nx

import data_loader as dl
import vaccination as vacc
import make_network as mn
import network_analysis as na
import models
import matplotlib.pyplot as plt
import models
import importlib
import igraph as ig

importlib.reload(models)

def create_transmat_from_adjmat_asym(adj_mat, vacc_dict=False, nodelist=[], R0=15., q=0.5):
    if vacc_dict == False:
        vacc_dict = dict(np.transpose([nodelist, np.zeros_like(nodelist)]))
    #print 'making fsvec'
    FSLIST = [max(0,models.op_fs(p, 15.))for p in np.arange(0, 1., 0.001)]
  
                         
    #fs_vec = np.array([ models.op_fs(1.-vacc_dict[n], R0) for n in nodelist])
    fs_vec = np.array([FSLIST[max(0, int((1.-vacc_dict[n])*1000) - 1)] for n in nodelist])
   # print 'making vacvec'
    vac_vec = np.array([vacc_dict[n] for n in nodelist])

 
   # print 'making Pstvec'
    Pst_vec = 1. - np.minimum(1.,1./(R0*(1.-vac_vec)))
    #Pst_mat = np.outer(np.ones(len(nodelist)), 1.-Pst_vec)
    #print 'making pairprob'
    pairprob = np.outer( q*fs_vec, (1.-vac_vec)*Pst_vec)
    #print 'making transmat'
    trans_mat = 1. - (1. - pairprob) ** adj_mat
    
    return trans_mat


def create_transmat_from_adjmat(adj_mat, vacc_dict=False, nodelist=[], R0=15., q=0.5):
    if vacc_dict == False:
        vacc_dict = dict(np.transpose([nodelist, np.zeros_like(nodelist)]))
    #print 'making fsvec'
    FSLIST = [max(0,models.op_fs(p, R0))for p in np.arange(0, 1., 0.001)]
  
                         
    #fs_vec = np.array([ models.op_fs(1.-vacc_dict[n], R0) for n in nodelist])
    fs_vec = np.array([FSLIST[max(0, int((1.-vacc_dict[n])*1000) - 1)] for n in nodelist])
   # print 'making vacvec'
    vac_vec = np.array([vacc_dict[n] for n in nodelist])

 
   # print 'making Pstvec'
    Pst_vec = fs_vec  # equivalent to the final size for gaussian offspring distribution
    #Pst_mat = np.outer(np.ones(len(nodelist)), 1.-Pst_vec)
    #print 'making pairprob'
    pairprob = np.outer( q*fs_vec, (1.-vac_vec)*Pst_vec)
    #print 'making transmat'
    trans_mat = 1. - (1. - pairprob) ** adj_mat
    
    return trans_mat



def create_network_from_transmat(trans_mat, nodelist):
    trans_net = nx.from_numpy_array(trans_mat, create_using=nx.MultiDiGraph())
    rename_dict = dict([[i, nodelist[i]] for i in range(len(nodelist))])
    trans_net = nx.relabel_nodes(trans_net, rename_dict)
    return trans_net


def create_network_from_transmat_undi(trans_mat, nodelist):
    trans_net = nx.from_numpy_array(trans_mat)
    rename_dict = dict([[i, nodelist[i]] for i in range(len(nodelist))])
    trans_net = nx.relabel_nodes(trans_net, rename_dict)
    return trans_net

def create_network_from_transmat_undi_ig(trans_mat, nodelist):

    # Create graph, A.astype(bool).tolist() or (A / A).tolist() can also be used.
    trans_net = ig.Graph.Adjacency((trans_mat > 0).tolist())
    
    # Add edge weights and node labels.
    trans_net.vs['label'] = nodelist 
    
    return trans_net
    

def make_geo_pos(school_data):
    
    school_data['xy'] = list(zip(school_data['x_coord'], school_data['y_coord']))
    
    pos_geo = school_data[['BRIN','xy']].set_index('BRIN').to_dict()
    pos_geo = pos_geo['xy']
    
    return pos_geo


def plot_transnet(trans_net, pos, highlight=[] , title='transmission network of the Nethelands', sizes='deg', ax=False):
    
    if ax == False:
    
        fig = plt.figure(figsize=[14,20])
        ax = fig.add_subplot(111)
    
    
    
    highlight_net = trans_net.subgraph(nodes=highlight)
    if sizes == 'deg':
        degree_dict = trans_net.degree(weight='weight')
        sizes = degree_dict
    else:
        sizes
    edge_dict = nx.get_edge_attributes(trans_net, 'weight')
    #return edge_dict
    #nx.draw_networkx_nodes(trans_net, weight='weight', pos=pos, node_size=10, node_color='Grey', alpha=0.3, linewidths=0.2, ax=ax)
    nx.draw_networkx_nodes(trans_net, pos=pos, node_size=np.array([sizes[n] for n in trans_net.nodes()])*1., node_color='DodgerBlue', alpha=0.7, linewidths=0.2, ax=ax)
    nx.draw_networkx_nodes(highlight_net, pos=pos, node_size=np.array([sizes[n] for n in highlight_net.nodes()])*1., node_color='Red', alpha=0.7, linewidths=0.2, ax=ax)
    nx.draw_networkx_edges(trans_net,  pos=pos, edgelist=trans_net.edges(), width=np.array([edge_dict[e + (0,) ] for e in trans_net.edges()]), edge_color='Grey', ax=ax, arrows=False)
    ax.set_aspect('equal')
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.set_title(title)
    
    return ax

school_data = dl.load_schools_data()

data_net = mn.build_network_alinks()

vacc_dict_samples = vacc.vacc_sample_dict()

nodelist = np.intersect1d(data_net.nodes(), 
                          list(map(str,vacc_dict_samples['sample 0'].keys())))

adj_mat_dt = nx.adjacency_matrix(data_net, weight='weight', 
                                 nodelist=nodelist).toarray()

nodes=np.intersect1d(nodelist, school_data.BRIN)



trans_mat = create_transmat_from_adjmat(adj_mat_dt, vacc_dict_samples['sample 0'], nodelist)

trans_net = create_network_from_transmat(trans_mat, nodelist)

reformed_codes = np.array(school_data.query('Denomination == "Reformatorisch"').BRIN)

ref_trans_net = trans_net.subgraph(nodes=reformed_codes)

trans_net_loc = trans_net.subgraph(nodes=nodes)
pos_geo = make_geo_pos(school_data)
degree_dict = trans_net.degree(weight='weight')

#plot_transnet(trans_net_loc, pos_geo, reformed_codes)


def create_yso_dict(years_since_ob, vacc_dict, school_data, years_in_pri=8., years_in_sec=5., post_outbreak=0.98):
    vacc_dict_ = vacc_dict.copy()
    ps = min(years_since_ob, years_in_pri)
    
    ss = min(max(0, years_since_ob - years_in_pri), years_in_sec)
    
    for k in nodes:
        if school_data.query('BRIN == @k').kind.values[0] == 'pri': 
            vacc_dict_[k] = (vacc_dict_[k] * ps/years_in_pri + post_outbreak * (1-ps/years_in_pri)) 
            
        else:
            vacc_dict_[k] = (vacc_dict_[k] * ss/years_in_sec + post_outbreak * (1-ss/years_in_sec))
        
    return vacc_dict_
    

def make_yso_set(adj_mat, vacc_dict, years_since_ob, school_data, nodelist, years_in_pri=7., years_in_sec=7., post_outbreak=0.98):
    vacc_dict_corr = create_yso_dict(years_since_ob, dict(vacc_dict), school_data, post_outbreak=post_outbreak)
    trans_mat_yso = create_transmat_from_adjmat(adj_mat, vacc_dict_corr, nodelist)
    trans_net_yso = create_network_from_transmat(trans_mat_yso, nodelist)
    nodes = np.intersect1d(nodelist, school_data.BRIN)
    trans_net_yso = trans_net_yso.subgraph(nodes=nodes)
    return trans_net_yso


def calc_probable_paths(net, nodes, path_weight='inv_weight', return_weight='weight'):

    distances_full = []
    
    weights = nx.get_edge_attributes(net, return_weight)
    dap = lambda n : distances_full.append(n)
    
    for i, node1 in enumerate(nodes): 
        for node2 in  nodes[i:]:
            try: 
                network_dist = nx.shortest_path(net, node1, node2, weight=path_weight)
                dap(network_dist)
            except: 
                foo = 0
    
    
    
    
    sum_distance = []
    prod_distance = []
    
    for d in distances_full: 
        if len(d) > 1:
            d_w = []
            for i in range(len(d)-1):
                try:
                    d_w.append(weights[(d[i+1],d[i])])
                    #print weights[(d[i+1],d[i])]
                except:
                    d_w.append(weights[(d[i],d[i+1])])
            sum_distance.append(sum(d_w))
            prod_distance.append(np.product(d_w))
    
    
    return  sum_distance, prod_distance
    

def calculate_connections(adj_mat, vacc_dict, school_data, nodes):
    
    yso=[]
    for n in range(0,15, 1):
        yso.append(make_yso_set(adj_mat, vacc_dict, n, school_data, nodes))
    
    
    sum_dist = []
    prod_dist = []
    for y in yso:
        weights = nx.get_edge_attributes(y, 'weight')
        inv_weights = dict([[k, 1. - weights[k]] for k in weights.keys()])
        nx.set_edge_attributes(y, 'inv_weight', inv_weights)
        sum_distance_no, prod_distance = calc_probable_paths(y, reformed_codes, return_weight='weight')
        sum_distance, prod_distance_no = calc_probable_paths(y, reformed_codes, return_weight='inv_weight')
        sum_dist.append(sum_distance)
        prod_dist.append(prod_distance)
    
    for i, p in enumerate(prod_dist[2:][::2]):
        h, b = np.histogram(p, bins=np.arange(0,1.0,0.1))
        plt.plot(b[:-1], h, '-o', color=plt.cm.viridis(40*i))
    plt.yscale('log')
    plt.legend(['{} years post outbreak'.format(n) for n in range(2,15,2)], bbox_to_anchor = (1.,1.0))
    plt.xlabel('Probability of completing path')
    
    
    
def yso_outbreak(adj_mat, years_since_ob, vacc_dict, school_data, years_in_pri, years_in_sec, post_outbreak, nodelist, seed_list='Vacc'): 
    res = []
    for n in xrange(20):
        n = np.random.randint(2499)
        vacc_dict_s = vacc_dict['sample {}'.format(n)]
        vacc_dict_yso = create_yso_dict(years_since_ob, vacc_dict_s, school_data, years_in_pri, years_in_sec, post_outbreak)
        res.append(models.model2(adj_mat, nodelist, vacc_dict_yso, nsteps=100, q=0.5, p_imp=0, R0=15., family_vacc=False, sample_from = seed_list))
    return res


def phase_change(adj_mat, vacc_dict, school_data, nodes):
    
    yso=[]
    for n in range(0,15, 1):
        yso.append(make_yso_set(adj_mat, vacc_dict, n, school_data, nodes))
    return yso


def find_outbreak(trans_mat, seed ='03LQ00', nodelist=nodelist):
    
    obnet = create_network_from_transmat(trans_mat > np.random.random(trans_mat.shape)*1., nodelist=nodelist)
    return obnet
    con_list = []
    nod_list = [seed]
    while len(nod_list) > 0:
        
        con_list.extend(nod_list)
        nod_list = list(np.hstack([obnet.predecessors(n) for n in nod_list]))
        nod_list = [n for n in nod_list if n not in con_list]
    
    return len(set(con_list))
    
    
    

