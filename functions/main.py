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
import models
import networkx as nx
import matplotlib.pyplot as plt



school_data = dl.load_schools_data()

data_net = mn.build_network_alinks()

vacc_dict_samples = vacc.vacc_sample_dict()

nodelist = np.intersect1d(data_net.nodes(), 
                          vacc_dict_samples['sample 0'].keys())

adj_mat_dt = nx.adjacency_matrix(data_net, weight='weight', 
                                 nodelist=nodelist).toarray()

nodes=np.intersect1d(nodelist, school_data.BRIN)



trans_mat = create_transmat_from_adjmat(adj_mat_dt, vacc_dict_samples['sample 1'], nodelist)

trans_net = create_network_from_transmat(trans_mat, nodelist)

reformed_codes = np.array(school_data.query('Denomination == "Reformatorisch"').BRIN)

ref_trans_net = trans_net.subgraph(nbunch=reformed_codes)

trans_net_loc = trans_net.subgraph(nbunch=nodes)
pos_geo = make_geo_pos(school_data)
degree_dict = trans_net.degree(weight='weight')

fig = plt.figure(figsize=[14,20])
ax = fig.add_subplot(111)
nx.draw_networkx_nodes(trans_net_loc, weight='weight', pos=pos_geo, node_size=np.array([degree_dict[n] for n in trans_net_loc.nodes()])*500., node_color='DodgerBlue', linewidths=0.2, ax=ax)
nx.draw_networkx_nodes(ref_trans_net, weight='weight', pos=pos_geo, node_size=np.array([degree_dict[n] for n in ref_trans_net.nodes()])*500., node_color='Red', linewidths=0.2, ax=ax)
nx.draw_networkx_edges(trans_net_loc, weight='weight', pos=pos_geo, width=np.log(nx.get_edge_attributes(trans_net_loc, 'weight').values())/5., edge_color='LightGrey', ax=ax)
ax.set_aspect('equal')

ax.set_xticks([])
ax.set_yticks([])

ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)


def create_transmat_from_adjmat(adj_mat, vacc_dict, nodelist):
    fs_vec = np.array([ models.op_fs(1.-vacc_dict[n], 15.) for n in nodelist])
    vac_vec = np.array([vacc_dict[n] for n in nodelist])
    fs_mat = np.outer(fs_vec, np.ones(len(fs_vec)))
    vacvact = np.outer( 1.-vac_vec, 1.-vac_vec)
    vac_mat = np.outer(np.ones(len(nodelist)), 1.-vac_vec)
    Pst_mat = 1. - np.minimum(1.,1./(15.*vac_mat))
    trans_mat = 1. - (1. - vacvact * fs_mat  * 0.5 * Pst_mat) ** adj_mat
    
    return trans_mat

def create_network_from_transmat(trans_mat, nodelist):
    trans_net = nx.from_numpy_matrix(trans_mat)
    rename_dict = dict([[i, nodelist[i]] for i in range(len(nodelist))])
    trans_net = nx.relabel_nodes(trans_net, rename_dict)
    return trans_net
    

def make_geo_pos(school_data):
    
    school_data['xy'] = list(zip(school_data['x_coord'], school_data['y_coord']))
    
    pos_geo = school_data[['BRIN','xy']].set_index('BRIN').to_dict()
    pos_geo = pos_geo['xy']
    
    

