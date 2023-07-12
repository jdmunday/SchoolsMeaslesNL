#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 12:00:13 2018

@author: LAPT0084
"""


import networkx as nx 
import pandas as pd
import numpy as np
import random 
import scipy as sp
from scipy import optimize

    
def model1(G, vacc_data, nsteps=100, q=0.5, R0=15., family_vacc=False, sample_from = 'ALL'):
    
    nodelist = G.nodes()
    
    ones_vec = np.ones(len(nodelist))
    
    
    adj_mat = nx.adjacency_matrix(G, weight='weight', nodelist=nodelist).toarray()
    vac_vec = np.array([1.-vacc_data[n] for n in nodelist])
    vac_mat = np.outer(ones_vec, vac_vec)
    vacvact = np.outer(vac_vec, vac_vec)
    
    
    
    Pst_mat = 1 - 1/(R0*vac_mat)
    
    op_fs_r = lambda s : op_fs(s, 15.)
    
    #fs_mat = np.transpose(np.array(map(op_fs_r, vac_vec)) * np.ones_like(adj_mat))
    fs_mat = np.outer(np.array([op_fs_r(v) for v in vac_vec]), ones_vec)
    
    
    
    if family_vacc == True: 
        vacvact = np.outer(vac_vec, np.ones(len(vac_vec)))
    
    trans_mat = 1. - (1. - vacvact * fs_mat  * q * Pst_mat) ** adj_mat
    #return trans_mat
    #trans_dec = np.reshape(map(decision, np.hstack(trans_mat)), (len(nodelist), len(nodelist)))
    
    trans_dec = np.random.random([len(nodelist),len(nodelist)]) < trans_mat
    
    
    
    
    init_vec = np.zeros_like(vac_vec)
    
    if sample_from == 'ALL': 
        seed_index = random.randrange(0,len(init_vec))
        
    elif sample_from == 'Vacc':
        cum_vac_vec = np.cumsum(vac_vec*vac_vec*(vac_vec > 0.1))
        ran = np.random.random()*max(cum_vac_vec)
        print(ran, max(vac_vec))
        seed_index = np.argmin(abs(cum_vac_vec - ran))
        
        print('seed susceptibility: ', vac_vec[seed_index])
        
    else: 
        random_index = random.randrange(0,len(sample_from))
        init_brin = sample_from[random_index]
        seed_index = nodelist.index(init_brin)

    init_vec[seed_index] = 1.
    
    
    
    
    
    s_vec = np.ones(len(nodelist)) - init_vec
    i_vec = init_vec
    r_vec = np.zeros(len(nodelist))
    
    s_vecs = [1.*s_vec]
    i_vecs = [1.*i_vec]
    r_vecs = [1.*r_vec]
    n = 0
    
    def append_vecs(S,I,R, s,i,r):
        S.append(1.*s)
        I.append(1.*i)
        R.append(1.*r)
    
    while n < 100 and sum(i_vec) > 0:
        print( n)
        n += 1
        r_vec += i_vec
        
        possible_trans = np.outer(i_vec, s_vec)*trans_dec
        
        #possible_trans = s_vec*np.transpose(i_vec*trans_dec.T)
        
        #plt.imshow(possible_trans)
        
        i_vec = np.array(sum(possible_trans) > 0 )
        s_vec -= i_vec
        
        append_vecs(s_vecs, i_vecs, r_vecs, s_vec, i_vec, r_vec)
        
        
        
    print(len(s_vec))
        
    return nodelist, s_vecs, i_vecs, r_vecs, init_vec
        
    
    
def decision(probability):
    return random.random() < probability


def op_fs(s, r0):
    rcalc = lambda s, r0, r_inf : r_inf - 1. +  np.exp(- s * r0 * r_inf)
    fun_op = lambda r : rcalc(s, r0, r) 
    fs = optimize.root(fun_op, 0.5).x[0]
    return fs * s
    




def model2(adj_mat, nodelist, vacc_data, nsteps=100, q=0.5, p_imp=0, R0=15., family_vacc=False, sample_from = 'ALL'):
    
    
    ones_vec = np.ones(len(nodelist))
    vac_vec = np.array([1.-vacc_data[n] for n in nodelist])
    
    if family_vacc == True:
        vac_vec_s = ones_vec
    else:
        vac_vec_s = vac_vec
    
    op_fs_r = lambda s : op_fs(s, 15.)
    
    
    
    init_vec = np.zeros_like(vac_vec)
    
    if sample_from == 'ALL': 
        seed_index = random.randrange(0,len(init_vec))
        
    elif sample_from == 'Vacc':
        cum_vac_vec = np.cumsum(vac_vec*vac_vec*(vac_vec > 0.1))
        ran = np.random.random()*max(cum_vac_vec)
        print( ran, max(vac_vec))
        seed_index = np.argmin(abs(cum_vac_vec - ran))
        
        print( 'seed susceptibility: ', vac_vec[seed_index])
    
    elif isinstance(sample_from, str):
        init_brin = sample_from
        print( 'seed set: ', init_brin)
        seed_index = np.where(nodelist == init_brin)[0]
        
        
    else: 
        random_index = random.randrange(0,len(sample_from))
        init_brin = sample_from[random_index]
        seed_index = np.where(nodelist == init_brin)[0]
        print( seed_index)

    init_vec[seed_index] = 1.
    
    s_vec = np.ones(len(nodelist)) - init_vec
    i_vec = init_vec
    r_vec = np.zeros(len(nodelist))
    
    s_vecs = [1.*s_vec]
    i_vecs = [1.*i_vec]
    r_vecs = [1.*r_vec]
    n = 0
    
    def append_vecs(S,I,R, s,i,r):
        S.append(1.*s)
        I.append(1.*i)
        R.append(1.*r)
    
    while n < 100 and sum(i_vec) > 0:
        print( n)
        n += 1
        r_vec += i_vec
        
        trans_bin = np.zeros([len(nodelist), len(nodelist)])
        
        for i in np.where(i_vec == 1.)[0]: 
            fs_val = 1.*op_fs_r(vac_vec[i]) 
            for j in np.where(s_vec == 1.)[0]:
                prob_not_link = 1.
                prob_not_impo = 1.
                if adj_mat[i,j] == 0 and p_imp != 0:
                    prob_not_impo = 1. - int(random.random() < vac_vec[j]) * (1 - 1/(R0*vac_vec[j])) * p_imp
                    
                if adj_mat[i,j] !=0:
                    prob_not_link = (1. - vac_vec_s[j] * vac_vec[i] * fs_val  * q * (1 - 1/(R0*vac_vec[j]))) ** adj_mat[i,j]
                
                trans_bin[i,j] = 1. - prob_not_link*prob_not_impo > np.random.random() 
        
        i_vec = np.array(sum(trans_bin) > 0 )
        s_vec -= i_vec
        
        append_vecs(s_vecs, i_vecs, r_vecs, s_vec, i_vec, r_vec)
                
    return nodelist, s_vecs, i_vecs, r_vecs, init_vec

def model3(adj_mat, nodelist, vacc_data, nsteps=100, q=0.5, p_imp=0, R0=15., family_vacc=False, sample_from = 'ALL'):
    
    
    ones_vec = np.ones(len(nodelist))
    vac_rel = np.random.randint(0, len(vacc_data)-1)
    vac_vec = np.array([1.-vacc_data['sample ' + str(vac_rel)][n] for n in nodelist])
    
    if family_vacc == True:
        vac_vec_s = ones_vec
    else:
        vac_vec_s = vac_vec
    
    op_fs_r = lambda s : op_fs(s, 15.)
    
    
    
    init_vec = np.zeros_like(vac_vec)
    
    if sample_from == 'ALL': 
        seed_index = random.randrange(0,len(init_vec))
        
    elif sample_from == 'Vacc':
        cum_vac_vec = np.cumsum(vac_vec*vac_vec*(vac_vec > 0.1))
        ran = np.random.random()*max(cum_vac_vec)
        print( ran, max(vac_vec))
        seed_index = np.argmin(abs(cum_vac_vec - ran))
        
        print( 'seed susceptibility: ', vac_vec[seed_index])
    
    elif isinstance(sample_from, str):
        init_brin = sample_from
        print( 'seed set: ', init_brin)
        seed_index = np.where(nodelist == init_brin)[0]
        
        
    else: 
        random_index = random.randrange(0,len(sample_from))
        init_brin = sample_from[random_index]
        seed_index = np.where(nodelist == init_brin)[0]
        print( seed_index)

    init_vec[seed_index] = 1.
    
    s_vec = np.ones(len(nodelist)) - init_vec
    i_vec = init_vec
    r_vec = np.zeros(len(nodelist))
    
    s_vecs = [1.*s_vec]
    i_vecs = [1.*i_vec]
    r_vecs = [1.*r_vec]
    n = 0
    
    def append_vecs(S,I,R, s,i,r):
        S.append(1.*s)
        I.append(1.*i)
        R.append(1.*r)
    
    while n < 100 and sum(i_vec) > 0:
        print( n)
        n += 1
        r_vec += i_vec
        
        trans_bin = np.zeros([len(nodelist), len(nodelist)])
        
        for i in np.where(i_vec == 1.)[0]: 
            fs_val = 1.*op_fs_r(vac_vec[i]) 
            for j in np.where(s_vec == 1.)[0]:
                prob_not_link = 1.
                prob_not_impo = 1.
                if adj_mat[i,j] == 0 and p_imp != 0:
                    prob_not_impo = 1. - int(random.random() < vac_vec[j]) * (1 - 1/(R0*vac_vec[j])) * p_imp
                    
                if adj_mat[i,j] !=0:
                    prob_not_link = (1. - vac_vec_s[j] * vac_vec[i] * fs_val  * q * (1 - 1/(R0*vac_vec[j]))) ** adj_mat[i,j]
                
                trans_bin[i,j] = 1. - prob_not_link*prob_not_impo > np.random.random() 
        
        i_vec = np.array(sum(trans_bin) > 0 )
        s_vec -= i_vec
        
        append_vecs(s_vecs, i_vecs, r_vecs, s_vec, i_vec, r_vec)
                
    return nodelist, s_vecs, i_vecs, r_vecs, init_vec
    
    

      
def monty(model_choice, n, params):
    res = []
    for n in xrange(n):
        res.append(model_choice(*params))
    return res
    
  
