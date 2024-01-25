#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 20:06:19 2018

@author: LAPT0084
"""

import pandas as pd 
import numpy as np
import networkx as nx
import data_loader as dl
import os

base_path = os.getcwd()

VACC_DATA_PATH = base_path + r'/data/Vaccination data'
PO_VACC_FILE = r'POposterior_infF_denomT_JAGS.csv'
VO_VACC_FILE = r'VOposterior_infF_denomT_JAGS.csv'

def convert_objects(data_frame, convert_numeric=True):
    new_dataframe = pd.DataFrame([pd.to_numeric(data_frame[col], errors='ignore') for col in data_frame.columns]).T
    new_dataframe.columns = data_frame.columns
    return new_dataframe

def vacc_sample_dict():
    
    PO_data = pd.read_csv(VACC_DATA_PATH + r'/' + PO_VACC_FILE)
    PO_data.columns=['school'] + ['sample ' + str(n) for n in range(2500)]
    VO_data = pd.read_csv(VACC_DATA_PATH + r'/' + VO_VACC_FILE)
    VO_data.columns=['school'] + ['sample ' + str(n) for n in range(2500)]
    
    PO_data.school = [s.replace('_', '0') if s[-2] == '_' else  
                      s.replace('_', '') for s in PO_data.school]
    VO_data.school = [s.replace('_', '0') if s[-2] == '_' else  
                      s.replace('_', '') for s in VO_data.school]
    
    ALL_data = pd.concat([PO_data,VO_data])
    
    vacc_dict_samples = ALL_data.set_index('school').to_dict()
    
    return vacc_dict_samples


def vacc_pc4_dict(nodelist, vacc_dict, n=20., by='school'):
    
    '''
    
    requires: 
        nodelist - list of schools to be included in dictionary
        vacc_dict - school level vaccination dictionary to calculate pc4 vacc
                    e.g. vacc_dict_samples['sample 1']
                    
    returns: 
        dictionary of vaccination coverage in schools infered from pc4 level 
        vaccination coverage
        
    '''
    

    schools_pc4 = dl.load_schools_2_pc4()
    
    vdf = pd.DataFrame.from_dict(vacc_dict, orient='Index').reset_index()
    vdf.columns = ['BRIN', 'vacc']
    
    schools_pc4_vacc = schools_pc4.merge(vdf, on='BRIN')
    schools_pc4_vacc['vaccpc'] = schools_pc4_vacc.vacc * schools_pc4_vacc.TOTAAL
    schools_pc4_vacc_sums = schools_pc4_vacc.groupby('PC4').sum()
    schools_pc4_vacc_sums['overallvacc'] = schools_pc4_vacc_sums.vaccpc / schools_pc4_vacc_sums.TOTAAL
    
    
    
    PC_Vacc = schools_pc4_vacc_sums.reset_index()[['PC4', 'overallvacc']]
    
    if by == 'pc4':
        return PC_Vacc
    
    schools_pc4_vacc = schools_pc4_vacc.merge(PC_Vacc, on='PC4')
    schools_pc4_vacc = convert_objects(schools_pc4_vacc, convert_numeric=True)
    
    schools_pc4_vacc['vaccpc_2'] = schools_pc4_vacc.overallvacc * schools_pc4_vacc.TOTAAL
    schools_pc4_vacc_sums = schools_pc4_vacc.groupby('BRIN').sum()
    schools_pc4_vacc_sums['overallvacc'] = schools_pc4_vacc_sums.vaccpc_2 / schools_pc4_vacc_sums.TOTAAL
    school_Vacc = schools_pc4_vacc_sums.reset_index()[['BRIN', 'overallvacc']]
    
    intersection = list(filter(lambda x: x  not in list(school_Vacc.BRIN), nodelist))
    
    vacc_dict_pc4 = pd.concat([school_Vacc, pd.DataFrame(np.transpose([intersection, [0.98]*len(intersection)]), columns=['BRIN', 'overallvacc'])], ignore_index=True).set_index('BRIN')
    vacc_dict_pc4 = convert_objects(vacc_dict_pc4, convert_numeric=True)
    vacc_dict_pc4 = vacc_dict_pc4.fillna(0.98).to_dict()['overallvacc']
    return vacc_dict_pc4





# --------- generate vaccination dictionary from mean primary school vacc------
#---------- status - superceeded by vacc_sample_dict---------------------------





def find_neighbours(schools_net, schools_data_gdf):
    
    neighbours = []
    nodes_op = schools_data_gdf.query('Denomination == "Reformatorisch"').BRIN
    
    for node in nodes_op:
        try:
            neighbours.extend(nx.neighbors(schools_net, node))
        except: 
            print(node + "  missing")
        
    return list(set(neighbours))


def assign_vacc_status(schools_net, schools, vacc_data, schools_data_gdf,  average_uptake):
    
    vacc_stats = []
    
    vacc_data.school = [s.replace('_', '0') for s in vacc_data.school]
    
    
    for n in schools_net.nodes():
        if n in np.array(vacc_data.school): 
            vacc_stats.append([n, vacc_data.query('school == @n')['50%'].values[0]])
        else: 
            vacc_stats.append([n, average_uptake])
    
    
    vacc_dict = dict(vacc_stats)
    missing = 0
    vacc_imputed = []
    for n in schools_net.nodes():
        try:
            if schools_data_gdf.query("BRIN == @n").kind.values[0] == 'sec':
                n_nbrs = nx.neighbors(schools_net, n)
                vacc = sum([vacc_dict[n_nbr]*schools_net.get_edge_data(n, n_nbr)['pup_count'] for n_nbr in n_nbrs])/sum([schools_net.get_edge_data(n, n_nbr)['pup_count'] for n_nbr in n_nbrs])
                vacc_imputed.append([n, vacc])
            else: 
                vacc_imputed.append([n, vacc_dict[n]])
        except:
            missing += 1
            vacc_imputed.append([n, average_uptake])
                
    print(missing)
        
    vacc_imputed = pd.DataFrame(vacc_imputed, columns=['BRIN', 'vacc'])
    
    schools_data_v_gdf = pd.merge(schools_data_gdf, vacc_imputed, on = 'BRIN')
    
    return vacc_imputed, schools_data_v_gdf
    
