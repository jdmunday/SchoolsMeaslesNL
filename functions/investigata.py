#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 17:08:53 2018

@author: LAPT0084
"""

import numpy as np

import network_analysis as nal
import infer_vacc as iv
import Build_Network_alinks as bna
import schools_at_risk as sar

import make_network as mn

import data_loader as dl

import matplotlib.pyplot as plt


schools_data = dl.load_schools_data('/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data')

schools_pc4  = iv.schools_2_pc4()

data_net = mn.build_network_alinks()

vacc_dict_samples = sar.vacc_sample_dict()


sdat_BRIN = schools_data.BRIN.unique()
spc4_BRIN = schools_pc4.BRIN.unique()
dnet_BRIN = data_net.nodes()
vacc_BRIN = vacc_dict_samples['sample 1'].keys()


BRIN_in_all = np.intersect1d(dnet_BRIN, np.intersect1d(vacc_BRIN, np.intersect1d(spc4_BRIN, sdat_BRIN)))

BRIN_only_in_net = np.setdiff1d(dnet_BRIN, np.concatenate([vacc_BRIN, spc4_BRIN, sdat_BRIN]))

BRIN_only_in_vac = np.setdiff1d(vacc_BRIN, np.concatenate([dnet_BRIN, spc4_BRIN, sdat_BRIN]))

BRIN_only_in_dat = np.setdiff1d(sdat_BRIN, np.concatenate([dnet_BRIN, spc4_BRIN, vacc_BRIN]))

BRIN_only_in_pc4 = np.setdiff1d(spc4_BRIN, np.concatenate([dnet_BRIN, sdat_BRIN, vacc_BRIN]))

brinsets = [sdat_BRIN, spc4_BRIN, vacc_BRIN, dnet_BRIN]

header = ['school data', 'school pc4', 'vaccination', 'network']

exceptions = []
for b in brinsets:
    mid_section = []
    for s in brinsets:
        mid_section.append(len(list(np.setdiff1d(b, s))))
    exceptions.append(mid_section)
    
    
intersections = []
for b in brinsets:
    mid_section = []
    for s in brinsets:
        mid_section.append(len(list(np.intersect1d(b, s))))
    intersections.append(mid_section)
 
    
def write_table(table_data, header, table_name):
    f = open('/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/' + table_name + '.csv', 'w')
    f.write(',' + ','.join(header))            
    for i, ln in enumerate(table_data):
        f.write('\n'+ header[i] + ',' + ','.join(map(str, ln)))
    f.close()
    

write_table(intersections, header, 'intersections')
write_table(exceptions, header, 'exceptions')

plt.imshow(intersections)
plt.imshow(exceptions)
