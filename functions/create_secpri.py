#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 14:38:26 2018

@author: LAPT0084
"""


import pandas as pd
import numpy as np 
import os



basepath = os.path.dirname(__file__)
#filepath = path.abspath(path.join(basepath, "..", "..", "fileIwantToOpen.txt"))
schools_path = basepath + '/Dutch data'

def create_secpri_multi(years=range(2012,2017)):


    feeder_sets = []
    
    for year in years: 
        flow_sec_prim_file = 'Primary_to_secondary_' + str(year) + '.csv'
        secpri = pd.read_csv(schools_path + '//' +  flow_sec_prim_file, skiprows=[0])
        
        print( year)
        
        zeropad = lambda n: '%02d' % n
        
        vestigingsnummer_pri = np.array(secpri['BRIN NUMMER']) + np.array(map(zeropad, secpri['VESTIGINGSNUMMER']))
        vestigingsnummer_sec = np.array(secpri['BRIN NUMMER.1']) + np.array(map(zeropad, secpri['VESTIGINGSNUMMER.1']))
        
        
        secpri['Primary_School_URN'] = vestigingsnummer_pri
        secpri['SecondaryURN'] = vestigingsnummer_sec
        secpri['link'] = secpri['Primary_School_URN'] + secpri['SecondaryURN']
        secpri = secpri.rename(columns = {'AANTAL DOORSTROMERS':'pupil_count_'+str(year)})
        
        
        feeder_sets.append(secpri)
    
    secpri = feeder_sets[0][['link','Primary_School_URN' , 'SecondaryURN', 'pupil_count_2012']]
    
    for i in xrange(1,5):
        print( years[i])
        secpri = secpri.merge(feeder_sets[i][['link','pupil_count_' + str(years[i])]], on='link')
        
    secpri['pupil_count'] = secpri[[col for col in secpri.columns if 'pupil' in col]].sum(axis=1)
    
    return secpri
