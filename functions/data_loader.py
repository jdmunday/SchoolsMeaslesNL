#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 09:16:56 2018

@author: LAPT0084
"""

import pandas as pd
import numpy as np
import geopandas as gpd 
import os
from shapely.geometry import Point

basepath = os.getcwd()
print(basepath)
#filepath = path.abspath(path.join(basepath, "..", "..", "fileIwantToOpen.txt"))


schools_path = basepath + r'/data'


def load_schools_data(schools_path = schools_path):
    
    pri_school_locs1 = pd.read_csv(schools_path + r'/location data/bo_jan2013.csv', encoding='latin1')
    pri_school_locs2 = pd.read_csv(schools_path + r'/location data/so_jan2013.csv', encoding='latin1')
    sec_school_locs = pd.read_csv(schools_path + r'/location data/vo_jan2013.csv', encoding='latin1')
    
    pri_schools_data1 = pri_school_locs1[['vestigingsnummer', 'pc6', 'provincie', 'leerlingen',  'denominatie', 'xcoord', 'ycoord']].rename(columns={'vestigingsnummer':'BRIN', 'denominatie':'Denomination'})
    pri_schools_data1['kind'] = ['pri']*len(pri_schools_data1)
    pri_schools_data2 = pri_school_locs2[['vestigingsnummer', 'pc6', 'provincie', 'leerlingen',  'denominatie', 'xcoord', 'ycoord']].rename(columns={'vestigingsnummer':'BRIN', 'denominatie':'Denomination'})
    pri_schools_data2['kind'] = ['pri']*len(pri_schools_data2)
    sec_schools_data = sec_school_locs[['vestigingsnummer', 'pc6', 'provincie', 'leerlingen',  'denominatie', 'xcoord', 'ycoord']].rename(columns={'vestigingsnummer':'BRIN', 'denominatie':'Denomination'})
    sec_schools_data['kind'] = ['sec']*len(sec_schools_data)
    
    schools_data = pd.concat([pri_schools_data1, pri_schools_data2, sec_schools_data])
    
    schools_data = schools_data.rename(columns = {'xcoord':'x_coord', 'ycoord':'y_coord'})
    
    geometry = [Point(xy) for xy in zip(schools_data.x_coord, schools_data.y_coord)]
    
    schools_data['geometry'] = geometry
    
    schools_data_gdf = gpd.GeoDataFrame(schools_data)
    
    return schools_data_gdf


def load_schools_2_pc4():

    schools_pc4_vo = pd.read_csv(schools_path + r'/registration_data/02.-leerlingen-per-vestiging-naar-postcode-en-leerjaar-2012.csv', delimiter=';', encoding='latin1')
    zeropad = lambda n: '%02d' % n
    BRIN = np.array(schools_pc4_vo['BRIN NUMMER']) + np.array(list(map(zeropad, schools_pc4_vo['VESTIGINGSNUMMER'])))
    schools_pc4_vo['BRIN'] = BRIN
    
    schools_pc4_vo['TOTAAL'] = schools_pc4_vo[[u'LEER- OF VERBLIJFSJAAR 1',
       u'LEER- OF VERBLIJFSJAAR 2', u'LEER- OF VERBLIJFSJAAR 3',
       u'LEER- OF VERBLIJFSJAAR 4', u'LEER- OF VERBLIJFSJAAR 5',
       u'LEER- OF VERBLIJFSJAAR 6']].sum(axis=1)
    schools_pc4_vo = schools_pc4_vo.rename(columns = {'POSTCODE LEERLING':'PC4'})
    
    schools_pc4_vo = schools_pc4_vo[['BRIN', 'PC4', 'TOTAAL']]
    
    schools_pc4_po = pd.read_csv(schools_path + r'/registration_data/03.-leerlingen-po-totaaloverzicht-2012-2013.csv', delimiter=';', encoding='latin1')
    
    
    BRIN = np.array(schools_pc4_po['BRIN_NUMMER']) + np.array(list(map(zeropad, schools_pc4_po['VESTIGINGSNUMMER'])))
    schools_pc4_po['BRIN'] = BRIN
    schools_pc4_po = schools_pc4_po.rename(columns = {'POSTCODE_LEERLING':'PC4'})
    
    schools_pc4_po = schools_pc4_po[['BRIN', 'PC4', 'TOTAAL']]
    
    schools_pc4 = pd.concat((schools_pc4_po,schools_pc4_vo))
    
    return schools_pc4


def load_case_data_pc4():
    
    meas_data = pd.read_csv(schools_path + r'/Measles_per_PC4.csv')
    meas_data = meas_data.rename(columns={'Postcode':'pc4'})
    
    return meas_data
    
    
    
    
    
    