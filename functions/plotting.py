#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 11:12:00 2018

@author: LAPT0084
"""

import networkx as nx 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec
import geopandas as gpd
import imageio
import numpy.ma as ma
import os

import data_loader as dl


import models as sar
import vaccination as iv
import importlib as imp



imp.reload(iv)



def plot_vacc_cov(schools_data_v_gdf):
    
    fig = plt.figure(figsize = [40, 30])
    ax = fig.add_subplot(111)
    
    ax.patch.set_facecolor('white')
    
    #lond_shapes.plot(color='MidnightBlue', ax=ax, edgecolor='white')
    #ad_bounds.plot(ax=ax)
    
    schools_data_v_gdf.plot(column='vacc', cmap='Blues_r', ax=ax, markersize=10, alpha=0.8)
        
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_aspect('equal')
    


def plot_schools_gen(selections, color=['Grey']):
    
    fig = plt.figure(figsize = [40, 30])
    ax = fig.add_subplot(111)
    
    ax.patch.set_facecolor('white')
    
    #lond_shapes.plot(color='MidnightBlue', ax=ax, edgecolor='white')
    #ad_bounds.plot(ax=ax)
    
    for i, sel in enumerate(selections):
    
        sel.plot(color=color[i], ax=ax, markersize=10, alpha=0.8)
        
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_aspect('equal')
    
    
    
def make_network_gif(G, nodes, s, i, r, vacc_dict, output_file = r'/Users/LAPT0084/Desktop/gif_figs'):
    fnams=[]
    pos_rg = nx.spring_layout(G, weight = 'weight', iterations=200)
    for n in xrange(len(r)):
        fig = plt.figure(figsize=[20,16])
        
        spec = GridSpec(3, 4).new_subplotspec((0,0), rowspan=3, colspan=3)
        ax = fig.add_subplot(spec)
        
        
        spec1 = GridSpec(3, 4).new_subplotspec((0,3), rowspan=1, colspan=1)
        ax1 = fig.add_subplot(spec1)
        spec2 = GridSpec(3, 4).new_subplotspec((1,3), rowspan=1, colspan=1)
        ax2 = fig.add_subplot(spec2)
        spec3 = GridSpec(3, 4).new_subplotspec((2,3), rowspan=1, colspan=1)
        ax3 = fig.add_subplot(spec3)

        
        nx.draw(G, nodelist=list(np.array(nodes)[ma.make_mask(map(int, r[n]))]), pos=pos_rg, ax=ax, node_color='DodgerBlue', edge_color='Grey', alpha=0.7, node_size=[200.*(1-vacc_dict[v]) for v in list(np.array(nodes)[ma.make_mask(map(int, r[n]))])])
        nx.draw_networkx_nodes(G, nodelist=list(np.array(nodes)[ma.make_mask(map(int, i[n]))]), pos=pos_rg, ax=ax, node_size=[200.*(1-vacc_dict[v]) for v in list(np.array(nodes)[ma.make_mask(map(int, i[n]))])], node_color='r', alpha=0.7)
        nx.draw_networkx_nodes(G, nodelist=list(np.array(nodes)[ma.make_mask(map(int, s[n]))]), pos=pos_rg, ax=ax, node_size=[200.*(1-vacc_dict[v]) for v in list(np.array(nodes)[ma.make_mask(map(int, s[n]))])], node_color='SeaGreen', alpha=0.7)
        
        ax1.plot(sum(np.array(s[:n]).T), c='SeaGreen', alpha=0.7, linewidth=2.)
        ax1.set_ylabel('Susceptible schools')
        ax2.plot(sum(np.array(i[:n]).T), c='r', alpha=0.7, linewidth=2.)
        ax2.set_ylabel('Infected schools')
        ax3.plot(sum(np.array(r[:n]).T), c='DodgerBlue', alpha=0.7, linewidth=2.)
        ax3.set_ylabel('Recovered schools')
        
        for a in [ax1, ax2, ax3]:
            a.set_xlim([0,len(s)])
            a.set_ylim([0,len(s[0])])
        
        fig.savefig(output_file + r'/fig' + str(n) + '.png')
        fnams.append(output_file + r'/fig' + str(n) + '.png')
        
def make_map_gif(schools_data, nodes, s, i, r,  output_file = r'/Users/LAPT0084/Desktop/gif_figs'):
    fnams=[]
    for n in xrange(len(r)):
        fig = plt.figure(figsize=[20,20])
        
        spec = GridSpec(3, 4).new_subplotspec((0,0), rowspan=3, colspan=3)
        ax = fig.add_subplot(spec)
        ax.patch.set_facecolor('white')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_aspect('equal')
        
        
        spec1 = GridSpec(3, 4).new_subplotspec((0,3), rowspan=1, colspan=1)
        ax1 = fig.add_subplot(spec1)
        spec2 = GridSpec(3, 4).new_subplotspec((1,3), rowspan=1, colspan=1)
        ax2 = fig.add_subplot(spec2)
        spec3 = GridSpec(3, 4).new_subplotspec((2,3), rowspan=1, colspan=1)
        ax3 = fig.add_subplot(spec3)
        
        infsch = pd.DataFrame(np.transpose([nodes, s[n], i[n], r[n]]), columns = ['BRIN', 's', 'i', 'r'])
        sdvi = pd.DataFrame(schools_data).merge(infsch, on='BRIN')
        sdvi = sdvi.convert_objects(convert_numeric=True)
        i_nodes = np.array(nodes)[ma.make_mask(map(int, i[n]))]
        r_nodes = np.array(nodes)[ma.make_mask(map(int, r[n]))]
        
        
        
        sdvi.query('BRIN not in @r_nodes').plot.scatter(x='x_coord', y='y_coord', c='LightGrey', alpha=0.6, ax=ax)
        
        if len(r_nodes) > 0:
            sdvi.query('BRIN in @r_nodes').plot.scatter(x='x_coord', y='y_coord', c='DodgerBlue', ax=ax)
        if len(i_nodes) > 0:
            sdvi.query('BRIN in @i_nodes').plot.scatter(x='x_coord', y='y_coord', c='Red', ax=ax)

        
        ax1.plot(sum(np.array(s[:n]).T), c='SeaGreen', alpha=0.7, linewidth=2.)
        ax1.set_ylabel('Susceptible schools')
        ax2.plot(sum(np.array(i[:n]).T), c='r', alpha=0.7, linewidth=2.)
        ax2.set_ylabel('Infected schools')
        ax3.plot(sum(np.array(r[:n]).T), c='DodgerBlue', alpha=0.7, linewidth=2.)
        ax3.set_ylabel('Recovered schools')
        
        for a in [ax1, ax2, ax3]:
            a.set_xlim([0,len(s)])
        
        ax1.set_ylim([len(s[0])-200., len(s[0])])
        ax2.set_ylim([0,200])
        ax3.set_ylim([0,200])
        
        fig.savefig(output_file + r'/map' + str(n) + '.png')
        fnams.append(output_file + r'/map' + str(n) + '.png')
        images = []
        
        
        
        
    for filename in fnams:
        for frames in range(10):
            images.append(imageio.imread(filename))
    imageio.mimsave(r'/Users/LAPT0084/Desktop/gif_figs/mapgif.gif', images)
    


def make_map_gif_multi(schools_data, nodes, res, samplesize=5, output_file = r'/Users/LAPT0084/Desktop/gif_figs_lon'):
    
    fnams=[]
    finals = []
    finali = []
    finalr = []
    for re in np.random.choice(range(len(res)), samplesize,  replace=False):
        
        s, i, r = res[re][1:4]
        for n in xrange(len(r)):
            fig = plt.figure(figsize=[10,10])
            
            
            spec = GridSpec(3, 4).new_subplotspec((0,0), rowspan=3, colspan=3)
            ax = fig.add_subplot(spec)
            ax.patch.set_facecolor('white')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
    
            ax.set_aspect('equal')
            
            
            spec1 = GridSpec(3, 4).new_subplotspec((0,3), rowspan=1, colspan=1)
            ax1 = fig.add_subplot(spec1)
            spec2 = GridSpec(3, 4).new_subplotspec((1,3), rowspan=1, colspan=1)
            ax2 = fig.add_subplot(spec2)
            spec3 = GridSpec(3, 4).new_subplotspec((2,3), rowspan=1, colspan=1)
            ax3 = fig.add_subplot(spec3)
            
            if len(finals) > 0:
                for final in finals:
                    ax1.plot(final, c='Purple', alpha=0.2, linewidth=2.)
                
                for final in finalr:
                    ax3.plot(final, c='Purple', alpha=0.2, linewidth=2.)
                    
                for final in finali:
                    ax2.plot(final, c='Purple', alpha=0.2, linewidth=2.)
            
            infsch = pd.DataFrame(np.transpose([map(str,nodes), s[n], i[n], r[n]]), columns = ['BRIN', 's', 'i', 'r'])
            sdvi = pd.DataFrame(schools_data).merge(infsch, on='BRIN')
            sdvi = sdvi.convert_objects(convert_numeric=True)
            i_nodes = np.array(nodes)[ma.make_mask(map(int, i[n]))]
            r_nodes = np.array(nodes)[ma.make_mask(map(int, r[n]))]
            
            
            
            sdvi.plot.scatter(x='x_coord', y='y_coord', c='LightGrey', alpha=0.6, ax=ax)
            
            if len(r_nodes) > 0:
                try:
                    #return sdvi, r_nodes
                    sdvi.query('BRIN in @r_nodes').plot.scatter(x='x_coord', y='y_coord', c='DodgerBlue', ax=ax)
                except:
                    foo = 0
            if len(i_nodes) > 0:
                try:
                    sdvi.query('BRIN in @i_nodes').plot.scatter(x='x_coord', y='y_coord', c='Red', ax=ax)
                except:
                    foo = 0
    
            
            ax1.plot(sum(np.array(s[:n]).T), c='SeaGreen', alpha=0.7, linewidth=2.)
            ax1.set_ylabel('Susceptible schools')
            ax2.plot(sum(np.array(i[:n]).T), c='r', alpha=0.7, linewidth=2.)
            ax2.set_ylabel('Infected schools')
            ax3.plot(sum(np.array(r[:n]).T), c='DodgerBlue', alpha=0.7, linewidth=2.)
            ax3.set_ylabel('Recovered schools')
            
            
            #for a in [ax1, ax2, ax3]:
            #    a.set_xlim([0,20])
            
            #ax1.set_ylim([0, len(s[0])])
            #ax2.set_ylim([0,len(s[0])])
            #ax3.set_ylim([0,len(s[0])])
            
            ax.set_xlabel('')
            ax.set_ylabel('')
            
            fig.savefig(output_file + r'/map' + str(re) + '_' + str(n) + '.png')
            fnams.append(output_file + r'/map' + str(re) + '_' + str(n) + '.png')
            images = []
            plt.close()
        finals.append(list(sum(np.array(s).T)))
        finali.append(list(sum(np.array(i).T)))
        finalr.append(list(sum(np.array(r).T)))
        
    for filename in fnams:
        for frames in range(3):
            images.append(imageio.imread(filename))
    imageio.mimsave(r'/Users/LAPT0084/Desktop/gif_figs/mapgif.gif', images)
    
    return finals

def plot_risk(schools_data_v_gdf, res):
    rs = np.array([r[-1] for r in np.array(res)[:,-2]])
    infsch = pd.DataFrame(np.transpose([res[0][0], list(sum(rs))]), columns = ['BRIN', 'N_inf'])
    fig = plt.figure(figsize=[40,20])
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    schools_risk = pd.DataFrame(schools_data_v_gdf.merge(infsch, on='BRIN'))
    schools_risk = schools_risk.convert_objects(convert_numeric=True)
    schools_risk['risk'] = 1.*schools_risk.N_inf / len(rs)
    #schools_data_v_gdf.merge(infsch, on='BRIN').sort_values('N_inf').plot('N_inf', cmap='Reds', ax=ax1, markersize=10, alpha=0.8)
    #schools_data_v_gdf.merge(infsch, on='BRIN').sort_values('vacc', ascending=False).plot('vacc', cmap='Blues_r', ax=ax2, markersize=10, alpha=0.8)
    schools_risk.sort_values('risk').plot.scatter(x='x_coord', y='y_coord', 
                            c=schools_risk.sort_values('risk').N_inf, 
                            cmap='Reds', ax=ax1, alpha=0.8, colorbar=False, 
                            s=100)
    schools_risk.sort_values('vacc', ascending = False).plot.scatter(x='x_coord', 
                            y='y_coord', c=schools_risk.sort_values('vacc', ascending = False).vacc, 
                            cmap='Blues_r', ax=ax2, alpha=0.8, colorbar=False, 
                            s=100)
    
    for ax in [ax1, ax2]:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
    
        ax.set_aspect('equal')
    
    fig.tight_layout()
    
    
    return schools_risk

def plot_risk_lond(schools_data_v_gdf, res):
    rs = np.array([r[-1] for r in np.array(res)[:,-2]])
    infsch = pd.DataFrame(np.transpose([map(str,res[0][0]), list(sum(rs))]), columns = ['URN', 'N_inf'])
    fig = plt.figure(figsize=[20,20])
    ax1=fig.add_subplot(111)
    
    schools_risk = pd.DataFrame(schools_data_v_gdf.merge(infsch, on='URN'))
    schools_risk = schools_risk.convert_objects(convert_numeric=True)
    schools_risk['risk'] = 1.*schools_risk.N_inf / len(rs)

    schools_risk['x_coord'] = [g.x for g in schools_risk.geometry]
    schools_risk['y_coord'] = [g.y for g in schools_risk.geometry]
    schools_risk.sort_values('risk').plot.scatter(x='x_coord', y='y_coord', 
                            c=schools_risk.sort_values('risk').N_inf, 
                            cmap='Reds', ax=ax1, alpha=0.8, colorbar=False, 
                            s=100)
    
    
    for ax in [ax1]:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
    
        ax.set_aspect('equal')
    
    fig.tight_layout()
    
    return schools_risk 

            
def plot_risk_map(res, vacc_dict, n=20.):
    
    rs = np.array([r[-1] for r in np.array(res)[:,-2]])
    
    infsch = pd.DataFrame(np.transpose([res[0][0], list(sum(rs))]), columns = ['BRIN', 'N_inf'])
    infsch = infsch.convert_objects(convert_numeric=True)
    
    schools_pc4_vo = pd.read_csv(r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/02.-leerlingen-per-vestiging-naar-postcode-en-leerjaar-2017-2018.csv')
    zeropad = lambda n: '%02d' % n
    BRIN = np.array(schools_pc4_vo['BRIN NUMMER']) + np.array(map(zeropad, schools_pc4_vo['VESTIGINGSNUMMER']))
    schools_pc4_vo['BRIN'] = BRIN
    
    schools_pc4_vo['TOTAAL'] = schools_pc4_vo[[u'LEER- OF VERBLIJFSJAAR 1',
       u'LEER- OF VERBLIJFSJAAR 2', u'LEER- OF VERBLIJFSJAAR 3',
       u'LEER- OF VERBLIJFSJAAR 4', u'LEER- OF VERBLIJFSJAAR 5',
       u'LEER- OF VERBLIJFSJAAR 6']].sum(axis=1)
    schools_pc4_vo = schools_pc4_vo.rename(columns = {'POSTCODE LEERLING':'PC4'})
    
    schools_pc4_vo = schools_pc4_vo[['BRIN', 'PC4', 'TOTAAL']]
    
    PO_path = r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/03.-leerlingen-po-per-gemeente-postcode-leerling,-leeftijd-2017-2018'
    
    schools_pc4_po = pd.DataFrame()
    
    for pofile in os.listdir(PO_path):
        
        schools_pc4_po = schools_pc4_po.append(pd.read_csv(PO_path + '//' + pofile, delimiter=';'))
    
    
    BRIN = np.array(schools_pc4_po['BRIN_NUMMER']) + np.array(map(zeropad, schools_pc4_po['VESTIGINGSNUMMER']))
    schools_pc4_po['BRIN'] = BRIN
    schools_pc4_po = schools_pc4_po.rename(columns = {'POSTCODE_LEERLING':'PC4'})
    
    schools_pc4_po = schools_pc4_po[['BRIN', 'PC4', 'TOTAAL']]
    
    schools_pc4 = schools_pc4_po.append(schools_pc4_vo)
    
    #return schools_pc4

    
    fs = []
    for B in infsch.BRIN:
        fs.append([B, (1.-vacc_dict[B])*sar.op_fs(1.-vacc_dict[B], 15.)])
    
    fsdf = pd.DataFrame(fs, columns=['BRIN', 'fs'])
    schools_pc4_cases = schools_pc4.merge(infsch, on='BRIN')
    schools_pc4_cases = schools_pc4_cases.convert_objects(convert_numeric=True)
    #return infsch, schools_pc4_cases 
    
    schools_pc4_cases = schools_pc4_cases.merge(fsdf, on='BRIN')
    
    
    
    schools_pc4_cases['cases'] = schools_pc4_cases.N_inf * schools_pc4_cases.TOTAAL * schools_pc4_cases.fs
    
    PC_Cases = schools_pc4_cases.groupby('PC4').sum()[['cases']].reset_index()
    #PC_Cases_mean = schools_pc4_cases.groupby('PC4').mean()[['cases']].reset_index().rename(columns={'cases':'cmean'})
    #PC_Cases_std  = schools_pc4_cases.groupby('PC4').std()[['cases']].reset_index().rename(columns={'cases':'cstd'})
    
    #PC_Cases = PC_Cases.rename(columns = {'POSTCODE LEERLING':'PC4'})
    
    vdf = pd.DataFrame.from_dict(vacc_dict, orient='Index').reset_index()
    vdf.columns = ['BRIN', 'vacc']
    
    schools_pc4_vacc = schools_pc4.merge(vdf, on='BRIN')
    schools_pc4_vacc['vaccpc'] = schools_pc4_vacc.vacc * schools_pc4_vacc.TOTAAL
    schools_pc4_vacc_sums = schools_pc4_vacc.groupby('PC4').sum()
    schools_pc4_vacc_sums['overallvacc'] = schools_pc4_vacc_sums.vaccpc / schools_pc4_vacc_sums.TOTAAL
    
    
    
    PC_Vacc = schools_pc4_vacc_sums.reset_index()[['PC4', 'overallvacc']]
    #PC_Vacc = PC_Vacc.rename(columns={'POSTCODE LEERLING':'PC4'})
    
    #return schools_pc4, PC_Vacc
    
    pc4shapes = gpd.read_file(r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/esri-openpostcodevlakkenpc4/ESRI-PC4-2015R1.shp')
    
    pc4shapes = pc4shapes.merge(PC_Cases, on='PC4')
    #pc4shapes = pc4shapes.merge(PC_Cases_mean, on='PC4')
    #pc4shapes = pc4shapes.merge(PC_Cases_std, on='PC4')
    pc4shapes = pc4shapes.merge(PC_Vacc, on='PC4')
    
    
    pc4shapes['log_cases'] = np.log(abs(pc4shapes.cases))
    pc4shapes['sus'] = 1. - pc4shapes['overallvacc']
    
    #return schools_pc4, pc4shapes
    
    
    fig = plt.figure(figsize=[20, 8])
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    
    pc4shapes.plot('sus', linewidth=0.0, figsize=[20,20], cmap='YlOrRd', ax=ax1)
    
    
    
    pc4shapes.plot(linewidth=0.0, figsize=[20,20], color='LightGrey', ax=ax2)
    pc4shapes[pc4shapes.cases/n>0.0001].plot('log_cases', linewidth=0.0, figsize=[20,20], cmap='Blues', ax=ax2)
    
    ax3.scatter(pc4shapes[pc4shapes.cases/n > 0.00001].overallvacc, pc4shapes[pc4shapes.cases/n > 0.00001].cases/n, edgecolor='k', facecolor='none', s=8)
    ax3.set_yscale('log')
    ax3.grid()
    ax3.set_ylim([0.005, 200])
    
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlabel('')
    ax1.set_ylabel('')

    ax1.set_aspect('equal')
    
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)

    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xlabel('')
    ax2.set_ylabel('')

    ax2.set_aspect('equal')
    
    return schools_pc4_cases, pc4shapes



    
def calculate_pc4_risk(res, vacc_dict):
    
    no_runs = len(res)
    
    rs = np.array([r[-1] for r in np.array(res)[:,-2]])
    all_res = pd.DataFrame(np.transpose([res[0][0]] + map(list, list(rs))), columns = ['BRIN'] + ['run ' + str(run) for run in range(1,no_runs + 1)])
    
    
    schools_pc4_vo = pd.read_csv(r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/02.-leerlingen-per-vestiging-naar-postcode-en-leerjaar-2017-2018.csv')
    zeropad = lambda n: '%02d' % n
    BRIN = np.array(schools_pc4_vo['BRIN NUMMER']) + np.array(map(zeropad, schools_pc4_vo['VESTIGINGSNUMMER']))
    schools_pc4_vo['BRIN'] = BRIN
    
    schools_pc4_vo['TOTAAL'] = schools_pc4_vo[[u'LEER- OF VERBLIJFSJAAR 1',
       u'LEER- OF VERBLIJFSJAAR 2', u'LEER- OF VERBLIJFSJAAR 3',
       u'LEER- OF VERBLIJFSJAAR 4', u'LEER- OF VERBLIJFSJAAR 5',
       u'LEER- OF VERBLIJFSJAAR 6']].sum(axis=1)
    schools_pc4_vo = schools_pc4_vo.rename(columns = {'POSTCODE LEERLING':'PC4'})
    
    schools_pc4_vo = schools_pc4_vo[['BRIN', 'PC4', 'TOTAAL']]
    
    PO_path = r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/03.-leerlingen-po-per-gemeente-postcode-leerling,-leeftijd-2017-2018'
    
    schools_pc4_po = pd.DataFrame()
    
    for pofile in os.listdir(PO_path):
        
        schools_pc4_po = schools_pc4_po.append(pd.read_csv(PO_path + '//' + pofile, delimiter=';'))
    
    
    BRIN = np.array(schools_pc4_po['BRIN_NUMMER']) + np.array(map(zeropad, schools_pc4_po['VESTIGINGSNUMMER']))
    schools_pc4_po['BRIN'] = BRIN
    schools_pc4_po = schools_pc4_po.rename(columns = {'POSTCODE_LEERLING':'PC4'})
    
    schools_pc4_po = schools_pc4_po[['BRIN', 'PC4', 'TOTAAL']]
    
    schools_pc4 = schools_pc4_po.append(schools_pc4_vo)
    
    
    fs = []
    for B in all_res.BRIN:
        fs.append([B, max((1.-vacc_dict[B])*sar.op_fs(1.-vacc_dict[B], 15.),0.)])
    
    fsdf = pd.DataFrame(fs, columns=['BRIN', 'fs'])
    schools_pc4_cases = schools_pc4.merge(all_res, on='BRIN')
    schools_pc4_cases = schools_pc4_cases.convert_objects(convert_numeric=True)
    #return infsch, schools_pc4_cases 
    
    schools_pc4_cases = schools_pc4_cases.merge(fsdf, on='BRIN')
    schools_pc4_cases = schools_pc4_cases.convert_objects(convert_numeric=True)
    
    PC_Vacc = iv.vacc_pc4_dict(res[0][0], vacc_dict, by='pc4')
    
    PC_Vacc = PC_Vacc.convert_objects(convert_numeric=True)
    school_Vacc = pd.DataFrame(np.transpose([vacc_dict.keys(), vacc_dict.values()]), columns=['BRIN', 'vacc'])
    
    schools_pc4_cases = schools_pc4_cases.merge(PC_Vacc, on='PC4')
    
    schools_pc4_cases = schools_pc4_cases.merge(school_Vacc, on='BRIN')
    #return schools_pc4_cases
    
    for run in range(1, no_runs + 1):
    
        schools_pc4_cases['cases_' + str(run)] = schools_pc4_cases['run ' + str(run)] * schools_pc4_cases.TOTAAL * schools_pc4_cases.fs
    
    PC_Cases = schools_pc4_cases.groupby('PC4').sum()[['cases_' + str(run) for run in range(1, no_runs + 1)]].reset_index()
    
    return PC_Cases
#    
#    infsch = pd.DataFrame(np.transpose([res[0][0] + rs)]), columns = ['BRIN', 'N_inf'])
#    infsch = infsch.convert_objects(convert_numeric=True)
    
    
def FS_dists(res_s, vacc_dicts):
    FS = []
    for i, res in enumerate(res_s):    
        PC_cases = calculate_pc4_risk(res, vacc_dicts[i])
        FS.append(list(PC_cases.sum().drop('PC4')))
    
    plt.boxplot(FS)
    
    return FS


def roc_map(PC_cases):
    
    meas_data = dl.load_case_data_pc4()
    meas_data = meas_data.rename(columns = {'pc4':'PC4'})
    
    mean_cases = PC_cases.set_index('PC4').mean(axis=1).reset_index()
    
    mean_cases = mean_cases.merge(meas_data, on='PC4', how='outer')
    mean_cases = mean_cases.rename(columns={0:'mod'})
    mean_cases = mean_cases.fillna(0)
    
    
    pc4shapes = gpd.read_file(os.getcwd() + r'/data/Location data/openpc4nl2015landonly/PC4_Nederland_2015.shp')   
    pc4shapes = pc4shapes.merge(mean_cases, on='PC4', how='outer')
    
    pc4shapes['cim'] = pc4shapes['mod'] > 0
    pc4shapes['cir'] = pc4shapes['Count'] > 0
    
    #pc4shapes = pc4shapes.dropna(subset='Shape_Area')
    
    return pc4shapes
    
#    fig = plt.figure(figsize=[7, 10])
#    ax = fig.add_subplot(111)
#    pc4shapes.plot(color='grey', linewidth=0., ax=ax)
#    pc4shapes.query('cir == True and cim == True').plot('mod', cmap =  'Greens', linewidth=0., ax=ax, vmin=0, vmax=pc4shapes.Count.max(), alpha=1.)
#    pc4shapes.query('cir == False and cim == True').plot('mod', cmap ='Reds', linewidth=0., ax=ax, vmin=0, vmax=pc4shapes.Count.max(), alpha=1.)
#    #pc4shapes.query('cir == True and cim == False').plot(color='b', linewidth=0., ax=ax)
#    ax.set_xticklabels([])
#    ax.set_yticklabels([])
#    ax.set_aspect('equal')

    
    
def roc_values(PC_cases, f, runs=20):
    
    meas_data = dl.load_case_data_pc4()
    meas_data = meas_data.rename(columns = {'pc4':'PC4'})
    
    mean_cases = PC_cases
    
    mean_cases = mean_cases.merge(meas_data, on='PC4', how='outer')
    
    mean_cases = mean_cases.fillna(0)
    
    rocs = []
    for mod in ['cases_{}'.format(n) for n in range(1,runs+1)]:
    
        rocs.append( [sum(mean_cases.query('{} > @f'.format(mod)).Count > 0) / (1. * sum(mean_cases.Count > 0)), 
         sum(mean_cases.query('{} <= @f'.format(mod)).Count == 0) / (1. * sum(mean_cases.Count == 0)),
         sum(np.array(mean_cases.query('Count > 0')[mod] > 0) / (1. * sum(mean_cases[mod] > 0))),
         sum(mean_cases.query('{} > @f'.format(mod)).Count == 0) / (1. * sum(mean_cases[mod] == 0)), 
         sum(mean_cases.query('{} <= @f'.format(mod)).Count > 0) / (1. * sum(mean_cases.Count > 0))])
    
    return rocs
    
def roc_values_weighted(PC_cases, runs=20):
    
    meas_data = dl.load_case_data_pc4()
    meas_data = meas_data.rename(columns = {'pc4':'PC4'})
    
    mean_cases = PC_cases
    
    mean_cases = mean_cases.merge(meas_data, on='PC4', how='outer')
    
    mean_cases = mean_cases.fillna(0)
    
    
    
    rocs = []
    for mod in ['cases_{}'.format(n) for n in range(1,runs+1)]:
    
        rocs.append( [sum(np.array(mean_cases.query('{} > 0'.format(mod)).Count > 0)*mean_cases.query('{} > 0'.format(mod)).Count/mean_cases.Count.sum()), 
         sum(mean_cases.query('{} <= 0'.format(mod)).Count == 0) / (1. * sum(mean_cases.Count == 0)), 
         sum(np.array(mean_cases.query('Count > 0')[mod] > 0)*mean_cases.query('Count > 0')[mod]/mean_cases[mod].sum()),
         sum(np.array(mean_cases.query('{} > 0'.format(mod)).Count == 0) * mean_cases.query('{} > 0'.format(mod))[mod]/(mean_cases[mod].sum())), 
         sum(mean_cases.query('{} <= 0'.format(mod)).Count > 0) / (1. * sum(mean_cases.Count > 0))])
    
    return rocs



def schools_pc4():
        
    

    schools_pc4_vo = pd.read_csv(r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/02.-leerlingen-per-vestiging-naar-postcode-en-leerjaar-2017-2018.csv')
    zeropad = lambda n: '%02d' % n
    BRIN = np.array(schools_pc4_vo['BRIN NUMMER']) + np.array(map(zeropad, schools_pc4_vo['VESTIGINGSNUMMER']))
    schools_pc4_vo['BRIN'] = BRIN
    
    schools_pc4_vo['TOTAAL'] = schools_pc4_vo[[u'LEER- OF VERBLIJFSJAAR 1',
       u'LEER- OF VERBLIJFSJAAR 2', u'LEER- OF VERBLIJFSJAAR 3',
       u'LEER- OF VERBLIJFSJAAR 4', u'LEER- OF VERBLIJFSJAAR 5',
       u'LEER- OF VERBLIJFSJAAR 6']].sum(axis=1)
    schools_pc4_vo = schools_pc4_vo.rename(columns = {'POSTCODE LEERLING':'PC4'})
    
    schools_pc4_vo = schools_pc4_vo[['BRIN', 'PC4', 'TOTAAL']]
    
    PO_path = r'/Users/LAPT0084/Documents/PhD/Objective 4 -  Integration through schools/Dutch data/03.-leerlingen-po-per-gemeente-postcode-leerling,-leeftijd-2017-2018'
    
    schools_pc4_po = pd.DataFrame()
    
    for pofile in os.listdir(PO_path):
        
        schools_pc4_po = schools_pc4_po.append(pd.read_csv(PO_path + '//' + pofile, delimiter=';'))
    
    
    BRIN = np.array(schools_pc4_po['BRIN_NUMMER']) + np.array(map(zeropad, schools_pc4_po['VESTIGINGSNUMMER']))
    schools_pc4_po['BRIN'] = BRIN
    schools_pc4_po = schools_pc4_po.rename(columns = {'POSTCODE_LEERLING':'PC4'})
    
    schools_pc4_po = schools_pc4_po[['BRIN', 'PC4', 'TOTAAL']]
    
    schools_pc4 = schools_pc4_po.append(schools_pc4_vo)
    
    return schools_pc4

def plot_rocmap(pc4shapes, vmax = np.nan, legendon=False):
    if vmax == np.nan: 
        vmax = pc4shapes['mod'].max()

    
    fig = plt.figure(figsize=[7, 10])
    ax = fig.add_subplot(111)
    
    pc4shapes['logmod'] = np.log(pc4shapes['mod'])
    pc4shapes['logcount'] = np.log(pc4shapes['Count'])
    
    pc4shapes.plot(color='grey', linewidth=0., ax=ax)
    if legendon: 
        pc4shapes.query('cir == True and cim == True').plot('mod', cmap =  'Greens', legend=legendon, linewidth=0., ax=ax, vmin=0, vmax=vmax, alpha=1., legend_kwds={"location": "left", "shrink":0.3, 'label':'Mean infections in PC4s \nwhere cases reported'})
        pc4shapes.query('cir == False and cim == True').plot('mod', cmap ='Reds', legend=legendon, linewidth=0., ax=ax, vmin=0, vmax=vmax, alpha=1., legend_kwds={"location": "left", "shrink":0.3, 'label':'where cases not reported'})
    else:
        pc4shapes.query('cir == True and cim == True').plot('mod', cmap =  'Greens',  linewidth=0., ax=ax, vmin=0, vmax=vmax, alpha=1.)
        pc4shapes.query('cir == False and cim == True').plot('mod', cmap ='Reds', linewidth=0., ax=ax, vmin=0, vmax=vmax, alpha=1.)
    #pc4shapes.query('cir == True and cim == False').plot(color='b', linewidth=0., ax=ax)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_aspect('equal')
    
    
    