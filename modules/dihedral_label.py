#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math

dfgin={'BLAminus':(-128.64,178.67,61.15,81.21,-96.89,20.53,289.12),\
       'BLAplus':(-119.24,167.71,58.94,34.08,-89.42,-8.54,55.63),\
       'ABAminus':(-111.82,-7.64,-141.55,148.01,-127.79,23.32,296.17),\
       'BLBminus':(-134.79,175.48,60.44,65.35,-79.44,145.34,287.56),\
       'BLBplus':(-125.28,172.53,59.98,32.92,-85.51,145.28,49.01),\
       'BLBtrans':(-106.16,157.24,69.37,21.33,-61.73,134.56,215.23)}
dfginter={'BABtrans':(-80.20,128.22,-117.47,23.76,-85.16,133.21,181.42)}
dfgout={'BBAminus':(-138.56,-176.12,-144.35,103.66,-82.59,-9.03,290.59)}

def dihedral_label(index,conf_df):

    xdfg_phi=float(conf_df.at[index,'XDFG_Phi']);xdfg_psi=float(conf_df.at[index,'XDFG_Psi']);dfg_asp_phi=float(conf_df.at[index,'Asp_Phi'])
    dfg_asp_psi=float(conf_df.at[index,'Asp_Psi']);dfg_phe_phi=float(conf_df.at[index,'Phe_Phi']);dfg_phe_psi=float(conf_df.at[index,'Phe_Psi'])
    dfg_phe_chi1=float(conf_df.at[index,'Phe_Chi1'])
    
    if conf_df.at[index,'Spatial_label']=='Unassigned':
        conf_df.at[index,'Dihedral_label']='Unassigned'
        return conf_df
    if xdfg_phi==99999 or xdfg_psi==99999 or dfg_asp_phi==99999 or dfg_asp_psi==99999 or dfg_phe_phi==99999 or  dfg_phe_psi==99999 or dfg_phe_chi1==99999:
        conf_df.at[index,'Dihedral_label']='Unassigned'
        return conf_df
    
    if conf_df.at[index,'Spatial_label']=='DFGin':
        for clusters in dfgin:
            cosine_dis=(2/7)*((1-math.cos(math.radians(xdfg_phi-float(dfgin[clusters][0]))))+(1-math.cos(math.radians(xdfg_psi-float(dfgin[clusters][1]))))+\
                        (1-math.cos(math.radians(dfg_asp_phi-float(dfgin[clusters][2]))))+(1-math.cos(math.radians(dfg_asp_psi-float(dfgin[clusters][3]))))+\
                        (1-math.cos(math.radians(dfg_phe_phi-float(dfgin[clusters][4]))))+(1-math.cos(math.radians(dfg_phe_psi-float(dfgin[clusters][5]))))+\
                        (1-math.cos(math.radians(dfg_phe_chi1-float(dfgin[clusters][6])))))

            if cosine_dis<0.3:
                conf_df.at[index,'Dihedral_label']=clusters


    if conf_df.at[index,'Spatial_label']=='DFGinter':
        for clusters in dfginter:
            cosine_dis=(2/7)*((1-math.cos(math.radians(xdfg_phi-float(dfginter[clusters][0]))))+(1-math.cos(math.radians(xdfg_psi-float(dfginter[clusters][1]))))+\
                        (1-math.cos(math.radians(dfg_asp_phi-float(dfginter[clusters][2]))))+(1-math.cos(math.radians(dfg_asp_psi-float(dfginter[clusters][3]))))+\
                        (1-math.cos(math.radians(dfg_phe_phi-float(dfginter[clusters][4]))))+(1-math.cos(math.radians(dfg_phe_psi-float(dfginter[clusters][5]))))+\
                        (1-math.cos(math.radians(dfg_phe_chi1-float(dfginter[clusters][6])))))
            
            if cosine_dis<0.3:
                conf_df.at[index,'Dihedral_label']='BABtrans'

    if conf_df.at[index,'Spatial_label']=='DFGout':
        for clusters in dfgout:
            cosine_dis=(2/7)*((1-math.cos(math.radians(xdfg_phi-float(dfgout[clusters][0]))))+(1-math.cos(math.radians(xdfg_psi-float(dfgout[clusters][1]))))+\
                        (1-math.cos(math.radians(dfg_asp_phi-float(dfgout[clusters][2]))))+(1-math.cos(math.radians(dfg_asp_psi-float(dfgout[clusters][3]))))+\
                        (1-math.cos(math.radians(dfg_phe_phi-float(dfgout[clusters][4]))))+(1-math.cos(math.radians(dfg_phe_psi-float(dfgout[clusters][5]))))+\
                        (1-math.cos(math.radians(dfg_phe_chi1-float(dfgout[clusters][6])))))
            
            if cosine_dis<0.3:
                conf_df.at[index,'Dihedral_label']='BBAminus'
    
    return conf_df