#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 12:27:13 2020

@author: vivekmodi
"""

def assign_default_values(index,conf_df):
    conf_df.at[index,'Sequence']='X'
    conf_df.at[index,'First_res']=99999
    conf_df.at[index,'Group']='None'
    conf_df.at[index,'Lys_num']=99999
    conf_df.at[index,'Lys_restype']='X'
    conf_df.at[index,'Glu_num']=99999
    conf_df.at[index,'Glu_restype']='X'
    conf_df.at[index,'Glu4_num']=99999
    conf_df.at[index,'Glu4_restype']='X'
    conf_df.at[index,'XDFG_num']=99999
    conf_df.at[index,'XDFG_restype']='X'
    conf_df.at[index,'Asp_num']=99999
    conf_df.at[index,'Asp_restype']='X'
    conf_df.at[index,'Phe_num']=99999
    conf_df.at[index,'Phe_restype']='X'
    conf_df.at[index,'Glu4-Phe-dis']=99999
    conf_df.at[index,'Lys-Phe-dis']=99999
    conf_df.at[index,'Lys-Glu-dis']=99999
    conf_df.at[index,'XDFG_Phi']=99999
    conf_df.at[index,'XDFG_Psi']=99999
    conf_df.at[index,'Asp_Phi']=99999
    conf_df.at[index,'Asp_Psi']=99999
    conf_df.at[index,'Phe_Phi']=99999
    conf_df.at[index,'Phe_Psi']=99999
    conf_df.at[index,'Phe_Chi1']=99999
    conf_df.at[index,'Chelix']='Unassigned'
    conf_df.at[index,'Spatial_label']='Unassigned'
    conf_df.at[index,'Dihedral_label']='Unassigned'
    
    return conf_df
    