#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 09:14:31 2020

@author: vivekmodi
"""
from Bio import SearchIO

def identify_residues(pwd,pdbfilename,index,conf_df):
    model_id=conf_df.at[index,'Model_id']
    chain_id=conf_df.at[index,'Chain_id']
    group=conf_df.at[index,'Group']
    first_res=conf_df.at[index,'First_res']
    lys={'AGC':30,'CAMK':30,'CK1':30,'CMGC':30,'NEK':30,'RGC':29,'STE':30,'TKL':28,'TYR':34}
    glu={'AGC':49,'CAMK':47,'CK1':44,'CMGC':45,'NEK':48,'RGC':45,'STE':47,'TKL':46,'TYR':51}
    glu4={'AGC':53,'CAMK':51,'CK1':48,'CMGC':49,'NEK':52,'RGC':49,'STE':51,'TKL':50,'TYR':55}
    xdfg={'AGC':141,'CAMK':141,'CK1':140,'CMGC':149,'NEK':144,'RGC':139,'STE':140,'TKL':141,'TYR':145}
    asp={'AGC':142,'CAMK':142,'CK1':141,'CMGC':150,'NEK':145,'RGC':140,'STE':141,'TKL':142,'TYR':146}
    phe={'AGC':143,'CAMK':143,'CK1':142,'CMGC':151,'NEK':146,'RGC':141,'STE':142,'TKL':143,'TYR':147}

    hmm_result=SearchIO.read(f'{pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_{group}.hmmer.txt',format='hmmer3-text')
    for hits in hmm_result:     #extract hit from alignment in HMM output file
        for hsps in hits:
            col_num=0;hmm_index=hsps.query_start;hit_index=hsps.hit_start+first_res-1
            for hmm_res in hsps.aln[0]:
                col_num=col_num+1
                if hmm_res!='.':
                    hmm_index=hmm_index+1
                if hsps.aln[1][col_num-1]!='-':
                    hit_index=hit_index+1
                if hmm_index==lys[group] and hmm_res!='.':
                    conf_df.at[index,'Lys_num']=hit_index
                    conf_df.at[index,'Lys_restype']=list({hsps.aln[1][col_num-1]})[0]       #The residue is extracted using list otherwise it is printed in ''
                if hmm_index==glu[group] and hmm_res!='.':
                    conf_df.at[index,'Glu_num']=hit_index
                    conf_df.at[index,'Glu_restype']=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==glu4[group] and hmm_res!='.':
                    conf_df.at[index,'Glu4_num']=hit_index
                    conf_df.at[index,'Glu4_restype']=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==xdfg[group] and hmm_res!='.':
                    conf_df.at[index,'XDFG_num']=hit_index
                    conf_df.at[index,'XDFG_restype']=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==asp[group] and hmm_res!='.':
                    conf_df.at[index,'Asp_num']=hit_index
                    conf_df.at[index,'Asp_restype']=list({hsps.aln[1][col_num-1]})[0]
                if hmm_index==phe[group] and hmm_res!='.':
                    conf_df.at[index,'Phe_num']=hit_index
                    conf_df.at[index,'Phe_restype']=list({hsps.aln[1][col_num-1]})[0]
               
    
    return conf_df
