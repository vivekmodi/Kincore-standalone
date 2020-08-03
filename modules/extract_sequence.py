#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 07:49:05 2020

@author: vivekmodi
"""
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def extract_seq(pwd,pdbfilename,index,conf_df):
    model_id=conf_df.at[index,'Model_id']
    chain_id=conf_df.at[index,'Chain_id']
    fhandle_outputseq=open(f"{pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta",'w')
    fhandle_outputseq.write(f">{pdbfilename[0:-4]}_{model_id}_{chain_id}\n")
    
    if '.cif' in pdbfilename.lower():
        for record in SeqIO.parse(f'{pwd}/{pdbfilename}','cif-atom'):      #Extracts sequence from ATOM records
            if str(record.annotations['model'])==model_id and record.annotations['chain']==chain_id:
                fhandle_outputseq.write(f'{record.seq}')
                conf_df.at[index,'Sequence']=record.seq
                conf_df.at[index,'First_res']=record.annotations['start']               #Extract the residue number of first residue as this information is lost in fasta file and HMM output
    if '.pdb' in pdbfilename.lower():
        for record in SeqIO.parse(f'{pwd}/{pdbfilename}','pdb-atom'):      #Extracts sequence from ATOM records
            if str(record.annotations['model'])==model_id and record.annotations['chain']==chain_id:
                fhandle_outputseq.write(f'{record.seq}')
                conf_df.at[index,'Sequence']=record.seq
                conf_df.at[index,'First_res']=record.annotations['start']               #Extract the residue number of first residue as this information is lost in fasta file and HMM output
    fhandle_outputseq.close()
    return conf_df                      
