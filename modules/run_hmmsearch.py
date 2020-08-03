#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 08:29:49 2020

@author: vivekmodi
"""
import subprocess

def run_hmmsearch(pwd,pdbfilename,index,conf_df):
    model_id=conf_df.at[index,'Model_id']
    chain_id=conf_df.at[index,'Chain_id']
    
    cmd=(f'hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_AGC.hmmer.txt {pwd}/HMMs/AGC.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_CAMK.hmmer.txt {pwd}/HMMs/CAMK.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_CK1.hmmer.txt {pwd}/HMMs/CK1.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_CMGC.hmmer.txt {pwd}/HMMs/CMGC.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_NEK.hmmer.txt {pwd}/HMMs/NEK.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_RGC.hmmer.txt {pwd}/HMMs/RGC.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_STE.hmmer.txt {pwd}/HMMs/STE.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_TKL.hmmer.txt {pwd}/HMMs/TKL.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;\
         hmmsearch -o {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}_TYR.hmmer.txt {pwd}/HMMs/TYR.hmm {pwd}/{pdbfilename[0:-4]}_{model_id}_{chain_id}.fasta;')
    
    subprocess.call(cmd,shell=True)
    return
