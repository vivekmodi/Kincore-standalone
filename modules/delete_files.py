#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:57:13 2020

@author: vivekmodi
"""
import os, subprocess

def delete_files(pdbfilename,model_id,chain_id):
    if os.path.isfile(f'{pdbfilename[0:-4]}_{model_id}_{chain_id}_AGC.hmmer.txt'):
        cmd=f'rm {pdbfilename[0:-4]}*AGC.hmmer.txt;rm {pdbfilename[0:-4]}*CAMK.hmmer.txt;rm {pdbfilename[0:-4]}*CK1.hmmer.txt;rm {pdbfilename[0:-4]}*CMGC.hmmer.txt;\
        rm {pdbfilename[0:-4]}*NEK.hmmer.txt;rm {pdbfilename[0:-4]}*RGC.hmmer.txt;rm {pdbfilename[0:-4]}*STE.hmmer.txt;rm {pdbfilename[0:-4]}*TKL.hmmer.txt;\
        rm {pdbfilename[0:-4]}*TYR.hmmer.txt;rm {pdbfilename[0:-4]}*HASP.hmmer.txt;rm {pdbfilename[0:-4]}*WNK.hmmer.txt;rm {pdbfilename[0:-4]}*BUB.hmmer.txt;\
        rm {pdbfilename[0:-4]}*ULK.hmmer.txt;rm {pdbfilename[0:-4]}*_*.fasta;'
        subprocess.call(cmd,shell=True)
