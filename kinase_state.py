#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 31 17:35:32 2020

@author: vivekmodi
"""
 
import pandas as pd, sys, os
from Bio import PDB
from modules.assign_default_values import assign_default_values
from modules.extract_sequence import extract_seq
from modules.run_hmmsearch import run_hmmsearch
from modules.identify_group import identify_group
from modules.identify_residues import identify_residues
from modules.compute_distance import compute_distance
from modules.compute_dihedrals import compute_dihedrals
from modules.chelix import chelix_conformation
from modules.spatial_label import spatial_label
from modules.dihedral_label import dihedral_label

def identify_state(pwd,pdbfilename):
    conf_df=pd.DataFrame()
    chain_list=list()
    
    if '.cif' in pdbfilename.lower():
        parser=PDB.MMCIFParser(QUIET=True)
    if '.pdb' in pdbfilename.lower():
        parser=PDB.PDBParser(QUIET=True)
        
    structure=parser.get_structure(pdbfilename, f'{pdbfilename}')
    index=-1
    for model in structure:
        for chain in model:
            if len(chain.get_list())<=30:
                continue
            chain_list.append(chain.id)
            index+=1
            
            for residue in chain:
                if residue.id[2]!=' ':
                    print('Please enter a structure file without insertion codes.')
                    break    #How to exit program here?
            conf_df.at[index,'Model_id']=str(model.id)
            conf_df.at[index,'Chain_id']=chain.id
            conf_df=assign_default_values(index,conf_df)
            conf_df=extract_seq(pdbfilename,index,conf_df)
            run_hmmsearch(pwd,pdbfilename,index,conf_df)
            conf_df=identify_group(pdbfilename,index,conf_df)
            
            if conf_df.at[index,'Group']=='None':
                print('Probably not a protein kinase.\n')
            
            else:
                conf_df=identify_residues(pdbfilename,index,conf_df)
                if conf_df.at[index,'Lys_restype']=='X':
                    print(f"{pdbfilename}\t{conf_df.at[index,'Model_id']}\t{conf_df.at[index,'Chain_id']}\tB3-Lys residue missing.\n")
                elif conf_df.at[index,'Glu_restype']=='X':
                    print(f"{pdbfilename}\t{conf_df.at[index,'Model_id']}\t{conf_df.at[index,'Chain_id']}\tChelix-Glu residue missing.\n")
                elif conf_df.at[index,'Phe_restype']=='X':
                    print(f"{pdbfilename}\t{conf_df.at[index,'Model_id']}\t{conf_df.at[index,'Chain_id']}\tDFG-Phe residue missing.\n")
                elif conf_df.at[index,'XDFG_restype']=='X':
                    print(f"{pdbfilename}\t{conf_df.at[index,'Model_id']}\t{conf_df.at[index,'Chain_id']}\tXDFG residue missing.\n")
                elif conf_df.at[index,'Asp_restype']=='X':
                    print(f"{pdbfilename}\t{conf_df.at[index,'Model_id']}\t{conf_df.at[index,'Chain_id']}\tDFG-Asp residue missing.\n")
                else:
                    
                    conf_df=compute_distance(pdbfilename,index,conf_df)
                    conf_df=compute_dihedrals(pdbfilename,index,conf_df)
                    conf_df=chelix_conformation(index,conf_df)
                    conf_df=spatial_label(index,conf_df)
                    conf_df=dihedral_label(index,conf_df)
            
                    print(f"{pdbfilename}\t{conf_df.at[index,'Group']}\t{conf_df.at[index,'Lys_num']}\t{conf_df.at[index,'Glu_num']}\t{conf_df.at[index,'Phe_num']}\t{conf_df.at[index,'Model_id']}\t{conf_df.at[index,'Chain_id']}\t{conf_df.at[index,'Spatial_label']}\t{conf_df.at[index,'Dihedral_label']}\n")
                
    

if __name__ == '__main__':
    #pwd=os.getcwd()
    #print(pwd)
    hmm_loc=os.path.dirname(os.path.realpath(__file__))+'/HMMs'    #gets the original location of the file
    pdbfilename=sys.argv[1]
    identify_state(hmm_loc,pdbfilename)
