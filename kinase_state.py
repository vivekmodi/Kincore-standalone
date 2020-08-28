#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July 31 17:35:32 2020

@author: vivekmodi
"""
 
import pandas as pd, sys, os, gzip, argparse
from Bio import PDB
from modules.assign_default_values import assign_default_values
from modules.extract_sequence import extract_seq
from modules.run_hmmsearch import run_hmmsearch
from modules.identify_restypes import identify_restypes
from modules.identify_group import identify_group
from modules.identify_residues import identify_residues
from modules.compute_distance import compute_distance
from modules.compute_dihedrals import compute_dihedrals
from modules.chelix import chelix_conformation
from modules.spatial_label import spatial_label
from modules.dihedral_label import dihedral_label
from modules.delete_files import delete_files

def identify_state(pwd,pdbfilename,align,user_model,user_chain,user_lys,user_glu,user_phe):
    conf_df=pd.DataFrame()
    chain_list=list()
    
    if '.gz' in pdbfilename.lower():
        handle=gzip.open(pdbfilename,'rt')
    else:
        handle=open(pdbfilename,'r')
        
    if '.cif' in pdbfilename.lower():
            parser=PDB.MMCIFParser(QUIET=True)
    if '.pdb' in pdbfilename.lower():
            parser=PDB.PDBParser(QUIET=True)
            
    structure=parser.get_structure(pdbfilename, handle)
        
    if align.upper()=='FALSE':
        index=int(user_model)
        conf_df.at[index,'Model_id']=str(user_model)
        conf_df.at[index,'Chain_id']=str(user_chain)
        conf_df=assign_default_values(index,conf_df)
        conf_df.at[index,'Lys_num']=int(user_lys)
        conf_df.at[index,'Glu_num']=int(user_glu)
        conf_df.at[index,'Glu4_num']=int(user_glu)+4
        conf_df.at[index,'Phe_num']=int(user_phe)
        conf_df.at[index,'XDFG_num']=int(user_phe)-2
        conf_df.at[index,'Asp_num']=int(user_phe)-1
        conf_df=identify_restypes(pdbfilename,conf_df,index,structure)
        conf_df=compute_distance(pdbfilename,index,conf_df,structure)
        conf_df=compute_dihedrals(pdbfilename,index,conf_df,structure)
        conf_df=chelix_conformation(index,conf_df)
        conf_df=spatial_label(index,conf_df)
        conf_df=dihedral_label(index,conf_df,0.4)
        conf_df=chelix_conformation(index,conf_df)
        print('Input'.rjust(13)+'Model'.rjust(6)+'Chain'.rjust(6)+'B3-Lys'.rjust(7)+'C-helix-Glu'.rjust(12)+'DFG-Phe'.rjust(8)+'Spatial_label'.rjust(14)+\
                                  'Dihedral_label'.rjust(15)+'C-helix_label'.rjust(14))
        print(pdbfilename.rjust(13)+conf_df.at[index,'Model_id'].rjust(6)+conf_df.at[index,'Chain_id'].rjust(6)+\
                              str(int(conf_df.at[index,'Lys_num'])).rjust(7)+str(int(conf_df.at[index,'Glu_num'])).rjust(12)+str(int(conf_df.at[index,'Phe_num'])).rjust(8)+conf_df.at[index,'Spatial_label'].rjust(14)+\
                              conf_df.at[index,'Dihedral_label'].rjust(15)+conf_df.at[index,'Chelix'].rjust(14))
    
    elif align.upper()=='TRUE':                
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
                    print(f'Model {model.id}, Chain {chain.id} is probably not a protein kinase.\n')
                
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
                        
                        conf_df=compute_distance(pdbfilename,index,conf_df,structure)
                        conf_df=compute_dihedrals(pdbfilename,index,conf_df,structure)
                        conf_df=chelix_conformation(index,conf_df)
                        conf_df=spatial_label(index,conf_df)
                        conf_df=dihedral_label(index,conf_df,0.4)
                        conf_df=chelix_conformation(index,conf_df)
                        
                        if len(chain_list)==1:
                            print('Input'.rjust(13)+'Group'.rjust(6)+'Model'.rjust(6)+'Chain'.rjust(6)+'B3-Lys'.rjust(7)+'C-helix-Glu'.rjust(12)+'DFG-Phe'.rjust(8)+'Spatial_label'.rjust(14)+\
                                  'Dihedral_label'.rjust(15)+'C-helix_label'.rjust(14))
                        print(pdbfilename.rjust(13)+conf_df.at[index,'Group'].rjust(6)+conf_df.at[index,'Model_id'].rjust(6)+conf_df.at[index,'Chain_id'].rjust(6)+\
                              str(int(conf_df.at[index,'Lys_num'])).rjust(7)+str(int(conf_df.at[index,'Glu_num'])).rjust(12)+str(int(conf_df.at[index,'Phe_num'])).rjust(8)+conf_df.at[index,'Spatial_label'].rjust(14)+\
                              conf_df.at[index,'Dihedral_label'].rjust(15)+conf_df.at[index,'Chelix'].rjust(14))
                        delete_files(pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'])
                        
    

if __name__ == '__main__':
    hmm_loc=os.path.dirname(os.path.realpath(__file__))+'/HMMs'    #gets the original location of the file
    parser=argparse.ArgumentParser()
    parser.add_argument('PDB',help='PDB file in .cif or .pdb format. You can also provide compressed file as file.cif.gz or file.pdb.gz')
    parser.add_argument('Align',help='True/False; If "True" then the program will align the input protein sequence to pre-computed HMMs and identify conserved \
                        residues (Example $python kinase_state.py 1GAG.pdb True).\
                        If "False" then the please provide the optional arguments (Example $python kinase_state.py --Model 0 --Chain A --Lys 1030 --Glu 1047 --Phe 1151 1GAG.pdb False)')
    parser.add_argument('--Model',help='Model id of the protein kinase structure in the input file, if the file does not contain any Model id then specify 0')
    parser.add_argument('--Chain',help='Chain id of the protein kinase chain in the input file')
    parser.add_argument('--Lys',help='Residue number of the conserved B3-Lys in the input PDB file')
    parser.add_argument('--Glu',help='Residue number of the conserved C-helix-Glu in the input PDB file')
    parser.add_argument('--Phe',help='Residue number of the conserved DFG-Phe in the input PDB file')
    args=parser.parse_args()
   
    pdbfilename=args.PDB
    if args.Align=='False':
        user_model=str(args.Model);user_chain=str(args.Chain);user_lys=int(args.Lys);user_glu=int(args.Glu);user_phe=int(args.Phe)
    if args.Align=='True':    
        user_model=user_chain=user_lys=user_glu=user_phe=align=''
   
    identify_state(hmm_loc,pdbfilename,args.Align,user_model,user_chain,user_lys,user_glu,user_phe)
    