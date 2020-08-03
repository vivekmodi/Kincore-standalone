#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 09:58:48 2020

@author: vivekmodi
"""

import numpy as np
from Bio import PDB

def compute_distance(pwd,pdbfilename,index,conf_df):
    restype_atom_dict={'F':'CZ','L':'CG','P':'CG','M':'CE','S':'OG','V':'CB','W':'CZ3','Y':'OH','A':'CB'}   #Add all the twenty residues here
    phe_atom_type=restype_atom_dict[conf_df.at[index,'Phe_restype']]
    
    #Distance Glu4-Phe
    conf_df.at[index,'Glu4-Phe-dis']=distance_atoms(pwd,pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                                     conf_df.at[index,'Glu4_num'],conf_df.at[index,'Phe_num'],'CA',phe_atom_type)      #Change atom names for other residue types
        
    #Distance Lys-Phe
    conf_df.at[index,'Lys-Phe-dis']=distance_atoms(pwd,pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                                    conf_df.at[index,'Lys_num'],conf_df.at[index,'Phe_num'],'CA',phe_atom_type)      #Change atom names for other residue types
    
    #Distance Lys-Glu
    conf_df.at[index,'Lys-Glu-dis']=distance_atoms(pwd,pdbfilename,conf_df.at[index,'Model_id'],conf_df.at[index,'Chain_id'],\
                                    conf_df.at[index,'Lys_num'],conf_df.at[index,'Glu_num'],'CB','CB')      #Change atom names for other residue types
    
    return conf_df


def distance_atoms(pwd,pdbfilename,model_id,chain_id,res1,res2,atm1,atm2):
    
    if '.pdb' in pdbfilename.lower():
        parser=PDB.PDBParser(QUIET=True)
    if '.cif' in pdbfilename.lower():
        parser=PDB.MMCIFParser(QUIET=True)
    structure=parser.get_structure('PDB',f'{pwd}/{pdbfilename}')
    atom_present=0; res1=int(res1); res2=int(res2)
   
    for model in structure:
        for chain in model:
            insertion_num=0    #Count residues with insertion codes and skip them
            if str(model.id)==model_id and chain.id==chain_id:
                
                for residue in chain:
                    if residue.get_id()[0]==' ' and residue.get_id()[2]!=' ':      #Insertion code present
                        insertion_num+=1
                    
                    if int(residue.id[1])==(int(res1)-insertion_num) and residue.get_id()[0]==' ':
                        if residue.has_id(atm1):
                            residue1=chain[res1-insertion_num]
                            atom_present=atom_present+1
                            
                    if int(residue.id[1])==(int(res2)-insertion_num) and residue.get_id()[0]==' ':
                        if residue.has_id(atm2):
                            residue2=chain[res2-insertion_num]
                            atom_present=atom_present+1
                            

    if atom_present==2:
        distance=np.round((residue1[atm1]-residue2[atm2]),2)
        return distance
    else:
        return 99999
