#!/usr/bin/env python3
from Bio.PDB.DSSP import *
import sys
import os
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import BiopythonWarning
from collections import Counter
warnings.simplefilter('ignore', BiopythonWarning)

'''
This script is used to get the RASA values for the amino acids and 
Output a line graph with the RASA values for each protein
calculate fraction of residues below a certain threshold (0.15 - 0.25)
Find the residues that have less than 0.18 threshold and calculate the occurrences
'''

fileExe = r"./dssp-2.0.4-win32.exe"
#folderRead = r"./barebones_pdb/RIFT"
folderRead = r"./pdb_by_fold/Zn"

asaName = sys.argv[1]
resList = []

def get_DSSPList(file):
	p = PDBParser() 
	'''
	# parses the pdb file
	'''
	s = p.get_structure('X', file) 
	'''
	# getting the structure 
	'''
	model = s[0]
	d = DSSP(model, file,dssp=fileExe,acc_array=asaName) 
	'''
	# DSSP executable
	'''
	dssp_dict,dssp_keys = dssp_dict_from_pdb_file(file,DSSP=fileExe) 
	'''
	# Create a dssp dictionary from a PDB file
	'''
	#print (file)
	
	dssp_list = []
	a_keys =list(d.keys())
	#print(a_keys)
	for v in a_keys:
		rasa_values = (d[v]) 
		'''
		#values of the dictionary that gives the RASA values
		'''
		acc_values = dssp_dict[v] 
		
		'''# from the dictionary that provides the residue number, aa, and the acc value'''
		
		x = v[1][1],acc_values[0],acc_values[2],rasa_values[3] 
		
		'''# residue number, amino acid, acc value and RASA value'''
		
		dssp_list.append(x) 
		
		'''# creating a list of these values'''
		
	return dssp_list

def do_Graph2(dssp_list):
	
	''' creating two lists x and y to input the residue number in x and the RASA value in y'''
	
	x_axis = [] 
	y_axis = []
	residue = []
	
	for z,a,b,y in dssp_list:
		
		
		x_axis.append(z)
		y_axis.append(y)
		res = z, y, a
		residue.append(res)
	

	''' Determining the threshold RASA value '''
	q = len(y_axis) ''' Total number of amino acids '''
	# alist = [0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25]
	# for e in alist:
	p = 0
	aaList = []
	
	for spot, rasa, aa in residue:
		#print(rasa, file)
		if rasa != 'NA':
			#if rasa < e: ''' Using e as the RASA threshold '''
			if rasa < 0.18:  '''Using 0.18 as the threshold '''
				p+=1 ''' Number of amino acids that have an RASA value below 0.18 '''
				temp = spot, aa
				aaList.append (temp)
			fraction_below = p/q  '''fraction of amino acids that have RASA value below the 0.18 threshold'''
		else:
			fraction_below = 'NA'
	#print (aaList)
	return (aaList)

	

''' summing up all residues in a particular fold (from multiple files)'''
def sumResidues(resList):
    
    all_resList = []		
    for sublist in resList:
            for item in sublist:
                    all_resList.append(item)
    #print (all_resList)
    c = Counter(all_resList) ''' Gives the count of each residue '''
    for x,y in sorted(c.items()):
            print (x, y)
			
			
def residueProtein(dssp_list):
	res_aa=[]
	res_Dict = {}
	
	for z,a,b,y in dssp_list:
		
		res_aa.append(a)
	c = Counter(res_aa)
	#print (c)
	#print(file)
	for k, v in c.items():
	#print(k, v)
		res_Dict[k] = v/(len(res_aa))
	for k in res_Dict:
		
		print (k, res_Dict[k])
	print("\n")
			

def main():
    for root, subdirs, files in os.walk(folderRead):
            for file in files:
                    # #print (file)
            
                    if file.endswith('.pdb'):
                            output = file+".out"
                            outgraph = file+".png"
                            rasagraph = os.path.join(folderRead, outgraph)
                            outfile = os.path.join(folderRead, output)
                            #rasagraph = os.path.join(file, outgraph) '''for reiteration in folders '''
                            #outfile = os.path.join(root, output) ''' reiteration '''
                            dssp_file = open(outfile, "w")
                            #realpath=os.path.join(root, file)  ''' reiteration '''
                            realpath=os.path.join(folderRead, file)
                            dssp_file.write(realpath)
                            dssp_file.write("\n")
                            dssp_file.write("\n")
                            
                            dsspList = get_DSSPList(realpath)
                            dssplist=get_DSSPList(realpath)
                            sumList = do_Graph(dsspList)
                            resList.append(sumList)
                            #create_Graph(dssp_file)
                            dssp_file.close()
    ''' Calling Sum Residues '''
    sumResidues(resList)


if __name__ == "__main__":
    main()                           







