import hipergator
import textwrap
import os
import pickle
from MELD_NEF_Classes import *
import argparse

def parse_args():                              #in line argument parser with help 
    '''
    Parse arguments for the program from command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-nef', type=str, help='NEF input file or object')
    parser.add_argument('-directory', type=str, help='Where to place the files',default='.')
    return(parser.parse_args())
 

map_to_heavy = {
        ('ARG','HG%'):(['CG'],1),
        ('ARG','HGx'):(['CG'],1),
        ('ARG','HGy'):(['CG'],1),
        ('ARG','HD%'):(['CD'],1),
        ('ARG','HDx'):(['CD'],1),
        ('ARG','HDy'):(['CD'],1),
        ('ARG','HH1%'):(['NH1'],1),
        ('ARG','HH2%'):(['NH2'],1),
        ('ARG','HH%'):(['NH1','NH2'],1),
        ('ARG','NZ%'):(['NH1','NH2'],0),
        ('ASN','HD2%'):(['ND2'],1),
        ('ASN','HD2x'):(['ND2'],1),
        ('ASN','HD2y'):(['ND2'],1),
        ('ASP','OD%'):(['OD1','OD2'],0),
        ('GLN','HG%'):(['CG'],1),
        ('GLN','HE2%'):(['OE1'],1),
        ('GLN','HE2x'):(['OE1'],1),
        ('GLN','HE2y'):(['OE1'],1),
        ('GLN','HGx'):(['CG'],1),
        ('GLN','HGy'):(['CG'],1,),
        ('GLU','HG%'):(['CG'],1),
        ('GLU','HGx'):(['CG'],1),
        ('GLU','HGy'):(['CG'],1),
        ('GLU','OE%'):(['OE1','OE2'],0),
        ('GLY','HA%'):(['CA'],1),
        ('GLY','HAx'):(['CA'],1),
        ('GLY','HAy'):(['CA'],1),
        ('ILE','HG1%'):(['CG1'],1),
        ('ILE','HG2%'):(['CG2'],1),
        ('ILE','HD1%'):(['CD1'],1),
        ('ILE','HG1x'):(['CG1'],1),
        ('ILE','HG1y'):(['CG1'],1),
        ('LEU','HD1%'):(['CD1'],1),
        ('LEU','HD2%'):(['CD2'],1),
        ('LEU','HD%'):(['CD1','CD2'],1),
        ('LYS','HG%'):(['CG'],1),
        ('LYS','HD%'):(['CD'],1),
        ('LYS','HE%'):(['CE'],1),
        ('LYS','HZ%'):(['NZ'],1),
        ('LYS','HDx'):(['CD'],1),
        ('LYS','HDy'):(['CD'],1),
        ('LYS','HEx'):(['CE'],1),
        ('LYS','HEy'):(['CE'],1),
        ('LYS','HGx'):(['CG'],1),
        ('LYS','HGy'):(['CG'],1),
        ('MET','HG%'):(['CG'],1),
        ('MET','HE%'):(['CE'],1),
        ('MET','HGx'):(['CG'],1),
        ('MET','HGy'):(['CG'],1),
        ('PHE','HD%'):(['CD1','CD2'],1),
        ('PHE','HE%'):(['CE1','CE2'],1),
        ('PHE','CD%'):(['CD1','CD2'],0),
        ('PHE','CE%'):(['CE1','CE2'],0),
        ('PRO','HG%'):(['CG'],1),
        ('PRO','HD%'):(['CD'],1),
        ('PRO','HDx'):(['CD'],1),
        ('PRO','HDy'):(['CD'],1),
        ('PRO','HGx'):(['CG'],1),
        ('PRO','HGy'):(['CG'],1),
        ('THR','HG2%'):(['CG2'],1),
        ('TYR','HD%'):(['CD'],1),
        ('TYR','HE%'):(['CE'],1),
        ('TYR','CD%'):(['CD1','CD2'],0),
        ('TYR','CE%'):(['CE1','CE2'],0),
        ('VAL','HG1%'):(['CG1'],1),
        ('VAL','HG2%'):(['CG2'],1),
        ('VAL','HG%'):(['CG1','CG2'],1),
        }

aminoacids = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]

for aa in aminoacids:
    map_to_heavy[(aa,'HB%')] = (['CB'],1)
    map_to_heavy[(aa,'HBx')] = (['CB'],1)
    map_to_heavy[(aa,'HBy')] = (['CB'],1)



def write_peaks(peaks,min_CO=4):
    '''We get peaks list and write it out. If we did process_peaks before then the peaks
    will be mapped to heavy atoms. if we used process_sequence, atom ordering will be amber-like

    We cannot separate in chains, some restraints could be ambiguos and satisfied in either monomer 
    or dimer

    Peaks that are trivial to satisfy (e.g. same residue) should remove the same ambiguous peak that defines it.
    If any possible contact in a peak is trivial, remove the whole peak.
    '''
    output_noe = ''
    for one_peak in peaks.restraint_id.unique():
        selection =  peaks.loc[ peaks['restraint_id'] == one_peak, ['sequence_code_1','atom_name_1','sequence_code_2','atom_name_2','upper_limit'] ]
        trivial = selection.loc[ abs(selection['sequence_code_1'] - selection['sequence_code_2']) < min_CO]
        print(len(trivial))
        if len(trivial < 1):
            print('again',len(trivial))
            print(selection.to_string(index=False,header=False))
            output_noe += selection.to_string(index=False,header=False) 
            output_noe += "\n\n"

    print(output_noe)


def process_sequence(NEF,peaks):
    '''USe mapping of residues from molecular system block to pass into amber-like topology numbering.
    First residue starts at 1, new chains keep numbering intead of starting over
    Here we run into an issue: pandas has dataframe.index as their internal numbering...starts at 0
    NEF sets the name of the column index which would start at 1 as we would like
    '''
    data = NEF.block_content['molecular_system'].loop_type_data['_nef_sequence']
    numbering = {}
    print(data)
    for seq,chain,index in zip(data.sequence_code,data.chain_code,data.index):
        numbering[(seq,chain)] = index + 1

    for (seq,chain) in numbering.keys():
        peaks.loc[ (peaks['chain_code_1'] == chain) & (peaks['sequence_code_1'] == seq),['sequence_code_1'] ] = numbering[(seq,chain)]
        peaks.loc[ (peaks['chain_code_2'] == chain) & (peaks['sequence_code_2'] == seq),['sequence_code_2'] ] = numbering[(seq,chain)]

    return(peaks)



def process_peaks(peaks):
    ''' Translate the peaks into MELD restraints. Put ambiguous H distances on the heavy atom
    add one to the upper distance treshold
    We will use a dictionary with all possible cases and replace by the heavy atom counterparts
    in cases whereabmiguity leads to more than one heavy atom, we will duplicate the row in the 
    pandas dataframe'''

    #Correct first atom
    for (amino, atom) in map_to_heavy.keys():
        (heavy,distance_correction) = map_to_heavy[(amino,atom)]
        #If have more than one heavy atom need to duplicate (or more) data
        aa = peaks.loc[(peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom)]
        if len(aa) < 1:
            continue
        if len(heavy)== 1:
            peaks.loc[ (peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom),['upper_limit'] ] = peaks.upper_limit + distance_correction
            peaks.loc[ (peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom),['atom_name_1'] ] = heavy[0]
        else:
            '''First pass changes the first ambiguos atom to heavy. Then adds the original ones to teh data frame and substitutes by the second atom'''
            peaks.loc[ (peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom),['upper_limit'] ] = peaks.upper_limit + distance_correction
            peaks.loc[ (peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom),['atom_name_1'] ] = heavy[0]
            peaks.append(aa)            
            peaks.loc[ (peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom),['upper_limit'] ] = peaks.upper_limit + distance_correction
            peaks.loc[ (peaks['residue_name_1'] == amino) & (peaks['atom_name_1'] == atom),['atom_name_1'] ] = heavy[1]
    #Correct second atom
    #TO DO: This is duplicated code, should create routine to handle
    for (amino, atom) in map_to_heavy.keys():
        (heavy,distance_correction) = map_to_heavy[(amino,atom)]
        #If have more than one heavy atom need to duplicate (or more) data
        aa = peaks.loc[(peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom)]
        if len(aa) < 1:
            continue
        if len(heavy)== 1:
            peaks.loc[ (peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom),['upper_limit'] ] = peaks.upper_limit + distance_correction
            peaks.loc[ (peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom),['atom_name_2'] ] = heavy[0]
        else:
            '''First pass changes the first ambiguos atom to heavy. Then adds the original ones to teh data frame and substitutes by the second atom'''
            peaks.loc[ (peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom),['upper_limit'] ] = peaks.upper_limit + distance_correction
            peaks.loc[ (peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom),['atom_name_2'] ] = heavy[0]
            peaks.append(aa)            
            peaks.loc[ (peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom),['upper_limit'] ] = peaks.upper_limit + distance_correction
            peaks.loc[ (peaks['residue_name_2'] == amino) & (peaks['atom_name_2'] == atom),['atom_name_2'] ] = heavy[1]
    return(peaks)






'''Examples of how to use this module
'''

NEF = pickle.load(open('/ufrc/alberto.perezant/alberto.perezant/NEF/Forked_NEF/NEF/data_1_1/PDBStat_developers/Perez/trial.nef.pkl','rb'))
NEF.peaks = {}
for i in NEF.chains:
    for j in NEF.chains:
        NEF.peaks[(i,j)] = []

for NOE in NEF.block_types['distance_restraint_list']:
    print(dir(NOE))
    distances = NOE.loop_type_data['_nef_distance_restraint']
    #Now we want to go through chains I-I and I-J
    print(NEF.sequence_names)
    distances = process_peaks(distances)
    distances = process_sequence(NEF,distances)
    write_peaks(distances)
    die
    NOE.loop_type_data['_nef_distance_restraint'] = process_peaks(distances)
    distances = NOE.loop_type_data['_nef_distance_restraint']
    write_peaks(distances,NEF.chains)


def main():
    args = parse_args()
    #Work in a temporary directory
    if args.directory == '.':
        args.directory = os.getcwd()
    try:
        NEF = pickle.load(open(args.nef,'rb'))
    except:
        NEF = NEF_system(args.nef,args.directory)
        NEF.sequence()
    
    print(NEF.block_types.keys())
    #we are interested in 'distance_restraint_list', 'dihedral_restraint_list'
    for NOE in NEF.block_types['distance_restraint_list']:
        pass
    #distance2MELD(NEF


#if __name__ == '__main__': #Python way to execute main()
    #main()


