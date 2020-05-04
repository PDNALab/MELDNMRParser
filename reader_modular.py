#! /ufrc/alberto.perezant/alberto.perezant/NMR_Format/Xplor-NIH/xplor-nih-2.53/bin/pyXplor
import nefTools as nef
import cif
import protocol 
import numpy


header = '''
data_nef_MELD input

save_nef_nmr_meta_data
    _nef_nmr_meta_data.sf_category     nef_nmr_meta_data
    _nef_nmr_meta_data.sf_framecode    nef_nmr_meta_data
    _nef_nmr_meta_data.format_name     nmr_exchange_format
    _nef_nmr_meta_data.format_version  1.1
    _nef_nmr_meta_data.program_name    MELD / Xplor-NIH
    _nef_nmr_meta_data.program_version 1.0 / 2.53
    _nef_nmr_meta_data.creation_date   2020-03-31_20:56:19
    _nef_nmr_meta_data.uuid            Xplor-NIH_2020-03-31_20:56:19_EBBC74D8
save_
'''



#pdb = protocol.loadPDB('1nk2.pdb')

data = nef.readNEF('CCPN_Commented_Example.nef')

#if we wanted the least amount of processing.This just reads the nef file. and writes back the text
nef_string =  data.asString()
open('new.nef','w').write(nef_string)

#Types of data
dtypes = data.keys()
#['nef_my_nmr_project_1', 'nef_nmr_meta_data', 'nef_molecular_system', 'nef_chemical_shift_list_1', 'nef_chemical_shift_list_2', 'nef_distance_restraint_list_L1', 'nef_distance_restraint_list_hbond1', 'nef_dihedral_restraint_list_L2', 'nef_rdc_restraint_list_3', 'nef_nmr_spectrum_cnoesy1', 'nef_nmr_spectrum_dummy15d', 'nef_peak_restraint_links']
'''
Here each key tells us about a record:
nef_my_nmr_project_1 --> this is the project name. First and last line of the file
nef_nmr_meta_data --> where the data comes from
nef_molecular_system --> our sequence. Xplor uses PSF kind of building, so it tells us how the Aa are patched from respect reference
then we can have severa; of catPrefixes = {'dihedral': 'nef_dihedral_restraint_list', 'distance': 'nef_distance_restraint_list', 'shifts': 'nef_chemical_shift_list', 'spectrum': 'nef_nmr_spectrum'}

For example, in this case we have nef_chemical_shift_list_1 and nef_chemical_shift_list_2 as examples of shift.
If there is only one, its name is by default. Otherwise, have to get it from the string name (1 or 2 in this case)
now we could read that block as bb = nef.getBlock(data, blockType='shifts',nefRestraintName='1)

Same goes for distance --> we will take a look at how many strings with "nef_distance_restraint_list" in this case 
L1, hbond1, L2 

Next up is dihedral. We only have one, so no need to specify a name  (it would be L2).

Finally, we have the spectrum --> cnoesy1 dummy15d

Notice that we have some data that would not have an entry right now: nef_rdc_restraint_list_3 nef_peak_restraint_links

catPrefixes = {'dihedral': 'neftring +=  data.asString()dihedral_restraint_list', 'distance': 'nef_distance_restraint_list', 'shifts': 'nef_chemical_shift_list', 'spectrum': 'nef_nmr_spectrum'}
dataTypes = ['dihedral','distance','shifts','spectrum']


distance -> nef_chemical_distance_list 
    ['datablock', 'nef_distance_restraint_list', 'nef_distance_restraint']
shifts -> nef_chemical_shift_list 
dihedral -> nef_dihedral_restraint
             aa['nef_dihedral_restraint_list']['sf_category'] sf_framecode
             datablock -> name

spectrum --> ['datablock', 'nef_nmr_spectrum', 'nef_spectrum_dimension', 'nef_spectrum_dimension_transfer', 'nef_peak']
bb['nef_nmr_spectrum'].keys()
['sf_category', 'sf_framecode', 'num_dimensions', 'chemical_shift_list', 'experiment_classification', 'experiment_type']

bb['nef_spectrum_dimension_transfer'].keys()
['dimension_1', 'dimension_2', 'transfer_type', 'is_indirect']

bb['nef_peak'].keys()
['index', 'peak_id', 'volume', 'volume_uncertainty', 'height', 'height_uncertainty', 'position_1', 'position_uncertainty_1', 'position_2', 'position_uncertainty_2', 'position_3', 'position_uncertainty_3', 'chain_code_1', 'sequence_code_1', 'residue_name_1', 'atom_name_1', 'chain_code_2', 'sequence_code_2', 'residue_name_2', 'atom_name_2', 'chain_code_3', 'sequence_code_3', 'residue_name_3', 'atom_name_3']


'''


blockType = {'shifts':'nef_chemical_shift_list',
             'distance':'nef_distance_restraint_list',
             'dihedral':'nef_dihedral_restraint_list',
             'spectrum':'nef_nmr_spectrum'}

properties = {'shifts':['nef_chemical_shift'],
             'distance':['nef_distance_restraint'],
             'dihedral':['nef_dihedral_restraint'],
             'spectrum':['nef_spectrum_dimension','nef_spectrum_dimension_transfer','nef_peak']}


def process_blocks(data,block=None,my_key=None,properties=None,fhand=None):
    dtypes = data.keys()
    my_str = ''
    for i in dtypes:
        if my_key in i:
            print(i,my_key,properties)
            name =  i.split('{}_'.format(my_key))[1]
            temp = nef.getBlock(data, blockType=block,nefRestraintName=name)
            # ['chain_code', 'sequence_code', 'residue_name', 'atom_name', 'value', 'value_uncertainty', 'element', 'isotope_number']
            #Here in are all the arrays we might need if we were to process shift data
            #you have one such array for each property
            #We could add, remove, ... That would allow us to write our own data.
            #First we specify group to work with
            fo.write('save_{}\n'.format(i))
            #we go through all the categories in the group
            for category in temp[my_key].keys():
                print(my_key,category,temp[my_key][category][0])
                print('_{}.{}     {}\n'.format(my_key,category,temp[my_key][category][0]))
                aa = '_{}.{}     {}\n'.format(my_key,category,temp[my_key][category][0])
                print(aa)
                fo.write('_{}.{}     {}\n'.format(my_key,category,temp[my_key][category][0]))
                my_str += aa
            #now we loope through all the properties it has with its data
            tmp_prop = []
            for p in properties:
                if p in temp.keys():
                    tmp_prop.append(p)
            for my_property in tmp_prop:
                fo.write('loop_\n')
                my_data = []
                for label in temp[my_property].keys():
                    fo.write('_{}.{}\n'.format(my_property,label))
                    my_data.append( temp[my_property][label] )
                my_data = numpy.array(my_data)
                print(my_data)
                for j in range(len(my_data[0,:])):
                    fo.write("{}\n".format(" ".join(my_data[:,j]) ) )
                fo.write("stop_\n")
            fo.write("save_\n")

nef_string = nef.genHeader(datablockName='MELD input')
with open('out.cif','w') as fo:
    #fo.write(header)
    fo.write(nef_string)
    for block in blockType.keys():
        process_blocks(data,block=block,my_key=blockType[block],properties=properties[block],fhand=fo)

#block = 'dihedral'
#print(process_blocks(data,block=block,my_key=blockType[block],properties=properties[block]))

#with open('out.cif','w') as fo:
#    fo.write(nef_string)


