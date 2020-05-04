#! /ufrc/alberto.perezant/alberto.perezant/NMR_Format/Xplor-NIH/xplor-nih-2.53/bin/pyXplor
import nefTools as nef
import cif
import protocol 
import numpy

pdb = protocol.loadPDB('1nk2.pdb')

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
'''


nef_string = nef.genHeader(datablockName='MELD input')
print(nef_string)

#Let's process shifts:
my_key = 'nef_chemical_shift_list'
for i in dtypes:
    if my_key in i:
        name =  i.split('{}_'.format(my_key))[1]
        temp = nef.getBlock(data, blockType='shifts',nefRestraintName=name)
        print(temp['nef_chemical_shift'].keys())
        # ['chain_code', 'sequence_code', 'residue_name', 'atom_name', 'value', 'value_uncertainty', 'element', 'isotope_number']
        #Here in are all the arrays we might need if we were to process shift data
        print(temp['nef_chemical_shift']['atom_name'])
        #you have one such array for each property

        #We could add, remove, ... That would allow us to write our own data.
        #nef_string += nef.makeCif(temp) ### this does not work, just generates the sequence data and overall header!
        print('save_{}'.format(temp['datablock']['name']))
        print('_{}.{}     {}'.format(my_key,'sf_category',temp[my_key]['sf_category']))
        print('_{}.{}     {}'.format(my_key,'sf_framecode',temp[my_key]['sf_framecode']))
        print('loop_')

        my_data = []
        for label in temp['nef_chemical_shift'].keys():
            print('_{}.{}'.format('nef_chemical_shift',label))
            my_data.append( temp['nef_chemical_shift'][label] )

        my_data = numpy.array(my_data)
        #print(my_data)
        #print by rows
        for j in range(len(my_data[0,:])):
            print(" ".join(my_data[:,j]))

print("stop_")
print("save_")




#catPrefixes = {'dihedral': 'neftring +=  data.asString()dihedral_restraint_list', 'distance': 'nef_distance_restraint_list', 'shifts': 'nef_chemical_shift_list', 'spectrum': 'nef_nmr_spectrum'}
#dataTypes = ['dihedral','distance','shifts','spectrum']

#for d in dataTypes:
#    try:
#        print( nef.getBlock(data, blockType=d) )
#    except:
#        print( '{} type needs extra argument'.format(d) ) 



#nef_string += nef.getBlock(data, blockType='dihedral')

#open('outFile.nef','w').write(nef_string)

#nef_string += nef.makeCif(datablockName='dihedral').asString()

#open('outFile.nef','w').write(nef_string)

#output = cif.Cif()
#output.addDatablock(nef.genHeader(datablockName='MELD input') ) 
#output.addDatablock('dihedral',  nef.getBlock(data, blockType='dihedral') )
#with open('out.cif','w') as fo:
#    output.write(fo)
