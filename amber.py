import templates
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
 

leap = ''' set default pbradii mbondi3
source leaprc.protein.ff14SBonlysc
aa = sequence {{ {0} }}
saveamberparm aa {1}.top {1}.rst7
quit
'''

minimize_in = '''Minimize a protein using GB
&cntrl
 imin=1, maxcyc=1000, ncyc=500,
 cut=999., rgbmax=999.,igb=8, ntb=0,
 ntpr=100,
/
'''


minimize_gpu = '''source ~/.load_Amber
pmemd.cuda -O \
    -i min.in \
    -o min_{0}.out  \
    -c {0}.rst7  \
    -ref {0}.rst7  \
    -r min_{0}.rst  \
    -p {0}.top

#Now create a pdb from the restart
cpptraj {0}.top <<EOF
trajin min_{0}.rst 
trajout min_{0}.pdb pdb
go
    '''



def sequence2amber(sequence,sequence_name):
    '''Gets a text string with 3 letter code for each amino acid, 
    prepares it for leap input and runs tleap'''

    sequence = "\n".join(textwrap.wrap(sequence,80))
    with open('tleap.in','w') as fo:
        fo.write(leap.format(sequence,sequence_name))
    os.system('tleap -f tleap.in')

def minimizeGPU(sequence_name):
    '''Creates input minimize script, input in and submits job'''
    txt = templates.single_GPU_head.format(sequence_name)
    txt += minimize_gpu.format(sequence_name)
    with open('min.in','w') as fo:
        fo.write(minimize_in)
    with open('{}.sh'.format(sequence_name),'w') as fo:
        fo.write(txt)
    print('submit job\n sbatch {}.sh'.format(sequence_name))





'''Examples of how to use this module

sequence2amber('MET GLU PHE THR VAL SER THR THR GLU ASP LEU GLN ARG TYR ARG THR GLU CYS VAL SER SER LEU ASN ILE PRO ALA ASP TYR VAL GLU LYS PHE LYS LYS TRP GLU PHE PRO GLU ASP ASP THR THR MET CYS TYR ILE LYS CYS VAL PHE ASN LYS MET GLN LEU PHE ASP ASP THR GLU GLY PRO LEU VAL ASP ASN LEU VAL HIS GLN LEU ALA HIS GLY ARG ASP ALA GLU GLU VAL ARG THR GLU VAL LEU LYS CYS VAL ASP LYS ASN THR ASP ASN ASN ALA CYS HIS TRP ALA PHE ARG GLY PHE LYS CYS PHE GLN LYS ASN ASN LEU SER LEU ILE LYS ALA SER ILE LYS LYS ASP','sequence_A')
minimizeGPU('sequence_A')

or

NEF = pickle.load(open('/ufrc/alberto.perezant/alberto.perezant/NEF/Forked_NEF/NEF/data_1_1/PDBStat_developers/Perez/trial.nef.pkl','rb'))
for seq in NEF.sequence_names:
    sequence2amber(NEF.sequence[seq],seq)
    minimizeGPU(seq)
'''


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

    for seq in NEF.sequence_names:
        sequence2amber(NEF.sequence[seq],seq)
        minimizeGPU(seq)
    


if __name__ == '__main__': #Python way to execute main()
    main()


