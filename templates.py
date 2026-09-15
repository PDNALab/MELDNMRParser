
single_GPU_head = '''#!/bin/bash
#SBATCH --job-name={0}      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1      
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1           
#SBATCH --partition=gpu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --reservation=perez
#SBATCH --mem-per-cpu=2000mb          # Memory per processor
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output={0}_%j.log     # Standard output and error log
'''


meld_NMR_script = '''#!/usr/bin/env python
# encoding: utf-8

#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import meld
from meld.remd import ladder, adaptor, leader
import meld.system.montecarlo as mc
from meld import system
from meld.system import patchers
from meld import comm, vault
from meld import parse
from meld import remd
from meld.system import param_sampling
from openmm import unit as u

from meld.system.scalers import LinearRamp,ConstantRamp
from collections import namedtuple
import glob as glob

N_REPLICAS = 30
N_STEPS = 1000     # 50 ns (1000 x 3.5 X 14286)
BLOCK_SIZE = 100


hydrophobes = 'AILMFPWV'
hydrophobes_res = ['ALA','ILE','LEU','MET','PHE','PRO','TRP','VAL']


def gen_state(s, index):
    state = s.get_state_template()
    state.alpha = index / (N_REPLICAS - 1.0)
    return state

def make_ss_groups(subset=None):
           active = 0
           extended = 0
           sse = []
           ss = open('ss.dat','r').readlines()[0]
           for i,l in enumerate(ss.rstrip()):
               #print i,l
               if l not in "HE.":
                   continue
               if l not in 'E' and extended:
                   end = i
                   sse.append((start+1,end))
                   extended = 0
               if l in 'E':
                   if i+1 in subset:
                       active = active + 1
                   if extended:
                       continue
                   else:
                       start = i
                       extended = 1
           print(active,':number of E residues')
           print(sse,':E residue ranges')
           return sse,active

def create_hydrophobes(s,group_1=np.array([]),group_2=np.array([]),CO=True):
    hy_rest=open('hydrophobe.dat','w')
    atoms = {"ALA":['CA','CB'],
             "VAL":['CA','CB','CG1','CG2'],
             "LEU":['CA','CB','CG','CD1','CD2'],
             "ILE":['CA','CB','CG1','CG2','CD1'],
             "PHE":['CA','CB','CG','CD1','CE1','CZ','CE2','CD2'],
             "TRP":['CA','CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
             "MET":['CA','CB','CG','SD','CE'],
             "PRO":['CD','CG','CB','CA']}
    #Groups should be 1 centered
    n_res = s.residue_numbers[-1]
    print(n_res)
    group_1 = group_1 if group_1.size else np.array(list(range(n_res)))+1
    group_2 = group_2 if group_2.size else np.array(list(range(n_res)))+1

    #Get a list of names and residue numbers, if just use names might skip some residues that are two
 #times in a row
    #make list 1 centered
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    print(sequence)
    sequence = dict(sequence)

    print(group_1)
    print(group_2)
    group_1 = [ res for res in group_1 if (sequence[res-1] in hydrophobes_res) ]
    group_2 = [ res for res in group_2 if (sequence[res-1] in hydrophobes_res) ]

    print(group_1)
    print(group_2)
    pairs = []
    hydroph_restraints = []
    for i in group_1:
        for j in group_2:

            # don't put the same pair in more than once
            if ( (i,j) in pairs ) or ( (j,i) in pairs ):
                continue

            if ( i ==j ):
                continue

            if (abs(i-j)< 7):
                continue
            pairs.append( (i,j) )

            atoms_i = atoms[sequence[i-1]]  #atoms_i = atoms[sequence[i]]
            atoms_j = atoms[sequence[j-1]]  #atoms_j = atoms[sequence[j]]

            local_contact = []
            for a_i in atoms_i:
                for a_j in atoms_j:
                    hy_rest.write('{} {} {} {}\n'.format(i,a_i, j, a_j))
            hy_rest.write('\n')

def generate_strand_pairs(s,sse,subset=np.array([]),CO=True):
    f=open('strand_pair.dat','w')
    n_res = s.residue_numbers[-1]
    subset = subset if subset.size else np.array(list(range(n_res)))+1
    strand_pair = []
    for i in range(len(sse)):
        start_i,end_i = sse[i]
        for j in range(i+1,len(sse)):
            start_j,end_j = sse[j]

            for res_i in range(start_i,end_i+1):
                for res_j in range(start_j,end_j+1):
                    if res_i in subset or res_j in subset:
                        f.write('{} {} {} {}\n'.format(res_i, 'N', res_j, 'O'))
                        f.write('{} {} {} {}\n'.format(res_i, 'O', res_j, 'N'))
                        f.write('\n')

def get_dist_restraints_hydrophobe(filename, s, scaler, ramp, seq):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=0.5*u.nanometer, r4=0.7*u.nanometer,
                                                 k=250*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists

def get_dist_restraints(filename, s, scaler, linear_ramp, seq):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])-1
            name_i = cols[1]
            j = int(cols[2])-1
            name_j = cols[3]
            dist = float(cols[4])/10. #using variable distance

            rest = s.restraints.create_restraint('distance', scaler, linear_ramp,
                                                 r1=0.0*u.nanometer, r2=0.0*u.nanometer, r3=dist*u.nanometer, r4=(dist+0.2)*u.nanometer,
                                                 k=350*u.kilojoule_per_mole/u.nanometer **2,
                                                 atom1=s.index.atom(i,name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j,name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists

#read 'rotamers.dat' to generate torsion restraints
def get_torsion_restraints(filename, s, constant_scaler,linear_ramp_tor):
    torsion_rests = []
    rotamer_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:        #checks empty line, group ends
            torsion_rests.append(s.restraints.create_restraint_group(rotamer_group, 1))
            rotamer_group = []
        else:               #non-empty line, create torsion restraint
            (res1, at1, res2, at2, res3, at3, res4, at4, rotamer_min, rotamer_max) = line.split()
            rotamer_max = float(rotamer_max)
            rotamer_min = float(rotamer_min)
            rotamer_avg = (rotamer_max+rotamer_min)/2.
            rotamer_sd = abs(rotamer_max - rotamer_min)/2.
            rotamer_rest = s.restraints.create_restraint('torsion', constant_scaler,
                                                 #LinearRamp(0,100,0,1), - ramp start at timestep 0 - end 100, weight 0 to 1
                                                 linear_ramp_tor,
                                                 atom1=s.index.atom(int(res1),at1),
                                                 atom2=s.index.atom(int(res2),at2),
                                                 atom3=s.index.atom(int(res3),at3),
                                                 atom4=s.index.atom(int(res4),at4), #non-linear ramp defined at the bottom,
                                                 phi=rotamer_avg*u.degree, delta_phi=rotamer_sd*u.degree, k=2.5*u.kilojoule_per_mole/u.degree**2)
            rotamer_group.append(rotamer_rest)
    return torsion_rests



#######################################

def setup_system():
    
    # load the sequence (search dynamically in common locations)
    seq_file = glob.glob('*sequence*.dat')
    print('Using sequence file: {}'.format(seq_file))
    sequence = parse.get_sequence_from_AA1(filename=seq_file)
    n_res = len(sequence.split())

    # build the system (search for PDB in common locations)
    pdb_input = glob.glob('*.pdb')
    print('Using pdb file: {}'.format(pdb_input))
    p = meld.AmberSubSystemFromPdbFile(pdb_input)
    build_options = meld.AmberOptions(
      forcefield             ="ff14sbside",
      implicit_solvent_model = 'gbNeck2',
      use_big_timestep       = True,
      cutoff                 = 1.8*u.nanometers,
      remove_com             = False,
      #use_amap = False,
      enable_amap            = False,
      amap_alpha_bias = 1.0,
      amap_beta_bias = 1.0
    )


    builder = meld.AmberSystemBuilder(build_options)
    s = builder.build_system([p]).finalize()
    s.temperature_scaler = system.temperature.GeometricTemperatureScaler(0, 0.4, 300.*u.kelvin, 550.*u.kelvin)

##########################

    ramp = s.restraints.create_scaler('nonlinear_ramp', start_time=1, end_time=200,
                                      start_weight=1e-3, end_weight=1, factor=4.0)
    
    linear_ramp = s.restraints.create_scaler('linear_ramp', start_time=0, end_time=200,
                                      start_weight=0.005, end_weight=1)

    linear_ramp_tor = s.restraints.create_scaler('linear_ramp', start_time=0, end_time=100,
                                      start_weight=0.005, end_weight=0.1)

    seq = sequence.split()
    for i in range(len(seq)):
        if seq[i][-3:] =='HIE': seq[i]='HIS'
    print(seq)
    hydrophobic_res_in_protein=[]
    for i in seq:
        for j in hydrophobes_res:
            if i ==j:
                hydrophobic_res_in_protein.append(i)

    no_hy_res=len(hydrophobic_res_in_protein)
    print(no_hy_res,':number of hydrophobic residue')


    conf_scaler = s.restraints.create_scaler('constant')
    confinement_rests = []
    for index in range(n_res):
        rest = s.restraints.create_restraint('confine', conf_scaler, ramp=linear_ramp, atom_index=s.index.atom(index, 'CA', expected_resname=seq[index][-3:]),
                                             radius=4.5*u.nanometer, force_const=250.0*u.kilojoule_per_mole/u.nanometer **2)
        confinement_rests.append(rest)
    s.restraints.add_as_always_active_list(confinement_rests)

    #
    # Setup Scaler
    #
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    subset1= np.array(list(range(n_res))) + 1

#######------- different scalers ---------
    distance_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=0.8, factor=4.0)
    distance_scaler_short = s.restraints.create_scaler('nonlinear', alpha_min=0.8, alpha_max=1.0, factor=4.0)
    constant_scaler = s.restraints.create_scaler('constant')
##############


    
    ##### ---- NOE --------
    # Load global NOE files dynamically using glob. This will pick up files like
    # `NOE_*.dat` in the working directory or fallback to the example data folder.
    NOESY_files = []
    NOESY_files.extend(glob.glob('NOE_*.dat'))
    NOESY_files = sorted(list(dict.fromkeys(NOESY_files)))
    for noe_file in NOESY_files:
        print('loading {} ...'.format(noe_file))
        NOESY = get_dist_restraints(noe_file, s, distance_scaler, linear_ramp, seq)
        s.restraints.add_selectively_active_collection(NOESY, int(len(NOESY) * 1.00))

    # Load local NOE files (if any) using glob as well.
    NOESY_local_files = []
    NOESY_local_files.extend(glob.glob('local_NOE*.dat'))
    NOESY_local_files = sorted(list(dict.fromkeys(NOESY_local_files)))
    for noe_file in NOESY_local_files:
        print('loading {} ...'.format(noe_file))
        NOESY_loc = get_dist_restraints(noe_file, s, distance_scaler_short, ramp, seq)
        s.restraints.add_selectively_active_collection(NOESY_loc, int(len(NOESY_loc) * 1.00))
   
   
    #
    # Torsion Restraints
    #

    # Load rotamer/torsion restraint files dynamically (rotamers*.dat)
    # Comment out if you don't want to use rotamer data
    rot_files = glob.glob('rotamers*.dat')
    rot_files = sorted(list(dict.fromkeys([r for r in rot_files if r])))
    if rot_files:
        for rot_file in rot_files:
            print('loading torsion restraints from {}'.format(rot_file))
            TALOS = get_torsion_restraints(rot_file, s, constant_scaler, linear_ramp_tor)
            s.restraints.add_selectively_active_collection(TALOS, int(len(TALOS) * 1.00))
    else:
        print('No rotamers.dat files found; skipping torsion restraints')


    # setup mcmc at startup
    movers = []
    n_atoms = s.n_atoms
    for i in range(0, n_res):
        n = s.index.atom(i, 'N', expected_resname=seq[i][-3:])
        ca = s.index.atom(i, 'CA', expected_resname=seq[i][-3:])
        c = s.index.atom(i, 'C', expected_resname=seq[i][-3:])
 
        atom_indxs = list(system.indexing.AtomIndex(j) for j in range(ca,n_atoms))
        mover = mc.DoubleTorsionMover(index1a=n, index1b=ca, atom_indices1=list(system.indexing.AtomIndex(i) for i in range(ca, n_atoms)),
                                      index2a=ca, index2b=c, atom_indices2=list(system.indexing.AtomIndex(j) for j in range(c, n_atoms)))

        movers.append((mover, 1))

    sched = mc.MonteCarloScheduler(movers, n_res * 60)

############################################################
    # create the options
    options = meld.RunOptions(
        timesteps = 14286,
        minimize_steps = 2000,
        min_mc = sched,
        param_mcmc_steps=200
    )


    # create a store
    store = vault.DataStore(gen_state(s,0), N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)  # why i need gen_state(s,0)? doubtful
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)

    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS, timeout=60000)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

#################################################
    return s.n_atoms


setup_system()
'''

meld_gpu_job = '''#!/bin/bash
#SBATCH --job-name={0}      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=30      
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1
#SBATCH --mem-per-cpu=2000mb            
#SBATCH --partition=gpu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --reservation=perez
#SBATCH --mem-per-cpu=2000mb          # Memory per processor
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --output=meld_{0}_%j.log     # Standard output and error log
 


source ~/.load_OpenMM
export OPENMM_CUDA_COMPILER=''

if [ -e remd.log ]; then             #If there is a remd.log we are conitnuing a killed simulation
     prepare_restart --prepare-run  #so we need to prepare_restart
fi

srun --mpi=pmix_v1  launch_remd --debug
'''
