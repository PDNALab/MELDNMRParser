
single_GPU_head = '''#!/bin/bash
#SBATCH --job-name={0}      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1      
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1
#SBATCH --mem-per-cpu=2000mb            
#SBATCH --partition=gpu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --reservation=perez
#SBATCH --mem-per-cpu=2000mb          # Memory per processor
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --output={0}_%j.log     # Standard output and error log
'''


meld_NMR_script = '''#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, master_runner
from meld import comm, vault
from meld import system
from meld import parse
import meld.system.montecarlo as mc
from meld.system.restraints import LinearRamp,ConstantRamp
from collections import namedtuple
import glob as glob


N_REPLICAS = 30
N_STEPS =1000
BLOCK_SIZE = 100


hydrophobes = 'AILMFPWV'
hydrophobes_res = ['ALA','ILE','LEU','MET','PHE','PRO','TRP','VAL']

def gen_state_templates(index, templates):
    n_templates = len(templates)
    print((index,n_templates,index%n_templates))
    a = system.ProteinMoleculeFromPdbFile(templates[index%n_templates])
    #Note that it does not matter which forcefield we use here to build
    #as that information is not passed on, it is used for all the same as
    #in the setup part of the script
    b = system.SystemBuilder(forcefield="ff14sbside")
    c = b.build_system_from_molecules([a])
    pos = c._coordinates
    c._box_vectors=np.array([0.,0.,0.])
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0 
    return system.SystemState(pos, vel, alpha, energy,c._box_vectors)


def get_dist_restraints(filename, s, scaler):
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
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                  r1=0.0, r2=0.0, r3=dist, r4=dist+0.2, k=350,
                                                  atom_1_res_index=i, atom_2_res_index=j,
                                                  atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists

def get_torsion_restraints(filename, s, scaler):
    torsion_rests = []
    rotamer_group = []
    for line in open(filename,'r'):
        if not line:
            torsion_rests.append(s.restraints.create_restraint_group(rotamer_group, 1))
            rotamer_group = []
        else:
            (res1, at1, res2, at2, res3, at3, res4, at4, rotamer_min, rotamer_max) = line.split()
            rotamer_max = float(rotamer_max)
            rotamer_min = float(rotamer_min)
            rotamer_avg = (rotamer_max+rotamer_min)/2.
            rotamer_sd = (rotamer_max - rotamer_min)/2.
            rotamer_rest = s.restraints.create_restraint('torsion', scaler, LinearRamp(0,100,0,1),
                                                 phi=rotamer_avg, delta_phi=rotamer_sd, k=2.5,
                                                 atom_1_res_index=int(res1), atom_1_name=at1,
                                                 atom_2_res_index=int(res2), atom_2_name=at2,
                                                 atom_3_res_index=int(res3), atom_3_name=at3,
                                                 atom_4_res_index=int(res4), atom_4_name=at4)
            rotamer_group.append(rotamer_rest)
    return torsion_rests

def setup_system():
    #
    # Start system from minimized structure(s) deposited in TEMPLATES directory
    #

    templates = glob.glob('TEMPLATES/*.pdb')
    
    #
    # build the system with force field
    #

    p = system.ProteinMoleculeFromPdbFile(templates[0])
    b = system.SystemBuilder(forcefield="ff14sbside")
    s = b.build_system_from_molecules([p])
    n_res = s.residue_numbers[-1]

    #
    # Temperature ladder
    #
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 0.4, 300., 500.)

    #
    # Scalers
    #
    distance_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=0.8, factor=4.0)
    distance_scaler_short = s.restraints.create_scaler('nonlinear', alpha_min=0.8, alpha_max=1.0, factor=4.0)
    constant_scaler = s.restraints.create_scaler('constant')

    #
    # Distance Restraints
    #
    for noe in glob.glob('NOE_*.dat'):
        print('loading {} ...'.format(noe))
        NOESY = get_dist_restraints(noe,s,scaler=distance_scaler)
        s.restraints.add_selectively_active_collection(NOESY, int( len(NOESY)*1.00 ) )
    
    if len(glob.glob('NOE_*.dat')) < 1:
        print('WARNING: no highContactOrder NOE data loaded')

    for noe in glob.glob('local_NOE_*.dat'):
        print('loading {} ...'.format(noe))
        NOESY = get_dist_restraints(noe,s,scaler=distance_scaler_short)
        s.restraints.add_selectively_active_collection(NOESY, int( len(NOESY)*1.00 ) )

    if len(glob.glob('local_NOE_*.dat')) < 1:
        print('WARNING: no shortContactOrder NOE data loaded')

    #
    # Torsion Restraints
    #
    for rotamer in glob.glob('rotamer_*.dat'):
        print('loading {} ...'.format(rotamer))
        TALOS = get_torsion_restraints(rotamer, s, constant_scaler)
        s.restraints.add_selectively_active_collection(TALOS, int( len(TALOS) * 1.00) )

    if len(glob.glob('rotamer*.dat')) < 1:
        print('WARNING: no rotamer data found')


    #
    # setup mc minimizer moves
    #
    movers= []
    n_atoms = s.n_atoms
    for i in range(1, n_res +1):
        n=s.index_of_atom(i, 'N') -1
        ca = s.index_of_atom(i, 'CA') -1
        c= s.index_of_atom(i, 'C') -1
        mover= mc.DoubleTorsionMover(n, ca, list(range(ca,n_atoms)), 
                                     ca, c, list(range(c, n_atoms)))

        movers.append((mover,1))
    sched= mc.MonteCarloScheduler(movers, n_res *60)
 
    #
    # create the simulation options
    #
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.use_big_timestep = False
    options.use_bigger_timestep = True
    options.cutoff = 1.8

    options.use_amap = False
    options.amap_alpha_bias = 1.0
    options.amap_beta_bias = 1.0
    options.timesteps = 11111
    options.minimize_steps = 20000

    #
    # MC minimizer?
    #
    options.min_mc = None
    options.run_mc = None
    #options.min_mc = sched

    #
    # Handle how to store the data
    #
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    #
    # Adaptation policy and REMD protocol
    #
    l = ladder.NearestNeighborLadder(n_trials=100)
    policy = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy)

    #
    # create and store the remd_runner
    #
    remd_runner = master_runner.MasterReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    #
    # create and store the communicator
    #
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    #
    # create and save the initial states
    #
    states = [gen_state_templates(i,templates) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    #
    # save data_store
    #
    store.save_data_store()

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
