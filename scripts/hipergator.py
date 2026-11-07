
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

