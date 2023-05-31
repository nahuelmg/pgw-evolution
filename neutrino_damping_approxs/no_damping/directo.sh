#!/bin/bash
#SBATCH -N 1 # 1 nodes
#SBATCH -n 32 # 32 tasks
#SBATCH -c 1 # 1 core per task

module load Anaconda3

srun -N 1 -n 1 python main_integration.py 1 0 1000 10000 &
srun -N 1 -n 1 python main_integration.py 2 1000 2000 10000 &
srun -N 1 -n 1 python main_integration.py 3 2000 3000 10000 &
srun -N 1 -n 1 python main_integration.py 4 3000 4000 10000 &
srun -N 1 -n 1 python main_integration.py 5 4000 5000 10000 &
srun -N 1 -n 1 python main_integration.py 6 5000 6000 10000 &
srun -N 1 -n 1 python main_integration.py 7 6000 7000 10000 &
srun -N 1 -n 1 python main_integration.py 8 7000 8000 10000 &
srun -N 1 -n 1 python main_integration.py 9 8000 9000 10000 &
srun -N 1 -n 1 python main_integration.py 10 9000 10000 10000 &

wait
