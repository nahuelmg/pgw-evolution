#!/bin/bash
#SBATCH -N 1 # 1 nodes
#SBATCH -n 32 # 32 tasks
#SBATCH -c 1 # 1 core per task

module load Anaconda3

srun -N 1 -n 1 python main_integration.py 1 0 2500 50000 &
srun -N 1 -n 1 python main_integration.py 2 2500 5000 50000 &
srun -N 1 -n 1 python main_integration.py 3 5000 7500 50000 &
srun -N 1 -n 1 python main_integration.py 4 7500 10000 50000 &
srun -N 1 -n 1 python main_integration.py 5 10000 12500 50000 &
srun -N 1 -n 1 python main_integration.py 6 12500 15000 50000 &
srun -N 1 -n 1 python main_integration.py 7 15000 17500 50000 &
srun -N 1 -n 1 python main_integration.py 8 17500 20000 50000 &
srun -N 1 -n 1 python main_integration.py 9 20000 22500 50000 &
srun -N 1 -n 1 python main_integration.py 10 22500 25000 50000 &
srun -N 1 -n 1 python main_integration.py 11 25000 27500 50000 &
srun -N 1 -n 1 python main_integration.py 12 27500 30000 50000 &
srun -N 1 -n 1 python main_integration.py 13 30000 32500 50000 &
srun -N 1 -n 1 python main_integration.py 14 32500 35000 50000 &
srun -N 1 -n 1 python main_integration.py 15 35000 37500 50000 &
srun -N 1 -n 1 python main_integration.py 16 37500 40000 50000 &
srun -N 1 -n 1 python main_integration.py 17 40000 42500 50000 &
srun -N 1 -n 1 python main_integration.py 18 42500 45000 50000 &
srun -N 1 -n 1 python main_integration.py 19 45000 47500 50000 &
srun -N 1 -n 1 python main_integration.py 20 47500 50000 50000 &

wait
