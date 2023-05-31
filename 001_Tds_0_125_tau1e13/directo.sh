#!/bin/bash
#SBATCH -N 1 # 1 nodes
#SBATCH -n 32 # 32 tasks
#SBATCH -c 1 # 1 core per task

module load Anaconda3

srun -N 1 -n 1 python main_integration.py 1 0 500 10000 &
srun -N 1 -n 1 python main_integration.py 2 500 1000 10000 &
srun -N 1 -n 1 python main_integration.py 3 1000 1500 10000 &
srun -N 1 -n 1 python main_integration.py 4 1500 2000 10000 &
srun -N 1 -n 1 python main_integration.py 5 2000 2500 10000 &
srun -N 1 -n 1 python main_integration.py 6 2500 3000 10000 &
srun -N 1 -n 1 python main_integration.py 7 3000 3500 10000 &
srun -N 1 -n 1 python main_integration.py 8 3500 4000 10000 &
srun -N 1 -n 1 python main_integration.py 9 4000 4500 10000 &
srun -N 1 -n 1 python main_integration.py 10 4500 5000 10000 &
srun -N 1 -n 1 python main_integration.py 11 5000 5500 10000 &
srun -N 1 -n 1 python main_integration.py 12 5500 6000 10000 &
srun -N 1 -n 1 python main_integration.py 13 6000 6500 10000 &
srun -N 1 -n 1 python main_integration.py 14 6500 7000 10000 &
srun -N 1 -n 1 python main_integration.py 15 7000 7500 10000 &
srun -N 1 -n 1 python main_integration.py 16 7500 8000 10000 &
srun -N 1 -n 1 python main_integration.py 17 8000 8500 10000 &
srun -N 1 -n 1 python main_integration.py 18 8500 9000 10000 &
srun -N 1 -n 1 python main_integration.py 19 9000 9500 10000 &
srun -N 1 -n 1 python main_integration.py 20 9500 10000 10000 &

wait
