#!/bin/bash

#SBATCH --job-name=fusion-dyson		#Nombre del job
#SBATCH -p short			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=2000		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=24:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH --mail-user=gtellez@uniandes.edu.co
#SBATCH --mail-type=ALL			
#SBATCH -o fusion-dyson_job.o%j			#Nombre de archivo de salida

#
#usage: many-sims-fusion-dyson.py [-h]
#                                 N N_iter dt beta lf startfusion N_sim
#                                 N_sim_start
#
# positional arguments:
#  N            initial number of particles
#  N_iter       number of iterations per simulation
#  dt           virtual time step
#  beta         inverse temperature
#  lf           fusion length
#  startfusion  iteration where to start the fusion processes
#  N_sim        number of simulations to perform
#  N_sim_start  number label for first simulation
#
#
#
N=100
N_iter=2000000
dt=0.0001
beta=11
lf=0.01
startfusion=50
N_sim=50
N_sim_start=0

module load anaconda/python3
python ./many-sims-fusion-dyson.py $N $N_iter $dt $beta $lf $startfusion $N_sim $N_sim_start
