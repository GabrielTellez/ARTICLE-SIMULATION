#!/bin/bash

#SBATCH --job-name=one-fusion-dyson		#Nombre del job
#SBATCH -p medium			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=2048		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=7-00:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH --mail-user=gtellez@uniandes.edu.co
#SBATCH --mail-type=ALL			
#SBATCH -o fusion_dyson_job.o%A_%a			#Nombre de archivo de salida

# usage: many-sims-fusion-dyson-compressed.py [-h] [-s] [--freediffusion]
#                                             [-Nf NF] [--log LOG]
#                                             N Niter Nsim Nsimst beta dt lf stf

# positional arguments:
#   N                initial number of particles
#   Niter            number of iterations per simulation
#   Nsim             number of simulations to perform
#   Nsimst           number label for first simulation
#   beta             inverse temperature
#   dt               virtual time step
#   lf               fusion length
#   stf              iteration where to start the fusion processes

# optional arguments:
#   -h, --help       show this help message and exit
#   -s, --silent     do not print diagnose messages
#   --freediffusion  turn off the log repulsion interaction
#   -Nf NF           stops the simulation when the number of particles falls
#                    below NF. Default NF=0
#   --log LOG        sets the log level. Default=INFO. DEBUG will print
#                    detailled output of the simulation


#  someter con: sbatch  job-one-diagnose.sh


N=100
Niter=100000000
Nsim=1
Nsimst=999999
beta=1.0
dt=0.0000001
lf=1.0
stf=50

module load anaconda/python3
python ./many-sims-fusion-dyson-compressed.py $N $Niter $Nsim $Nsimst $beta $dt $lf $stf --log DEBUG

