#!/bin/bash

#SBATCH --job-name=analysis-fusion-dyson		#Nombre del job
#SBATCH -p short			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=2048		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=47:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH --mail-user=gtellez@uniandes.edu.co
#SBATCH --mail-type=ALL			
#SBATCH -o analysis_fusion_dyson_job.o%A_%a			#Nombre de archivo de salida

N=100
Niter=2000000
Nsim=100
Nsimst=0
beta=0.01
dt=0.0001
lf=0.01
stf=50

FILENAME=N${N}_Niter${Niter}_Nsim${Nsim}_Nsimst${Nsimst}_beta${beta}_dt${dt}_lf${lf}_stf${stf}_args.json

module load anaconda/python3
python ../analysis-noninteractive-decompress.py $FILENAME




