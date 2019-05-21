# Compressed codes:

Nueva version del codigo que escribe los resultados "comprimidos" reportando únicamente el numero de partículas y el tiempo en que cambió este número.  
## many-sims-fusion-dyson-optimized.py
Version optimizada de many-sims-fusion-dyson-compressed.py en la detección de los eventos de fusión.

## many-sims-fusion-dyson-compressed.py
Corre varias simulaciones secuencialmente.  
### usage:   
many-sims-fusion-dyson-compressed.py [-h] [-s] [--freediffusion]
                                            N Niter Nsim Nsimst beta dt lf stf  

positional arguments:  
  * N             :   initial number of particles
  * Niter         :   number of iterations per simulation
  * Nsim          :   number of simulations to perform
  * Nsimst        :   number label for first simulation
  * beta          :   inverse temperature
  * dt            :   virtual time step
  * lf            :   fusion length
  * stf           :   iteration where to start the fusion processes

optional arguments:  
  * -h, --help    :   show this help message and exit
  * -s, --silent  :   do not print diagnose messages
  * --freediffusion :  turn off the log repulsion interaction  

### Resultados
Los resultados se escriben en archivos CSV de dos columnas [t, N] en donde t es el tiempo en que que hubo una fusión y N el numero de partículas que quedan después de ese tiempo.

### Filenaming convention:  
* filename=Nxx_Niterxx_Nsimxx_Nsimstxx_betaxx_dtxx_lfxx_stfxx_  
* filename+'args.json': archivo JSON con los copia de los parametros usados en las simulaciones. Este archivo se usa para leer después los parametros para el análisis por el programa "_analysis-noninteractive-decompress.py_"  
* filename+'simXX.csv': resultados de la simulación # XX

### job-array-fusion-dyson.sh  
Script para lanzar "_many-sims-fusion-dyson-compressed.py_" en un cluster con manejador de colas SLURM (magnus uniandes). Editarlo para cambiar los parámetros de la simulación. Para lanzar un _array_ de 10 jobs:  
```console
sbatch --array=0-9 job-array-fusion-dyson.sh
```
Ajustar Nsim = el numero de simulaciones por job.  

## analysis-noninteractive-decompress.py  
Este programa analiza los datos de las simulaciones: calcula el promedio y desviación estandar de N para cada tiempo. Recibe como argumento el nombre del archivo JSON con los parámetros de las simulaciones y lee todas las simulaciones disponibles en el directorio actual que correspondan a estos parámetros. Pero no tiene en cuenta Nsim ni Nsimst ya que usa todas las simulaciones disponibles que pueden provenir de varias corridas de "_many-sims-fusion-dyson-compressed.py_"  

### Resultados:
Archivo CSV con columnas [t, log_t, N_avg, std_N, log_N]
* t: tiempo
* log_t: log natural de t
* N_avg: numero promedio de partículas en el tiempo t sobre todas las simulaciones analizadas.
* std_N: desviación estandar del numero de partículas en el tiempo t
* log_N: log natural de N_avg  
Nombre del archivo: resumen_NXX_NiterXX_betaXX_dtXX_lf_XXstfXX_NsimstotXX.csv.  
Nsimstot es el numero total de simulaciones encontradas y analizadas.  

### job-analysis-noninteractive.sh
Script para correr el analysis. Editar los parámetros en el archivo para seleccionar las simulaciones a analizar.  

