import numpy as np
import pandas as pd
import glob
import argparse
import json
import logging

# input JSON filename to get the data 
# TODO write code to read the args and the data and create the filename 

def decompress_data(compressed_df,dt,N_iter,mul=1):
    """ Decompresses a data frame that has only the time when the particles changes 
    to one with all times for each time step.

    Parameters
    ----------
    compressed_df : pd.DataFrame
        compressed dataframe with only the changes in N and the times
    dt : float
        original time step
    N_iter : int
        Total number of iterations (time steps) in the simulation
    mul : int
        Generate data for times that are multiples of dt.
        New time step is mul*dt 
        Default mul=1

    Returns
    -------
    pd.DataFrame
        decompressed data with all the times repeating the number of particles when it does not change

    """
    tf=N_iter*dt
    steps=int(round(N_iter/m,0))
    decomp_df=pd.DataFrame([(t, 0.0) for t in np.linspace(0.0,tf-dt*m, steps)], columns=['t','N'])
    for i in compressed_df.index:
        t_now=float(compressed_df.at[i,'t'])
        if(i!=compressed_df.index[-1]):
            tf=float(compressed_df.at[i+1,'t'])
        else:
            tf=N_iter*dt
        N_now=int(compressed_df.at[i,'N'])
        decomp_df.loc[(decomp_df['t']>=t_now) & (decomp_df['t']<=tf),'N']=float(N_now)
    return decomp_df
    
def argument_parsing():
    # Parse parameters filename for the analysis
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str, help="name of the JSON file with the simulations parameters")
    parser.add_argument("-m", "--mul", type=int, help="generate data for multiples the original time step: m * dt. Default = 1", default=1)
    parser.add_argument("-s", "--silent", help="do not print diagnose messages", action="store_true")
    return parser.parse_args()

def get_params(filename):
    # Read the parameters from the JSON file
    with open(filename, "r") as read_file:
        params = json.load(read_file)
    return params

def build_filenames(params_dict):
    # filenaming convention: 
    # filenamebase is  =
    # Nxx_Niterxx_betaxx_dtxx_lfxx_stfxx_
    # to parse all data files we use glob with :
    # filenameglob = Nxx_Niterxx_Nsim*_Nsimst*_betaxx_dtxx_lfxx_stfxx_sim*.csv
    # 
    # contruct the filename with the simulation parameters:
    filenamebase=''
    filenameglob=''
    for key, val in sorted(params_dict.items()):
        if key != 'silent' and key != 'Nsim' and key != 'Nsimst' and key != 'freediffusion' and key != 'Nf' and key !='log':
            filenamebase = filenamebase + key + str(val)+'_'
            filenameglob = filenameglob + key + str(val)+'_'
        elif key == 'Nsim' or key == 'Nsimst':
            filenameglob = filenameglob + key + '*_'
        
    if params["freediffusion"]:
        filenameglob = filenameglob+'freediffusion_sim*.csv'
        filenamebase = filenamebase + 'freediffusion_'
    else:
        filenameglob = filenameglob+'sim*.csv'
    return filenamebase, filenameglob

def build_finalfilename(filenamebase, Nsims, mul=1):
    filenamefinal = 'resumen_'+filenamebase+'Nsimstot'+str(Nsims)+'_mul'+str(mul)+'.csv'
    return filenamefinal

# Main program
args = argument_parsing()
if not args.silent :
    logging.basicConfig(format='%(asctime)s : %(message)s', level=logging.INFO)
params = get_params(args.filename)
logging.info("Analysis program called with arguments:\n%s", json.dumps(vars(args), indent=4, sort_keys=True))
logging.info("Starting analysis using the parameters:\n %s", json.dumps(params, indent=4, sort_keys=True))
filenamebase, filenameglob = build_filenames(params)
dt=float(params["dt"])
N_iter=int(params["Niter"])
m=args.mul
filenames=glob.glob(filenameglob)
Nsims=len(filenames)
logging.info('Found %s simulations files to analyse', Nsims)
tf=N_iter*dt
steps=int(round(N_iter/m,0))
time_df=pd.DataFrame([(t) for t in np.linspace(0.0, tf-dt*m,steps)], columns=['t'])
# time_df=pd.DataFrame([ n*m*dt for n in np.arange(N_iter/m)], columns=['t'])
N_df=pd.concat([time_df, pd.DataFrame(np.zeros((time_df.shape[0],4)))],axis=1, ignore_index=True)
N_df.columns=['t','log_t','N_avg','std_N','log_N']
count=1
for filename in filenames:
    logging.info("Analysing file # %s : %s", count, filename)
    tmp_comp_df=pd.read_csv(filename)
    tmp_df=decompress_data(tmp_comp_df,dt,N_iter,m)
    N_df['N_avg']=N_df['N_avg']+tmp_df['N']
    N_df['std_N']=N_df['std_N']+tmp_df['N']*tmp_df['N']
    count=count+1
N_df['N_avg']=np.float64(N_df['N_avg']/Nsims)
N_df['std_N']=np.float64(N_df['std_N']/Nsims)
N_df['std_N']=np.sqrt(N_df['std_N']-N_df['N_avg']*N_df['N_avg'])
N_df['log_N']=np.log(N_df['N_avg'])
N_df['log_N']=np.log(N_df['N_avg'])
N_df['log_t']=np.log(N_df['t'])
# N_df=N_df.drop([0])
filenamefinal = build_finalfilename(filenamebase,Nsims,m)
N_df.to_csv(filenamefinal,index=False)
logging.info('Wrote final data analysis to %s', filenamefinal)
