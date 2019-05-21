import numpy as np
import pandas as pd
import os
import argparse
import json
import logging


# Custom functions to use

# Function to wrap angle in order to have values between 0 and 2pi
def wrapAngle(ang):
    x=np.mod(ang,(2*np.pi))
    return x

def performfusions(thetas,fusion_angle):
    thetas=np.sort(thetas)
    i=0
    while i < len(thetas):
        no_fusion_event=True
        j=i+1
        if j == len(thetas):
            # i = N-1. Check fusion between particles 0 and N-1
            j=0
        diff=np.abs(thetas[i]-thetas[j])
        if(diff>np.pi):
            diff=np.pi-diff%np.pi
        if(diff<fusion_angle):
            thetas=np.delete(thetas,(i,j))
            no_fusion_event=False
        if(no_fusion_event):
            # if no fusion, increase i
            # but if fusion occurs the index i is now for the next particle, so no need to increase it
            i=i+1
    return thetas
def force(i,thetas):
    """ Computes the force on the i-th particle
    Parameters
    ----------
    i : int
        particle number
    thetas : ndarray
        array of particles positions
    Returns: float
        force on particle i
    """
    sum_i=0
    for j in np.arange(len(thetas)):
        if(i!=j):
            sum_i=sum_i+1/np.tan(0.5*(thetas[i]-thetas[j]))
    return sum_i

def one_simulation(N,N_iter,dt,beta,lf,startfusion,sim,filenamebase,freediffusion,N_final):
    """Performs one simulation of the Dyson brownian model on a circle with fusion events

    Parameters
    ----------
    N : int
        Initial number of particles
    N_iter : int
        Number of iterations (steps) of the simulation
    dt : float
        time step
    beta : float
        inverse temperature, coupling constant
    lf : float
        fusion length. Angle length is 0.1* lf * 2 pi/N
    startfusion : int
        iteration where to start the fusion processes
    sim : int
        simulation number used to save the data results
    filenamebase : str
        filenamebase to save the data. The data is saved to filenamebase+str(sim)
    N_final : int
        Final number of particles to stop the simulation if reached before N_iter

    """

    # friction coefficient (absorbed in the time unit)
    f=1
    # Calculate the mean separation depending on the initial number of particles
    sprom=2*np.pi/N
    # Set the critical fusion angle (\theta_f)
    fusion_tol=lf*0.1*sprom    
    # Reset time, fusion time and number of particles vector 
    t=0
    # Create vectors to evolve the simulation
    theta_now=[]
    theta_later=[]
    # Initialize the simulation with the first configuration
    # Create random noise to add to the first configuration
    ruido=np.random.uniform(-0.5*np.pi/N,0.5*np.pi/N,N)
    # Initialize charges uniformly around the circle and add noise 
    first=np.linspace(0.0,2*np.pi-2*np.pi/N,N)+ruido
    first=np.sort(first)
    theta_now=first
    # Data results dataframe iteration_time number_particles
    data_n=pd.DataFrame([[t, len(theta_now)]], columns=['t', 'N'])
    k=0 
    logging.info("Running simulation # %s", sim)
    logging.debug("time, N")
    logging.debug("%s, %s", 0, N)
    while k<N_iter and len(theta_now)>N_final:
        theta_later=[]
        for i in np.arange(len(theta_now)):
            if freediffusion:
                mu_i=0
            else:
                # Sum over all charges to get electric force (E(\theta))
                sum_i=force(i,theta_now)
                mu_i=0.5*sum_i*dt
            sigma_i=np.sqrt(2*dt/(f*beta))        
            # Generate the random step base on the latter values to evolve the system
            step=np.random.normal(mu_i,sigma_i)
            # Evolve the system and wrap the angles 
            theta_to_append=wrapAngle(theta_now[i]+step)             
            # Evolve the system
            theta_later.append(theta_to_append)
        # After the startfusion-th iteration check every step for fusion events
        if (k>=startfusion):
            theta_later=performfusions(theta_later,fusion_tol)        
       
        # Save only if there was a fusion (to save HDD space)
        if len(theta_later) < len(theta_now) :       
            data_n=data_n.append({'t': t, 'N': len(theta_later)},ignore_index=True)   
            logging.debug("%s, %s",t,len(theta_later))      
        # Evolve the simulation
        theta_now=theta_later 
        t=t+dt
        k=k+1        
    # simulation end
    logging.info("Simulation # %s completed",sim)
    # Save the data on a file
    filename=filenamebase+'sim'+str(sim)+'.csv'
    data_n.to_csv(filename, index=False)
    logging.info("Wrote results to %s",filename)

def argument_parsing():
    # Parse initial parameters for the simulation
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int, help="initial number of particles")
    parser.add_argument("Niter", type=int, help="number of iterations per simulation")
    parser.add_argument("Nsim", type=int, help="number of simulations to perform")
    parser.add_argument("Nsimst", type=int, help="number label for first simulation")
    parser.add_argument("beta", type=float, help="inverse temperature")
    parser.add_argument("dt", type=float, help="virtual time step" )
    parser.add_argument("lf",type=float, help="fusion length")
    parser.add_argument("stf", type=int, help="iteration where to start the fusion processes")
    parser.add_argument("-s", "--silent", help="do not print diagnose messages", action="store_true")
    parser.add_argument("--freediffusion", help="turn off the log repulsion interaction", action="store_true")
    parser.add_argument("-Nf", type=int, 
                        help="stops the simulation when the number of particles falls below NF. Default NF=0", 
                        default=0)
    parser.add_argument("--log", help="sets the log level. Default=INFO. DEBUG will print detailled output of the simulation",
                        default="INFO")
    return parser.parse_args()

def build_filenamebase(args):
    """ Builds the filenamebase.
    
    Paramemeters
    ------------
    args : namespace with the program arguments

    Filenaming convention: 
    filename base is filename =
    Nxx_Niterxx_Nsimxx_Nsimstxx_betaxx_dtxx_lfxx_stfxx_
    filename+'args.json': saves the simulations parameters
    filename+'simXX.csv': saves the data from simulation number XX
     
    """

    # transform arguments args to a dictionary
    args_dict=vars(args)
    # contruct the filename with the simulation parameters:

    filename=''
    for key, val in sorted(args_dict.items()):
        if key != 'silent' and key != 'freediffusion' and key != 'Nf' and key != 'log':
            filename = filename + key + str(val)+'_'
    if args.freediffusion:
        filename = filename + 'freediffusion_'
    return filename
    
def saveargs(args):
    # saves simulation parameters to a JSON file:
    args_dict=vars(args)
    paramsfilename = filenamebase + 'args.json'
    with open(paramsfilename, 'w') as paramsfile:
        json.dump(args_dict, paramsfile, indent=4, sort_keys=True)

# Main program
args = argument_parsing()
# Set up logging info
# assuming loglevel is bound to the string value obtained from the
# command line argument. Convert to upper case to allow the user to
# specify --log=DEBUG or --log=debug
loglevel=args.log
if not args.silent:
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(asctime)s : %(message)s', level=numeric_level)
filenamebase = build_filenamebase(args)
saveargs(args)
# run the simulations
logging.info("Starting %s simulations with parameters:\n %s",args.Nsim, json.dumps(vars(args), indent=4,sort_keys=True))
for sim in range(args.Nsimst,args.Nsimst+args.Nsim):
    one_simulation(args.N,args.Niter,args.dt,args.beta,args.lf,args.stf,
                   sim,filenamebase,args.freediffusion,args.Nf)
logging.info("Finished %s simulations", args.Nsim)


