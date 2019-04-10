import numpy as np
import os
import argparse


# Custom functions to use

# Function to wrap angle in order to have values between 0 and 2pi
def wrapAngle(ang):
    x=np.mod(ang,(2*np.pi))
    return x

def performfusions(thetas,fusion_angle):
    i=0
    while i < len(thetas)-1:
        no_fusion_event=True
        j=i+1
        while (j<len(thetas)) and (no_fusion_event):
            diff=np.abs(thetas[i]-thetas[j])
            if(diff>np.pi):
                diff=np.pi-diff%np.pi
            if(diff<fusion_angle):
                thetas=np.delete(thetas,(i,j))
                no_fusion_event=False
            j=j+1
        if(no_fusion_event):
            # if no fusion, increase i
            # but if fusion occurs the index i is now for the next particle, so no need to increase it
            i=i+1
    return thetas


def one_simulation(N,N_iter,dt,beta,lf,startfusion,sim):
    # ("N", type=int, help="initial number of particles")
    # ("N_iter", type=int, help="number of iterations per simulation")
    # ("dt", type=float, help="virtual time step" )
    # ("beta", type=float, help="inverse temperature")
    # ("lf",type=float, help="fusion length")
    # ("startfusion", type=int, help="iteration where to start the fusion processes")
    # ("sim", type=int, help="simulation number")
    # Calculate the mean separation depending on the initial number of particles
    sprom=2*np.pi/N
    
    # Set constant friction (f) coefficient
    # time scale is t/f
    f=1


    # Set the critical fusion angle (\theta_f)
    fusion_tol=lf*0.1*sprom    



    # Filename to save results
    patht='N_iter'+str(N_iter)+'N'+str(N)+'B'+str(beta)+'S'+str(lf)+'dt'+str(dt)+'sim'+str(sim)

    # Reset time, fusion time and number of particles vector 
    t=0
    t_file=[]
    Nt=[]
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
    # Array with iteration time number_particles
    data_n=[]
    k=0
    while k<N_iter:
        theta_later=[]
        for i in np.arange(len(theta_now)):
            sum_i=0
            # Sum over all charges to get electric force (E(\theta))
            for j in np.arange(len(theta_now)):
                if(i!=j):
                    sum_i=sum_i+1/np.tan(0.5*(theta_now[i]-theta_now[j]))
            # Calculate the mean and variance 
            mu_i=0.5*sum_i*dt/f
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
        # Evolve the system
        theta_now=theta_later   
        # Save the pertinent parameters of the time step
        Nt.append(len(theta_now))
        t_file.append(t)            
        data_n.append([t, len(theta_now)])            
        # Evolve the simulation
        t=t+dt
        k=k+1        
    # simulation end
    # Save the data on a file
    np.savetxt(patht+'.txt',data_n)


# Parse initial parameters for the simulation

parser = argparse.ArgumentParser()
parser.add_argument("N", type=int, help="initial number of particles")
parser.add_argument("N_iter", type=int, help="number of iterations per simulation")
parser.add_argument("dt", type=float, help="virtual time step" )
parser.add_argument("beta", type=float, help="inverse temperature")
parser.add_argument("lf",type=float, help="fusion length")
parser.add_argument("startfusion", type=int, help="iteration where to start the fusion processes")
parser.add_argument("N_sim", type=int, help="number of simulations to perform")
parser.add_argument("N_sim_start", type=int, help="number label for first simulation")
args=parser.parse_args()
# Initial number of particles 
N=args.N

# Number of iterations per simulation
N_iter=args.N_iter

# Set virtual time step
dt=args.dt

#simulation number
N_sim=args.N_sim


# time iteration to start fusion
startfusion=args.startfusion

lf=args.lf
beta=args.beta
N_sim_start=args.N_sim_start

for sim in range(N_sim_start,N_sim_start+N_sim):
    one_simulation(N,N_iter,dt,beta,lf,startfusion,sim)



