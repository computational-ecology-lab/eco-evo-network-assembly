## ---------------------------
##
## Script name: EcoEvoMicData.py
##
## Purpose of script: This code loads the dataframes of simulations, analyses the microbial timeseries data,
## generates stochastic simulations from the df's, then plots the MAD and correlation abundance patterns for both
## simulations and data. The code runs separately for each scenario (each df)
##
## Author: gui araujo
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
##
## Date Created: 2024
##
## Email: gui.dav.araujo@gmail.com
##
## ---------------------------


# Packages used
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
np.random.seed(555)

# this is the path to the file's directory
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))


# Set the parameters for figure
fsize=14
params = { 'figure.dpi': 600,
        'legend.fontsize': fsize-3,
          #'figure.figsize': (3.54, 3.54),
         'axes.labelsize': fsize,
         'axes.titlesize':fsize,
         'xtick.labelsize':fsize-2,
         'ytick.labelsize':fsize-2}
mpl.rcParams.update(params)


def NDynamics(y0,time):
# Goal: generate the stochastic solutions from the final point of simulations,
# adding the external noise.
# Inputs: final state of the model and total time
# Outputs: densities array and time array


    TT = np.linspace(0,time,100*time+1) # Time array. We integrate the function for these time values


    sol = [y0]
    
    # This function reads and modifies global parameters that are shared between functions
    # These will hold the submatrices for each type of interaction
    global MM # mutualism
    global PP # consumer
    global pp # resource
    global CC # competition
    
    # Creates positive types (MM,pp) from the matrix of positive interactions
    aux = np.zeros_like(pos)
    aux[pos != 0] = 1
    aux = np.logical_and(aux, aux.T)
    PP = np.where(aux == 0, pos, 0)
    MM = np.where(aux != 0, pos, 0)
    
    # This is for negative types
    bux = np.zeros_like(neg)
    bux[neg != 0] = 1
    bux = np.logical_and(bux, bux.T)
    pp = np.where(bux == 0, neg, 0)
    CC = np.where(bux != 0, neg, 0)
    
    # These will hold the denominators for type II
    global p_load # consumer-resource
    global mp_load # mutualism
    
    
    
    global cl # harvesting cost function (sum of lambdas for positive interactions)
    
    if(hasCost==1):
        # Calculate cost from positive types
        cl = cost * ( np.count_nonzero(PP,axis=1)+np.count_nonzero(MM,axis=1) )
    
    # Loop for each step of the solution
    for iii in range(len(TT)-1):
    
        # Euler-Maruyama update for the SDE
        y0s = y0.copy()
        y0s = np.array(y0) + (TT[1]-TT[0])*Sys(y0) + np.array(y0)*(sig*np.sqrt(TT[1]-TT[0])*np.random.randn(len(y0)))
          
        
        y0s[y0s<=exth]=0
        y0=y0s.copy()
        sol.append(y0)
        #test.append(y0)
    return np.array(sol).T, TT




def Sys(yy):
# Goal: to be called from NDynamics to update the SDE integration, this is the model's equation
# Inputs: solutions list
# Outputs: the differentials

    # The equation uses several global parameters defined by GenerateTimeSeries and NDynamics
    # Calculate type II denominators from positive types
    p_load = 1/(1 + hp*np.dot(PP,yy))      
    mp_load = 1/(1+hm*np.dot(MM, yy))
    
    
    dy = yy*np.dot(PP,yy)*p_load - yy*np.dot(pp,np.multiply(yy,p_load)) - yy*np.dot(CC,yy) + yy*np.dot(MM,yy)*mp_load + (r-cl)*yy  -yy*np.dot(ins,yy)  
 
    return dy





def GenerateTimeSeries(model,rk=[]):
# Goal: generate timeseries from the stochastic solutions, to later be sampled
# Inputs: a list of simulations (model), the std of the white noise,
# a list of model variables and parameters (rk)
    
    # This function reads and modifies global parameters to use between functions
    global S # richness
    global pos # matrix of positive interactions
    global neg # matrix of negative interactions
    global ins # matrix of intraspecific competition values
    global r # vector of individual growth rates
    global y0 # vector of densities
    global sig # standard deviation of the noise
    
    solList=[] # list of solutions
    iii=0
    for net in model: # loop for each simulation in the list

        tt=-1 # this will take only the last step of simulations
        
        # load values from the simulation list
        pos=np.where(net[0][tt]>0,net[0][tt],0)
        neg=np.where(net[0][tt]<0,-net[0][tt],0)
        y0=net[1][tt]
        ins=rk[iii][0][1]
        r=rk[iii][0][0]
        S=np.array(r).shape[0]
        y0=rk[iii][0][2]
        pos=rk[iii][0][3]
        neg=rk[iii][0][4]
        S=y0.shape[0]
        #S=0
        
        poss = np.zeros((S,S))
        negg = np.zeros((S,S))
        y00 = np.zeros(S)
        inss = np.zeros((S,S))
        rr = np.zeros(S)
        for ii in range(S):
            for jj in range(S):
                poss[ii][jj]=pos[ii][jj]
                negg[ii][jj]=neg[ii][jj]
        for ii in range(S):
            y00[ii] = y0[ii]
            rr[ii] = r[ii]
            inss[ii][ii] = ins[ii][ii]
        pos = poss.copy()
        neg = negg.copy()
        y0 = y00.copy()
        r = rr.copy()
        ins = inss.copy()
            
        
        
        # total simulation time
        maxtime=1000#0#20
        
        print('-') # log the start of the loop
        
        
        # Call the model
        sol,solT = NDynamics(y0.reshape(y0.shape[0],),maxtime)
        
        
        print('oo') # log the end of the model's run
        
        sol_ = [] # list of points in the timeseries
        cutoff = 0#int(50*maxtime)#int(25*maxtime) # determine a starting point for the time-series
        
        for ii in range(sol.shape[0]): # loop for each species

            if(np.mean(sol[ii][cutoff:])>0): # if the species didn't become extinct
                ab = sol[ii]
                sol_.append(ab[cutoff::10])
        solList.append(sol_)
        
        iii+=1
                
    print('!!!!!') # log the end
    
    return solList




# global parameters
r=0
p_load=0
mp_load=0
PP=0
MM=0
pp=0
CC=0
pos=0
neg=0
ins=0
S=0
y0=0
cl = 0

# 1 to enable interaction cost (delta). this is always on
hasCost=1

# parameters
exth = 0.000001 # extinction threshold
hm=0.1 # handling time for mutualism
hp=0.1 # handling time for consumer-resource
sig=0.1#.05 # standard deviation of noise
n_counts=20000
#test=[]


# This block loads the dataframes with model simulations (generated by Dataframe code)
########################################################
evo_history0=[] # This list will have only one element here. This code loads one df at a time (one kind).


############################### TO CHANGE ######################################################
kind=2 # Change this according to the name of the df to load: 0 - Evo 5, 1 - Evo 30, 2 - Inv Rand, 3 - Inv Prop
n_nets=20 # Change this according to the number of simulations in the df to load (it's 20 in the paper)
model_type=2  # Change this according to type 1 or type 2
isTest=0 # Change to 1 if running custom simulations
################################################################################################


# If using custom costs, change this manually to the cost used
cost = 0.01
if(model_type==1):
    cost=2.5*cost
    
    
    
filecodes=[[kind,0]] # Don't change. Second number is not being used anymore

if(isTest==0):
    path = dir_path+'/df/t'+str(model_type)+'_df_1/'+str(filecodes[0][0])+'-df-'+str(filecodes[0][1])+'.pkl'
if(isTest==1):
    path = dir_path+'/test_df_1/'+str(filecodes[0][0])+'-df-'+str(filecodes[0][1])+'.pkl'

with open(path, 'rb') as file:
    sollist=pickle.load(file)
    evo_history0.append([])
    evo_history0[-1].extend(sollist[1])
    evo_history0[0] = evo_history0[0][:n_nets]
ts = GenerateTimeSeries(evo_history0[0],rk=sollist[2])

# Depending on the kind, this chooses the colors of graphics (evo5 is 0 - blue, evo30 is 1 - green, invR is 2 - orange, invP is 3 - red)
if(kind==2):
    for i in range(n_nets):
        col='#ff7f00'
if(kind==1):
    for i in range(n_nets):
        col='#4daf4a'
if(kind==0):
    for i in range(n_nets):
        col='#377eb8'
if(kind==3):
    for i in range(n_nets):
        col='#f781bf'
        

Ldf=[]
for i in range(n_nets):
    net = ts[i]
    Ldf.append(pd.DataFrame(np.array(net).T))
    

print('MODEL PROCESSED')


# This is the data-processing code. It reads and transforms the data
####################################################################
import pandas as pd


# Load the data
path = dir_path+'/timeseries_data/timeseries.csv'
df = pd.read_csv(path)

# Drop columns that won't be used
df.drop(columns=['Unnamed: 0','project_id','sample_id','run_id','classification','host_age','date_collection'],inplace=True)



# This block just reshapes the data to produce timeseries tables of each type of sample
df['experiment_day'] = df['experiment_day'].astype(str)
df['sampleday'] = df['host_id']+df['experiment_day']+df['samplesite']
df['sample'] = df['host_id']+df['samplesite']
df.drop(columns=['samplesite','host_id'],inplace=True)
df_info = df[['sampleday','sample', 'experiment_day']].drop_duplicates()
df = df.pivot_table(index='sampleday',columns='otu_id', values='count', aggfunc='sum').reset_index()
df = df.fillna(0)
df = pd.merge(df_info,df,on='sampleday')


# This drops samples with less than 10k reads
df = df.drop(df[df.iloc[:,3:].sum(axis=1)<10000].index)


# Produce table of relative abundances
df_r = df.iloc[:,3:].div(df.iloc[:,3:].sum(axis=1),axis=0)
df_r['sample'] = df['sample']
df_r['experiment_day'] = df['experiment_day']
sample_counts = df_r['sample'].value_counts()


# These create a different df for each type of data
# Some of them will be empty or almost empty, because of the filter on reads, so they'll not be used
df_mtongue = df_r[df['sample']=='M3Tongue']
df_mtongue=df_mtongue.loc[:, (df_mtongue != 0).any(axis=0)]

df_mlpalm = df_r[df['sample']=='M3L_palm']
df_mlpalm=df_mlpalm.loc[:, (df_mlpalm != 0).any(axis=0)]

df_mrpalm = df_r[df['sample']=='M3R_palm']
df_mrpalm=df_mrpalm.loc[:, (df_mrpalm != 0).any(axis=0)]

df_mfeces = df_r[df['sample']=='M3feces']
df_mfeces=df_mfeces.loc[:, (df_mfeces != 0).any(axis=0)]

df_ftongue = df_r[df['sample']=='F4Tongue']
df_ftongue=df_ftongue.loc[:, (df_ftongue != 0).any(axis=0)]

df_flpalm = df_r[df['sample']=='F4L_palm']
df_flpalm=df_flpalm.loc[:, (df_flpalm != 0).any(axis=0)]

df_frpalm = df_r[df['sample']=='F4R_palm']
df_frpalm=df_frpalm.loc[:, (df_frpalm != 0).any(axis=0)]

df_ffeces = df_r[df['sample']=='F4feces']
df_ffeces=df_ffeces.loc[:, (df_ffeces != 0).any(axis=0)]



############################### Pattern: MAD
i=0
plt.figure(figsize=(4,4))
Ldf_net = Ldf.copy()
binsList = []
nList = []
#for df_i in [df_mtongue,df_mrpalm,df_mfeces,df_ffeces,*Ldf_net]:
for df_i in [df_mfeces,df_ffeces,df_mtongue,df_mrpalm,df_flpalm,df_mlpalm,*Ldf_net]:
    if(i<3+3):
        mean = df_i.iloc[:,:-2].astype('float')
        mean = mean.mean(axis=0)
        mean = (np.log(mean)-np.log(mean).mean())/np.log(mean).std()
        #mean=mean[mean>0]
        n, bins, patches = plt.hist(mean, bins=15, density=True)
        binsList.append(bins)
        nList.append(n)
    else:
        mean = df_i.astype('float')
        mean=mean[mean>0]
        mean = mean.mean(axis=0)
        n, bins, patches = plt.hist((np.log(mean)-np.log(mean).mean())/np.log(mean).std(), bins=15,density=True)
        binsList.append(bins)
        nList.append(n)
    i+=1
plt.show()


#x = np.linspace(-5,5,1000000)
#y = np.exp(-x**2/2)/np.sqrt(2*np.pi)
#plt.plot(x,y,color='gray')
i=0
plt.figure(figsize=(4,4))
for df_i in [df_mfeces,df_ffeces,df_mtongue,df_mrpalm,df_flpalm,df_mlpalm,*Ldf_net]:
    if(i<3+3):


        n = nList[i]
        bins = binsList[i]
        bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
        plt.scatter(bins_mean, n,alpha=0.5,c='k')
        plt.yscale('log')
        plt.ylim(0.005,1)
        #plt.xlim(-7,7)
        plt.ylabel('Probability density')
        plt.xlabel('Rescaled log average relative abundance')
    else:
        n = nList[i]
        bins = binsList[i]
        bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
        plt.scatter(bins_mean, n,s=10,c=col,alpha=0.5)
        plt.yscale('log')
        plt.ylim(0.005,1)
        #plt.xlim(-7,7)
        plt.ylabel('Probability density')
        plt.xlabel('Rescaled log average rel. abundances')


    i+=1

plt.savefig('MAD.png', dpi=300, bbox_inches='tight')
plt.show()





############################### Pattern: correlations


d4_binsList = []
d4_nList = []
plt.figure(figsize=(4,4))
for df_i in [df_mfeces,df_ffeces,df_mtongue,df_mrpalm,df_flpalm,df_mlpalm]:
            df_i = df_i.iloc[:,:-2].astype('float')
            pears = df_i.corr()
            pears=pears.stack()
            pears = pears[pears<1]
            
            n, bins, patches = plt.hist(pears, bins=15, density=True)
            d4_binsList.append(bins)
            d4_nList.append(n)
plt.show()
            
 
Ldf_net = Ldf.copy()
binsList=[]
nList=[]
plt.figure(figsize=(4,4))
for df_i in [*Ldf_net]:
        df_i = df_i.astype('float')
        pears = df_i.corr()
        pears=pears.stack()
        pears = pears[pears<1]
        
        n, bins, patches = plt.hist(pears, bins=15, density=True)
        binsList.append(bins)
        nList.append(n)

        
plt.show()
        
i=0
plt.figure(figsize=(4,4))
for df_i in [df_mfeces,df_ffeces,df_mtongue,df_mrpalm,df_flpalm,df_mlpalm,*Ldf_net]:
    if(i<3+3):

        n = d4_nList[i]
        bins = d4_binsList[i]
        bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
        plt.scatter(bins_mean, n,alpha=0.5,c='k')
        plt.yscale('log')
        plt.ylim(1e-7,1e2)
        plt.vlines(x=0,linestyles='--',ymin=1e-7,ymax=1e2,color='gray')
        plt.ylabel('Probability density')
        plt.xlabel('Correlation coefficients of relative abundances')
    else:
        n = nList[i-3-3]
        bins = binsList[i-3-3]
        bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(n))]
        plt.scatter(bins_mean, n,s=10,c=col,alpha=0.5)
        plt.yscale('log')
        plt.ylim(1e-7,1e2)
        plt.vlines(x=0,linestyles='--',ymin=1e-7,ymax=1e2,color='gray')
        plt.ylabel('Probability density')
        plt.xlabel('Correlation coef. of rel. abundances')

    i+=1

plt.savefig('CORR.png', dpi=300, bbox_inches='tight')
plt.show()




