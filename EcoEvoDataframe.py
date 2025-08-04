## ---------------------------
##
## Script name: EcoEvoDataframe.py
##
## Purpose of script: This code takes the models' direct outputs and generates dataframes for
## each batch of simulations, naming them with a 'kind' number.
## These dataframes contain several variables calculated throughout the history of the simulation.
## If you use our generated dataframes, this script is not necessary.
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
import networkx as nx
import pandas as pd
import pickle

# this is the path to the file's directory
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))



def RebuildMatrix(snaps):
# Goal: rebuild interaction matrices from the compact form stored by the simulation code
# Inputs: an array with encoded interactions at a given time
# Outputs: the matrix of interactions
    new_snaps = []
    for snap in snaps:
        rm = np.zeros((int(snap[0][0]),int(snap[0][0])))

        for el in snap[1:]:
            rm[int(el[0])][int(el[1])]=el[2]

        new_snaps.append(rm)
    return new_snaps






def NetworkSC(G):
# Goal: calculate connectance and number of species
# Inputs: network (object from networkx)
# Outputs: connectance and number of nodes

    u_G = G.to_undirected()

    if(u_G.number_of_nodes()!=0):
        Cr = 2*u_G.number_of_edges()/u_G.number_of_nodes()**2
    else:
        Cr=0
    n_nodes = G.number_of_nodes()

    
    return Cr, n_nodes


def GenerateNetworks(aa):
# Goal: generate the newtorkx network
# Inputs: matrix with all interactions (both positive and negative)
# Outputs: the network object from networkx    

    aa = np.absolute(aa) # it doesn't care about interaction type

    RN = nx.from_numpy_array(aa)
    
    RN.remove_nodes_from(list(nx.isolates(RN))) # makes sure there are no isolated nodes left

    return RN


def InteractionTypes(a,species):
# Goal: counts the numbers of each interaction types and the total density
# Inputs: matrix of interactions, vector of abundances
# Outputs: # mutualisms, # consumer-resource, # competition, # total, total density
    


    # This code separates positive and negative interactions,
    # then separates each type and counts them
    pos = np.where(a<0, 0, a)
    neg = np.absolute(np.where(a>0, 0, a))
    
    aux = np.zeros_like(pos)
    aux[pos != 0] = 1
    aux = np.logical_and(aux, aux.T)
    PP = np.where(aux == 0, pos, 0) # consumer
    MM = np.where(aux != 0, pos, 0) # mutualism
    
    bux = np.zeros_like(neg)
    bux[neg != 0] = 1
    bux = np.logical_and(bux, bux.T)
    pp = np.where(bux == 0, neg, 0) # resource
    CC = np.where(bux != 0, neg, 0) # competition
    
    # counts
    countM = np.count_nonzero(MM)/2 # each two entries are 1 mutualism interaction
    countP = np.count_nonzero(PP) # each consumer (or resource) entry is 1 consumer-resource interaction
    countC = np.count_nonzero(CC)/2 # each two entries are 1 competition interaction
    countT = countM+countP+countC
    

    
    total_density = np.sum(species)
    
    
           
    return countM,countP,countC,countT, total_density






def CalculateHistory(net,tvec):
# Goal: call InteractionTypes for all timesteps in the series and store the results
# Inputs: list with interaction matrix and vector of abundances, list of times
# Outputs: lists with values for all times for: original times, mutualism counts, consumer-resource counts, competition counts, total counts, total densities

    m_snaps,y_snaps,_ = net
    times = np.arange(len(y_snaps))
    
    
    M = []
    P = []
    C = []
    T = []

    total_d = []
    
    for t in range(len(tvec)):
        
        

        countM,countP,countC,countT, total_density = InteractionTypes(m_snaps[t],y_snaps[t])

        
        
        M.append(countM)
        P.append(countP)
        C.append(countC)
        T.append(countT)
        
        total_d.append(total_density)
            
        
    return times, M, P, C, T, total_d




def SaveDF():
    
    df = pd.DataFrame() # initialize the df
    evo_history=[] # this will store history of variables
    evo_historyX=[] # this stores the history to be saved together with the df
    rk_list=[] # this will store the list of additional variables defining the model
    tvec = list(np.array([1,50,100,150,200,250,300,350,400,450,499])) # this list chooses the assembly times in which we calculate values (here chosen at each 50 assembly events)

    # helper code to read the files in the particular way we stored
    if(istest==0): # isTest=0 when reproducing the paper simulations
        if(isHeatmap==0): # when generating the 8 scenarios (figs 2, 3, 4), this chooses 20 samples
            frange = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
        
        else: # when generating the heatmaps, this chooses the 3 batches of 5 samples
            if(kind%10==0):
                frange = [1,2,3,4,5]
            if(kind%10==1):
                frange = [6,7,8,9,10]
            if(kind%10==2):
                frange = [11,12,13,14,15]
                
    if(istest==1): # isTest=1 when generating custom simulations
        frange = [1,2] # choose the number of samples to include in the df
        

    for j in frange: # These will run for all sample codes, the same as the num_sample in the Model code
        path = dir_path+'/data/'+folder+str(j)+'FEVO-O1.pkl'
        with open(path, 'rb') as file:
          

            sollist=[]
            sollist.append(pickle.load(file)) # receive the loaded file from the model
            
        for i in range(len(sollist[0])):

            # Storing variables
            numbers=sollist[0][i][1]
            rk_list.append(sollist[0][i][2])
            mat = RebuildMatrix(numbers[2]) # the compacted matrix is rebuilt
            numbers[2]=mat
            
            
            
              

       
            # Calculate the history of variables for all times in tvec
            evo_history.append([[numbers[2][ii] for ii in tvec],[numbers[3][ii] for ii in tvec],[]])
            times, M, P, C, T, total_d = CalculateHistory(evo_history[-1], tvec)
            
            
            dic = {} # Initialise the sample dic that'll entry the df

            
            
       
        
            # This block will compose the dic variables for each time in tvec
            # The final dic will have all variables for all chosen times, like a timeseries
            k=0
            for ss in range(len(tvec)):
                
                # calculate connectance and richness
                RN = GenerateNetworks(evo_history[-1][0][ss])
                Cr, n_nodes = NetworkSC(RN)
                
                
                dic['kind']=kind # store kind of simulations


                # store total density, richness, and connectance
                dic['tot_density'+str(tvec[ss])]=total_d[k]
                dic['n_species'+str(tvec[ss])]=n_nodes
                try:
                    dic['Con'+str(tvec[ss])] = 2*RN.number_of_edges()/n_nodes**2
                except:
                    dic['Con'+str(tvec[ss])] = 0
                
                
                 
                
                # store all proportions of interaction types
                jj=0
                for nn in [M,P,C,T]:

                    name = ['M','P','C','T']
                    
                    
                    dic[name[jj]+str(tvec[ss])] = nn[k]
                    
                    if(jj<3):
                        try:
                            dic['prop'+name[jj]+str(tvec[ss])] = nn[k]/T[k]
                        except:
                            dic['prop'+name[jj]+str(tvec[ss])] = 0
            
                    jj+=1
                
                
                
                k+=1
            
            

            
            # include the dic into the df
            df = pd.concat([df,pd.DataFrame([dic])], ignore_index = True)
            
            print(j) # log the processing of a sample
            
            
            evo_historyX.append(evo_history[-1])

       

                               
            
    # save the df with additional model info
    path = dir_path+'/df/'+saveFolder+str(kind)+'-df-0.pkl'
    with open(path, 'wb') as file:
          
        pickle.dump([df,evo_historyX,rk_list], file)










##################################################

isHeatmap=0 # This value is 1 for heatmaps and 0 for other graphics
istest=1 # This value is 1 for custom simulations and 0 for the simulations in the figures

if(isHeatmap==0 and istest==0): # This generates the df's for the 8 scenarios (figs 2, 3, 4)
    communityType = ['1','2']
    modelType = ['evo5','evo30','invR','invP']
    for i in [0,1]: # for types 1 and 2
        for j in [0,1,2,3]: # kinds for evo5, evo30, invR, and invP
            
            kind = j # This labels the type for all simulations in the dataframe that will be saved, for classification purposes among dfs
            folder = 'modelData_1/t'+communityType[i]+'_'+modelType[j]+'_1/' # specific folder where the simulations are stored
    
            saveFolder='t'+str(i+1)+'_df_1/'
    
            print(kind,folder)
            
            SaveDF()

if(isHeatmap==1 and istest==0): # This generates the df's for the heatmaps
    levelDic = {0.01:'',0.05:'1',0.1:'2',0.15:'3',0.2:'4',0.25:'5'}
    levelDic2 = {0.01:'6',0.05:'7',0.1:'8',0.15:'9',0.2:'10',0.25:'11'}
    for heatBatch in ['','2','3']:
        for level in ['','2']:
            for sl in [0.01,0.05,0.1,0.15,0.2,0.25]:
                for j in ['0','1','2']:
                    
                    if(level==''):
                        dec=levelDic[sl]
                    if(level=='2'):
                        dec=levelDic2[sl]
            
                    folder = 'heatData_1/'+heatBatch+'heat'+level+'_'+str(sl)+'/' # specific folder where the simulations are stored
                    kind=int( heatBatch+heatBatch+dec+j ) # This labels the type for all simulations in the dataframe that will be saved, for classification purposes among dfs
                    saveFolder='heat_df_1/'
                    
                    print(kind,folder)
                    
                    SaveDF()
                  
                

if(istest==1): # This generates the df's for custom simulations
    
            
    kind = 0 # This labels the type for all simulations in the dataframe that will be saved, for classification purposes among dfs
    folder = 'testData_1/' # specific folder where the simulations are stored, change it accordingly

    saveFolder='test_df_1/' # specific folder where simulations will be saved, change it accordingly

    print(kind,folder)
    
    SaveDF()