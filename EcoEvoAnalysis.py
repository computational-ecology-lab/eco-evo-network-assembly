## ---------------------------
##
## Script name: EcoEvoAnalysis.py
##
## Purpose of script: This code loads the dataframes of simulations and creates the figures from the paper (apart from the data analysis).
## after dataframes are produced by EcoEvoDataframe, they can be loaded in this code to be processed for figures.
## there are two modes: for heatmaps, and for other figures. For heatmaps, there should be a dataframe for each
## cell of the heatmap. For other figures, there should be a dataframe for each type of model (each color).
## Each dataframe should contain a collection of simulations (samples).
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
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import math
import seaborn as sns
import community
import matplotlib as mpl

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


def Modularity(Gr):
# Goal: calculate unweighted modularity of network
# Inputs: network (object from networkx)
# Outputs: modularity value

    G = Gr.to_undirected()
    part = community.best_partition(G)
    mod = community.modularity(part, G)
    return mod


def Entropy(G,degrees):
# Goal: calculate degree entropy of network
# Inputs: network (object from networkx), list of degrees
# Outputs: entropy value
    
    vk = degrees
    try:
        maxk = np.max(vk)
    
        #mink = np.min(vk)
        kvalues= np.arange(0,maxk+1) # possible values of k
        Pk = np.zeros(maxk+1) # P(k)
        for k in vk:
            Pk[k] = Pk[k] + 1
        Pk = Pk/sum(Pk) # the sum of the elements of P(k) must to be equal to one
        
        H = 0
        for p in Pk:
            if(p > 0):
                H = H - p*math.log(p, 2)
    except:
        H=0
    
    return H


def NetworkMetrics(G,degrees):
# Goal: calculate connectance, number of nodes (richness), and entropy (using the Entropy function)
# Inputs: network (object from networkx), list of degrees
# Outputs: values of connectance, number of nodes, and entropy

    u_G = G.to_undirected()
    
    if(u_G.number_of_nodes()!=0):
        Cr = 2*u_G.number_of_edges()/u_G.number_of_nodes()**2 # C = L / (S^2/2)
    else:
        Cr=0
    n_nodes = G.number_of_nodes()
    

    Hr = Entropy(u_G,degrees)

    
    return Cr, n_nodes, Hr


def GenerateNetworks(aa):
# Goal: generate the newtorkx network
# Inputs: matrix with all interactions (both positive and negative)
# Outputs: the network object from networkx    

    aa = np.absolute(aa) # it doesn't care about interaction type

    RN = nx.from_numpy_array(aa)
    
    RN.remove_nodes_from(list(nx.isolates(RN))) # makes sure there are no isolated nodes left

    return RN









isHeatmap=0 # Choose 1 if this run is for heatmaps or 0 if for 8 scenarios
isTest=0 # Choose 1 to run df's for custom simulations and 0 for simulations from the paper

communityType = 2 # If the run is for the 8 scenarios (isHeatmap=isTest=0), choose if it's type 1 or type 2

tvec = [1,50,100,150,200,250,300,350,400,450,499] # this list chooses the assembly times in which we calculate values (has to match the list used in EcoEvoDataframe)


if(isTest==1): # Parameters for custom simulations
    folder = 'analysis_df_1/' # Choose this folder as the one where you stored the df's
    filecodes = [[0,0],[1,0],[2,0],[3,0]] # The second value is always 0, choose the first to match the kind of your df's and the quantity of entries 
    num_samples=20 # Number of samples in each df (should be the all same, but this is only for calculating standard errors)

if(isHeatmap==0 and isTest==0):

    # This lists the types of simulations to load (numbers on the left; leave the ones on the right as 0)
    # 0: evo5, 1: evo30, 2: invR, 3: invP
    filecodes = [[0,0],[1,0],[2,0],[3,0]]
    num_samples=20 # Number of samples in each df, for standard errors
    
    
if(isHeatmap==1 and isTest==0): # The heatmaps have a sort of convoluted way to sort out the file numbers
    
    filecodes = [[0,0],[1,0],[2,0],[10,0],[11,0],[12,0],[20,0],[21,0],[22,0],[30,0],[31,0],[32,0],[40,0],[41,0],[42,0],
                 [50,0],[51,0],[52,0],[60,0],[61,0],[62,0],[70,0],[71,0],[72,0],[80,0],[81,0],[82,0],[90,0],[91,0],[92,0],
                 [100,0],[101,0],[102,0],[110,0],[111,0],[112,0],
                 [220,0],[221,0],[222,0],[2210,0],[2211,0],[2212,0],[2220,0],[2221,0],[2222,0],[2230,0],[2231,0],[2232,0],[2240,0],[2241,0],[2242,0],
                 [2250,0],[2251,0],[2252,0],[2260,0],[2261,0],[2262,0],[2270,0],[2271,0],[2272,0],[2280,0],[2281,0],[2282,0],[2290,0],[2291,0],[2292,0],
                 [22100,0],[22101,0],[22102,0],[22110,0],[22111,0],[22112,0],
                 [330,0],[331,0],[332,0],[3310,0],[3311,0],[3312,0],[3320,0],[3321,0],[3322,0],[3330,0],[3331,0],[3332,0],[3340,0],[3341,0],[3342,0],
                              [3350,0],[3351,0],[3352,0],[3360,0],[3361,0],[3362,0],[3370,0],[3371,0],[3372,0],[3380,0],[3381,0],[3382,0],[3390,0],[3391,0],[3392,0],
                              [33100,0],[33101,0],[33102,0],[33110,0],[33111,0],[33112,0]]



# Read the dataframes listed in filecodes
df = pd.DataFrame()
evo_history = []
if(isHeatmap==1 and isTest==0):
    folder = 'heat'
elif(isHeatmap==0 and isTest==0):
    folder = 't'+str(communityType)
for i in filecodes:
    path=dir_path+'/df/'+folder+'_df_1/'+str(i[0])+'-df-'+str(i[1])+'.pkl'
    with open(path, 'rb') as file:
          
        sollist=pickle.load(file)
    
        df = pd.concat([df,sollist[0]], ignore_index = True)
        evo_history.extend(sollist[1])
        
        



hdf=df.copy()
if(isHeatmap==1 and isTest==0):                
    # This block is to adjust the dataframe with the correct sigmas and costs    
    hdf=hdf.iloc[:,101:] # This should be calibrated to get only the variables of the last recorded timestep (we use 499)... to do so, look in the dataframes for the right position
    sigma = [0.01,0.05,0.1,0.15,0.2,0.25]
    cost0=100*0.01
    samples=5
    hdf['Complexity499']=hdf['n_species499']*hdf['Con499'] # Include complexity in the df
    
    # This block includes sigmas and costs for all combinations
    cost = np.round([0.5*cost0,cost0,1.5*cost0],2)
    for s in range(len(sigma)):
        for c in range(len(cost)):
            for r in range(samples):
                hdf.loc[s*len(cost)*samples+c*samples+r,'sigma']=sigma[s]
                hdf.loc[s*len(cost)*samples+c*samples+r,'cost']=cost[c]
    cost = np.round([2*cost0,2.5*cost0,3*cost0],2)
    for s in range(len(sigma)):
        for c in range(len(cost)):
            for r in range(samples):
                hdf.loc[len(sigma)*len(cost)*samples+ s*len(cost)*samples+c*samples+r,'sigma']=sigma[s]
                hdf.loc[len(sigma)*len(cost)*samples+ s*len(cost)*samples+c*samples+r,'cost']=cost[c]
    
    cost = np.round([0.5*cost0,cost0,1.5*cost0],2)
    for s in range(len(sigma)):
        for c in range(len(cost)):
            for r in range(samples):
                hdf.loc[2*len(sigma)*len(cost)*samples+s*len(cost)*samples+c*samples+r,'sigma']=sigma[s]
                hdf.loc[2*len(sigma)*len(cost)*samples+s*len(cost)*samples+c*samples+r,'cost']=cost[c]
    cost = np.round([2*cost0,2.5*cost0,3*cost0],2)
    for s in range(len(sigma)):
        for c in range(len(cost)):
            for r in range(samples):
                hdf.loc[2*len(sigma)*len(cost)*samples+len(sigma)*len(cost)*samples+ s*len(cost)*samples+c*samples+r,'sigma']=sigma[s]
                hdf.loc[2*len(sigma)*len(cost)*samples+len(sigma)*len(cost)*samples+ s*len(cost)*samples+c*samples+r,'cost']=cost[c]
    
    cost = np.round([0.5*cost0,cost0,1.5*cost0],2)
    for s in range(len(sigma)):
        for c in range(len(cost)):
            for r in range(samples):
                hdf.loc[4*len(sigma)*len(cost)*samples+s*len(cost)*samples+c*samples+r,'sigma']=sigma[s]
                hdf.loc[4*len(sigma)*len(cost)*samples+s*len(cost)*samples+c*samples+r,'cost']=cost[c]
    cost = np.round([2*cost0,2.5*cost0,3*cost0],2)
    for s in range(len(sigma)):
        for c in range(len(cost)):
            for r in range(samples):
                hdf.loc[4*len(sigma)*len(cost)*samples+len(sigma)*len(cost)*samples+ s*len(cost)*samples+c*samples+r,'sigma']=sigma[s]
                hdf.loc[4*len(sigma)*len(cost)*samples+len(sigma)*len(cost)*samples+ s*len(cost)*samples+c*samples+r,'cost']=cost[c]

    
    

    
    
    ####### HEATMAPS - Plots
    mean_df = hdf.groupby(['sigma', 'cost']).mean().reset_index()
    
    value = 'propM'
    sns.heatmap(mean_df.pivot('sigma','cost',str(value)+'499')).invert_yaxis()
    plt.title('Proportion of mutualism')
    plt.show()
    value = 'propP'
    sns.heatmap(mean_df.pivot('sigma','cost',str(value)+'499')).invert_yaxis()
    plt.title('Proportion of consumer-resource')
    plt.show()
    value = 'propC'
    sns.heatmap(mean_df.pivot('sigma','cost',str(value)+'499')).invert_yaxis()
    plt.title('Proportion of competition')
    plt.show()
    value = 'Con'
    sns.heatmap(mean_df.pivot('sigma','cost',str(value)+'499')).invert_yaxis()
    plt.title('Connectance')
    plt.show()
    value = 'n_species'
    sns.heatmap(mean_df.pivot('sigma','cost',str(value)+'499')).invert_yaxis()
    plt.title('Species richness')
    plt.show()
    value = 'Complexity'
    sns.heatmap(mean_df.pivot('sigma','cost',str(value)+'499')).invert_yaxis()
    plt.title('Complexity (richness x connectance)')
    plt.show()



if(isHeatmap==0): # This block generates plots for scenarios besides the heatmaps
    
    if(isTest==0): # For the 8 scenarios, plots in figures 3, 4, S1, S2
        # Set the names of networks based on kinds and their respective colors
        dicKind = {0:'Evo. ($\Delta=5$)', 1:'Evo. ($\Delta=30$)',2:'Inv. Rand.',3:'Inv. Prop.' } # Scenario names
        colors=['#377eb8','#4daf4a','#ff7f00','#f781bf'] # Colours used
    
    if(isTest==1): # This is for custom simulations: choose the names and colours accordingly
        #dicKind = {0:'Evo. ($\Delta=5$)', 1:'Evo. ($\Delta=30$)',2:'Inv. Rand.',3:'Inv. Prop.' }
        dicKind = {0:'T1 Evo. ($\Delta=5$)', 1:'T2 Evo. ($\Delta=5$)',2:'K-T1 Evo. ($\Delta=5$)',3:'K-T2 Evo. ($\Delta=5$)' }
        colors=['#377eb8','#4daf4a','#ff7f00','#f781bf']
    
    
    
    # This generates the graphics of proportions of interaction types
    # The codes don't work properly if there's just one sample in a dataframe,
    # because they need to calculate variances
    ############################# PROPORTIONS ##################################
    #########################################################################
    df_types = []
    for k in df['kind'].unique():
        df_types.append(df[df['kind']==k].copy())
    
    
    for name in ['propM','propP','propC']:
        parl = []
        parls=[]
        times = tvec
        plt.figure(figsize=(4,4))
        for k in df['kind'].unique():
            parl.append([])
            parls.append([])
        for i in times:
            par = name+str(i)
    
            for k in range(len(df_types)):
                parl[k].append(df_types[k][par].mean())
                parls[k].append(df_types[k][par].std())
    
    
    
        for k in range(len(df_types)):
            lab = dicKind
            plt.plot(times,parl[k],'-o',c=colors[k],label=lab[k])
            plt.errorbar(times[1:],np.array(parl[k])[1:],yerr=np.array(parls[k])[1:]/np.sqrt(num_samples),fmt='o',c=colors[k],alpha=0.3)
    
    
        
        plt.xlabel('Assembly event')
        if(name=='propM'):
            plt.legend(prop={'size': 10}) # We are using legend just for the first graphic, the rest is the same
            plt.ylabel('Proportion of mutualism')
        if(name=='propP'):
            plt.ylabel('Proportion of consumer-resource')
        if(name=='propC'):
            plt.ylabel('Proportion of competition')
        plt.show()
    ########################################################################
    
    
    # This is the same as above, but for the other variables
    ############################# OTHER VARIABLES ##################################
    #########################################################################
    df_types = []
    for k in df['kind'].unique():
        df_types.append(df[df['kind']==k].copy())
    
    for name in ['n_species','Con','r_dens','CS','r_P']: # Richness, connectance, average density, complexity, consumer-resource per species
        parl = []
        parls=[]
        times = tvec
        plt.figure(figsize=(4,4))
        for k in df['kind'].unique():
            parl.append([])
            parls.append([])
        for i in times:
            if(name=='r_dens'):
                par2= 'n_species'+str(i)
                par1= 'tot_density'+str(i)
    
                for k in range(len(df_types)):
                    parl[k].append((df_types[k][par1]/df_types[k][par2]).mean())
                    parls[k].append((df_types[k][par1]/df_types[k][par2]).std())
            
            elif(name=='CS'):
                par2= 'n_species'+str(i)
                par1= 'Con'+str(i)
    
                for k in range(len(df_types)):
                    parl[k].append((df_types[k][par1]*df_types[k][par2]).mean())
                    parls[k].append((df_types[k][par1]*df_types[k][par2]).std())
            
            elif(name=='r_P'):
                par2= 'n_species'+str(i)
                par1= 'P'+str(i)
    
                for k in range(len(df_types)):
                    parl[k].append((df_types[k][par1]/df_types[k][par2]).mean())
                    parls[k].append((df_types[k][par1]/df_types[k][par2]).std())
                    
                    
            else:    
                par = name+str(i)
            
    
                for k in range(len(df_types)):
                    parl[k].append(df_types[k][par].mean())
                    parls[k].append(df_types[k][par].std())
    
    
    
        for k in range(len(df_types)):
            lab = dicKind
            plt.plot(times,parl[k],'-o',c=colors[k],label=lab[k])
            plt.errorbar(times[1:],np.array(parl[k])[1:],yerr=np.array(parls[k])[1:]/np.sqrt(num_samples),fmt='o',c=colors[k],alpha=0.3)
    
    
        if(name=='n_species'):
            plt.ylabel('Species richness')
        if(name=='Con'):
            plt.ylabel('Connectance')
        if(name=='r_dens'):
            plt.ylabel('Average abundance')
        if(name=='CS'):
            plt.ylabel('Complexity (richness x connectance)')
        if(name=='r_P'):
            plt.ylabel('Consumer-resource per species')
        plt.xlabel('Assembly event')
        #plt.legend(prop={'size': 6})
        plt.show()
    
    plt.figure(figsize=(4,4))   
    for kind_idx in df['kind'].unique():
        for k in df[df['kind']==kind_idx].index: 
            if(k==df[df['kind']==kind_idx].index[-1]):
                plt.scatter(df.loc[k,'n_species'+str(tvec[-1])],df.loc[k,'Con'+str(tvec[-1])],color=colors[df.loc[k,'kind']],label=dicKind[df.loc[k,'kind']])
            else:
                plt.scatter(df.loc[k,'n_species'+str(tvec[-1])],df.loc[k,'Con'+str(tvec[-1])],color=colors[df.loc[k,'kind']])
    plt.xlabel('Species richness')
    plt.ylabel('Connectance')
    #plt.legend(prop={'size': 6})
    if(communityType==1):
        plt.legend(prop={'size': 10})
    plt.show()
    
    
    
    
    
    
    
    # This is to calculate values of entropy and modularity, then generate the graphics
    ######### DEGREE ENTROPY AND MODULARITY
    ######################################################################################
    
    
    
    # Create networks for each sample
    nets=[]
    for k in (df.index):
    
    
        RNm = GenerateNetworks(evo_history[k][0][-1])
    
        nets.append(RNm)
    
    
    # Define a sublist for each kind of sample
    SEntropy=[]
    Modul=[]
    for i in df['kind'].unique():
        SEntropy.append([])
        Modul.append([])
        
        
    
    for i in df['kind'].unique():
        for k in df.index:
            if(i==df.loc[k,'kind']):
                net = nets[k]
                
                degrees = np.array([val for (node, val) in net.degree()])
                Cr, n_nodes, H = NetworkMetrics(net,degrees)
                
                
                mod = Modularity(net)
    
                
                # This part calculates the values of the random network with same C and S, then their average, then subtracts it from the original value
                randomList = []
                for z in range(50):
                    randomList.append(nx.erdos_renyi_graph(n_nodes,Cr))
    
                r_sentropy=[]
                r_modul=[]
                for z in range(len(randomList)):
                    rdegrees = np.array([val for (node, val) in randomList[z].degree()])
                    rCr, rn_nodes, rH = NetworkMetrics(randomList[z],rdegrees)
                    r_sentropy.append(rH)
                    r_modul.append(Modularity(randomList[z]))
                
                SEntropy[i].append((H-np.array(r_sentropy).mean())/np.array(r_sentropy).mean())
                Modul[i].append((mod-np.array(r_modul).mean())/np.array(r_modul).mean())
                print(i,k)
            
    
    dicKind = {0:'Evo.\n($\Delta=5$)', 1:'Evo.\n($\Delta=30$)',2:'Inv.\nRand.',3:'Inv.\nProp.' }
    
    # Generate entropy graphic
    metric = 'Entropy increase'
    mdata=SEntropy.copy()
          
    rows_list=[]
    for kind in range(len(mdata)):
        for i in range(len(mdata[kind])):
            rows_list.append({metric:mdata[kind][i],'kind':dicKind[kind]})
    mdatadf = pd.DataFrame(rows_list)
    palette = sns.color_palette(colors)
    f = plt.figure(figsize=[4,4])
    ax = f.add_subplot(111)
    sns.boxplot(x='kind',y=metric,data=mdatadf,palette=palette,ax=ax, width=0.4,linewidth=1)
    plt.xlabel('Scenario')
    f.tight_layout()
    
    
    
    
    # Generate modularity graphic
    metric = 'Modularity increase'
    mdata=Modul.copy() 
    
    rows_list=[]
    for kind in range(len(mdata)):
        for i in range(len(mdata[kind])):
            rows_list.append({metric:mdata[kind][i],'kind':dicKind[kind]})
    mdatadf = pd.DataFrame(rows_list)
    palette = sns.color_palette(colors)
    f = plt.figure(figsize=[4,4])
    ax = f.add_subplot(111)
    sns.boxplot(x='kind',y=metric,data=mdatadf,palette=palette,ax=ax, width=0.4,linewidth=1)
    plt.xlabel('Scenario')
    f.tight_layout()