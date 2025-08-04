## ---------------------------
##
## Script name: EcoEvoAnalysis.py
##
## Purpose of script: This code loads the dataframes of simulations and creates the ternary plot in the paper
## It combines both type 1 and type 2 simulations together into the same plot
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





isTest=0 # Choose one to run custom simulations




df = pd.DataFrame()

if(isTest==1): # This block reads custom df's, choose folders, names, and parameters accordingly
    folder = 'test'
    filecodes = [[0,0]]
    num_samples=2
    communityType=1

    for i in filecodes:
        path=dir_path+'/df/'+folder+'_df_1/'+str(i[0])+'-df-'+str(i[1])+'.pkl'
        with open(path, 'rb') as file:
              
            sollist=pickle.load(file)
        
            df = pd.concat([df,sollist[0]], ignore_index = True)


    df.loc[df['kind']==0,'kind']='Test 0'
    df = df.loc[:,['kind','propC4','propP4','propM4','n_species4','Con4']]
    df.rename(columns={'propC4': 'C', 'propP4': 'P', 'propM4': 'M'}, inplace=True)
    dic = {'Test 0': 'blue'}
    
    
     
if(isTest==0): # This block reads df's for the 8 scenarios

    # This lists the types of simulations to load (numbers on the left; leave the ones on the right as 0)
    # 0: evo5, 1: evo30, 2: invR, 3: invP
    filecodes = [[0,0],[1,0],[2,0],[3,0]]
    num_samples=20
    communityType = 1 # This is to load type 1 scenarios first

    folder = 't'+str(communityType)

    for i in filecodes:
        path=dir_path+'/df/'+folder+'_df_1/'+str(i[0])+'-df-'+str(i[1])+'.pkl'
        with open(path, 'rb') as file:
              
            sollist=pickle.load(file)
        
            df = pd.concat([df,sollist[0]], ignore_index = True)
            
    
    df.loc[df['kind']==0,'kind']='T1 Evo. 5'
    df.loc[df['kind']==1,'kind']='T1 Evo. 30'
    df.loc[df['kind']==2,'kind']='T1 Inv. Rand.'
    df.loc[df['kind']==3,'kind']='T1 Inv. Prop.'
            
    communityType = 2 # Then it loads type 2 scenarios

    folder = 't'+str(communityType)

    for i in filecodes:
        path=dir_path+'/df/'+folder+'_df_1/'+str(i[0])+'-df-'+str(i[1])+'.pkl'
        with open(path, 'rb') as file:
              
            sollist=pickle.load(file)
        
            df = pd.concat([df,sollist[0]], ignore_index = True)
            
    df.loc[df['kind']==0,'kind']='T2 Evo. 5'
    df.loc[df['kind']==1,'kind']='T2 Evo. 30'
    df.loc[df['kind']==2,'kind']='T2 Inv. Rand.'
    df.loc[df['kind']==3,'kind']='T2 Inv. Prop.'


    
    df = df.loc[:,['kind','propC499','propP499','propM499','n_species499','Con499']] # This selects only relevant variables in the last time step




    # This chooses more appropriate names for columns and chooses all names and colours
    df.rename(columns={'propC499': 'C', 'propP499': 'P', 'propM499': 'M'}, inplace=True)
    dic = {'T1 Evo. 5': 'cyan', 'T1 Evo. 30': 'lime', 'T1 Inv. Rand.': 'burlywood', 'T1 Inv. Prop.': 'magenta', 'T2 Evo. 5': 'blue', 'T2 Evo. 30': 'green', 'T2 Inv. Rand.': 'gold', 'T2 Inv. Prop.': 'red'}
    
   
    

# This block creates the ternary plot used in the paper
import ternary
import matplotlib.pyplot as plt

# Set up figure and axes
scale = 1  # determines the size of the triangle
figure, tax = ternary.figure(scale=scale)
figure.set_size_inches(5, 4)



# Custom styling
tax.boundary(linewidth=2)
tax.gridlines(color="black", multiple=.10, linewidth=0.5)
tax.ticks(axis='lbr', linewidth=1, multiple=.20, offset=0.02, fontsize=12, tick_formats="%.1f")
tax.clear_matplotlib_ticks()
#tax.legend()

# Plot points grouped by 'kind'
for kind, group in df.groupby('kind'):
    points = group[['M', 'C', 'P']].values.tolist()
    color = dic.get(kind, 'gray')
    tax.scatter(points, label=kind, color=color, s=60, alpha=0.8)

# Labels
tax.left_axis_label("prop. C-R", fontsize=12, offset=0.15)
tax.right_axis_label("prop. C", fontsize=12, offset=0.15)
tax.bottom_axis_label("prop. M", fontsize=12, offset=0.1)
#tax.set_background_color("lightgrey")


#tax.legend(fontsize=9, loc='upper right', handletextpad=0.4, frameon=False)
#tax.get_axes().get_legend().remove()

tax.get_axes().axis('off')  
plt.tight_layout()
plt.savefig('ternary.png', dpi=300, bbox_inches='tight')
plt.show()