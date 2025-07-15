This repository contains the code used in the paper: **The eco-evolutionary assembly of complex communities with multiple interaction types**

Be sure to adjust parameters as explained in comments and to check the paths for saving and loading files

**#1. Model:** fevomodel.p
Run this code to sample nd save simulations.

**#2. Dataframe:** EcovoDataframe.py
Run this code to load simulations and generate a dataframe with the assembly histories of variables

**#3. Plots:** EcoEvoAnalysis.py
Run this code to load dataframes and produce the plots in Figures 1, 3, S1, S2

**#4. Data:** EcoEvoMicData.py
Run this code to load dataframes and timeeries data to simulate and generate plots of macroecological patterns in Figure 4

**#5. Ternary:** EcoEvoTernary.py
Run this code to produce the ternary plot in Figure 2

**spec-file.txt** -> use this file to recreate the conda environment with packages and versions used in this project

The idea is to generate a series of samples for each scenario using #1, separating the files by folders. Then using #2, read all samples of the same scenario into a single dataframe, for all types. Then using #3/#4/#5, read dataframes to generate the plots. Read headings and comments in the scripts for more in-depth explanations.

The main 2 options to choose are the variables isTest and isHeatmap. Choose isTest=1 to test custom parameters and run your own simulations. Use isHeatmap=1 to generate plots in Figure 1. 

If you choose isTest=0, then you will run pre-made scenarios to reproduce the results of the paper.

All the folders are ready to store the generated data, inside the folders 'data' and 'df'. For example, to run the 8 scenarios in the paper, just choose isTest=0 and isHeatmap=0 in all scripts and run first #1 then #2. In this case, after running #1, the files inside the data folder will have the simulation data. After running #2, the folder df will have the dataframes. Then, running #3/#4/#5 will produce the plots.

To change any variable, follow the comments explaining what they are. Be sure to change the affected variables in other scripts. For example, to change the number of assembly events, change the variable n_evolutions in fevomodel. Then, you will have to change accordingly the tvec variable in the other scripts.

This material provides the original df's used in the paper. Istead of running #1 and #2, it is possible to generate plots with the same df's from the main text. By running #1 and #2 with new simulations, new df's will overwrite the old ones (these df's are larger files not included in the github repository, only in zenodo:  https://doi.org/10.5281/zenodo.15621028). 

When downloading from github, first include the contents from **data_folders.zip** into the same folder containing the scripts.

This project was supported by the Leverhulme Trust through Research Project Grant **\# RPG-2022-114**.
