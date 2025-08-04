This repository contains the code used in the paper:

**The eco-evolutionary assembly of complex communities with multiple interaction types**, by Gui Araujo and Miguel Lurgi
================

Be sure to adjust parameters as explained in comments and to check the paths for saving and loading files

**#1**. Model: fevomodel.py
Run this code to sample and save simulations.

**#2**. Dataframe: EcovoDataframe.py
Run this code to load simulations and generate a dataframe with the assembly histories of variables

**#3**. Plots: EcoEvoAnalysis.py
Run this code to load dataframes and produce the plots in Figures 1, 3, S1, S2

**#4**. Data: EcoEvoMicData.py
Run this code to load dataframes and timeseries data to simulate and generate plots of macroecological patterns in Figure 4

**#5**. Ternary: EcoEvoTernary.py
Run this code to produce the ternary plot in Figure 2

**spec-file.txt** -> use this file to recreate the conda environment with packages and versions used in this project

The idea is to generate a series of samples for each scenario using #1, separating the files by folders. Then using #2, read all samples of the same scenario into a single dataframe, for all types. Then using #3/#4/#5, read dataframes to generate the plots. Read headings and comments in the scripts for more in-depth explanations.

The main 2 options to choose are the variables isTest and isHeatmap. Choose isTest=1 to test custom parameters and run your own simulations. Use isHeatmap=1 to generate plots in Figure 1. 

If you choose isTest=0, then you will run pre-made scenarios to reproduce the results of the paper.

All the folders are ready to store the generated data, inside the folders 'data' and 'df'. For example, to run the 8 scenarios in the paper, just choose isTest=0 and isHeatmap=0 in all scripts and run first #1 then #2. In this case, after running #1, the files inside the data folder will have the simulation data. After running #2, the folder df will have the dataframes. Then, running #3/#4/#5 will produce the plots.

To change any variable, follow the comments explaining what they are. Be sure to change the affected variables in other scripts. For example, to change the number of assembly events, change the variable n_evolutions in fevomodel. Then, you will have to change accordingly the tvec variable in the other scripts.

This material provides the original df's used in the paper. Instead of running #1 and #2, it is possible to generate plots with the same df's from the main text. By running #1 and #2 with new simulations, new df's will overwrite the old ones (these df's are larger files not included in the github repository, only in zenodo: **https://doi.org/10.5281/zenodo.15621027**).

When downloading from github, first include the contents from data_folders.zip into the same folder containing the scripts.

This project was supported by the Leverhulme Trust through Research Project Grant **# RPG-2022-114**.



**CREATING THE ENVIRONMENT TO RUN THE CODE ON WINDOWS**
------------

Installation should take less than 10 minutes

1. Install the latest version of miniconda (for minimal setup): https://docs.conda.io/en/latest/miniconda.html

2. Open anaconda command prompt and go to the extracted folder: cd your-local-path/ecoevo_codes

3. Run the specs file to create your environment (accept all questions asked): conda create --EcoEvoEnv --file spec-file.txt

4. Activate the created environment: conda activate EcoEvoEnv

5. Run codes. As a simple demo, run: python EcoEvoAnalysis.py. This will quickly create plots for Type 2 scenarios, as shown in figures 1, 3, S1, and S2 (S figures take about 10 minutes to generate)


6. Running the entire simulations in fevomodel.py can take several hours per simulation. To quickly test the simulation code, use isTest=1 in all files you want to run, then run fevomodel.py with n_evolutions=10 for both num_sample=1 and num_sample=2. Then, change all tvec variables to [1,2,3,4,5,6,7,8,9] and run EcoEvoDataframe.py. Now you cam run other files to generate any desired plots. These simulations are simple and quick test running for only 10 assembly events.
