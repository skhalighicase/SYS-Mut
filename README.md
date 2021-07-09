# SYS-Mut
A model to Decode the Functional Significance of Rare Somatic Mutations in Cancer

##### For Academic/Non-Profit use
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>

##### For Technical Support:
Send us email at varadanlab@gmail.com, skhalighi@gmail.com

### Introduction
________________________________________________________________
<p> Large-scale tumor profiling studies and databases had helped us in discovering frequently mutated driver genes across multiple cancers also with the set of large and heterogeneous set of relatively rare somatic mutations whose bioligical significance remains largely unknown and less studied. SYS-Mut is a system biology tool that integrates multi-omic tumor data and uses our novel hierarchical Bayesian regression model to infer the functional impact of rare mutations on the transcriptional patterns of downstream targets while accounting for other confounding cis- and trans-regulatory factors.
</p>

### System Requirements
________________________________________________________________
### Platforms and System Requirements

Hard Disk Space : The module needs around 1GB of Hard Disk space.

The performance of the tool depends on the processing capabilities of the machine. 

The module performance was tested on the following platform

  * 16-core  2.10 GHz Intel(R) Xeon(R) CPU E5-2450 with 64 GB linux and Windows machines

### Software Architecture

The script was developed on the R platform with version 3.4, and it is assumed that the user also uses the same for running the script. 

### Software Requirements

 * [**R**](http://www.r-project.org/) : Free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms.
<br> preferred Version : >= 3.4

* [**JAGS**]JAGS: Just Another Gibbs Sampler (https://sourceforge.net/projects/mcmc-jags/) : It is a program for the statistical analysis of Bayesian hierarchical models by Markov Chain Monte Carlo.
<br> preferred Version : >= 3.4

________________________________________________________________
### Downloading SYS-Mut
SYS-Mut can be downloaded from github by recursively cloning the git repository
   
   $ git clone https://github.com/VaradanLab/SYS-Mut.git
   
### Preparing Data for the run

##### Input Files : All the files should be Independent matrices(Genes as rows and samples as columns)
 * Gene Expression File
 * Somatic Copy Number File
 * Methylation Data File
 * Somatic Mutation File
#### Users must note that for running the SYS-Mut All the four Input files should have same rows and columns in the same order.

 * SAMPLE_INFORMATION : This File consisting of information if sample is "Tumor" or "Normal"

________________________________________________________________
#### Network File
 SYS-Mut requires the the network information to be provided in a tab delimited file. This file consists of the Nodes Information followed  by the Interactions. 

________________________________________________________________
### Running SYS-Mut

You will need to run the "SYS-Mut.R" script included in the SYS-Mut package.

    $ EXPORT num_of_cores==8 #User have to provide the number of cores available for the process. 
    $ RScript SYS-Mut.R [-R] [location to SYS-Mut_PROJ_Config.txt]
for help :

    $ RScript SYS-Mut.R [-H]

To run SYS-Mut on HPC. Slurm Script should look the following. 

________________________________________________________________
    #SBATCH -J SYS-Mut_Job
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=8
    #SBATCH --output=RD_1ST.out
    #SBATCH --mem=20gb

    module load gcc
    export Number_of_Cores==cpus-per-task
    RScript SYS-Mut.R [-R] [location to SYS-Mut_PROJ_Config.txt]

_________________________________________________________________
### For any queries kindly mail us at varadanlab@gmail.com

