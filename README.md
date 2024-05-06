Example code of computing associations between gene values and clinical outcomes across CIDE cohorts.    

# Pre-requisites:  
1, Python >= 3.9 (you may install [Anaconda](https://www.anaconda.com/download))  
2, [R >= 4.3.0](https://cran.r-project.org)  
3, [logis_batch](https://github.com/data2intelligence/logis_batch)    

# Step 1: download open access data  
CD into the src folder and run "./download.py" to download files for open access data. For restricted data, users can directly unzip the data folders into "data/data_open" after getting them from the authors (with your approval notice for each individual dataset).   
  
# Step 2: Analysis for discrete outcomes  
Run "./run.py" to finish the analysis for all cohorts with binary or RECIST outcomes.  
  
# Step 3: Analysis for survival outcomes  
Run "./run.py inx 69" where inx is a number between 0 and 69 to compute the Cox-PH regressions. If you also include restricted datasets, please use the total number of cohorts (69 + number of restricted cohorts).  
Instead of using Python, we used R for Cox-PH regression as we found that R implementation of Cox-PH is more stable. This is a computational intensive step, because we use a sliding-value approach to find the best cutoff for Kaplan-Meier plots, thus involving lots of Cox-PH regressions. We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 69". You may re-write this file for your local HPC.

The output will be available in each cohort folder, named as [cohort name].[data type].response[endpoint type].    
