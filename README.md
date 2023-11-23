Example code of run CIDE analysis and reproduce main results  

# Step 1: download open access data  
Please CD into the src folder and run "./download.py" to download files for open access data. For restricted data, users can directly unzip the data folders into "data/data_open" after getting them from the authors (with your approval notice from data managers of each individual study).   
  
  
# Step 2: You may run "./run.py inx 31" where inx is a number between 0 and 31 to compute the Tres results for each dataset.    
We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 31". You may re-write this file for your local HPC.

The output will be available in data/output, named with the cancer type followed with database accession ID.  
