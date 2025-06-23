#!/usr/bin/env Rscript
library(tools)
source(file.path("~/workspace/BioAnalysis/R/", "survival_util.R"))

NSTEP = 100

commands = commandArgs(trailingOnly=T)

# three inputs here
data = commands[1]
survival = commands[2]
output = commands[3]

# check if data is a gzip file
if(file_ext(data) == 'gz'){data = gzfile(data)}

data = as.matrix(read.table(data, sep='\t', header=T, check.names=F, quote=NULL))
data = t(data) # R matrix is column major

survival = as.matrix(read.table(survival, sep='\t', header=T, check.names=F, quote=NULL))

margin = NULL
if(length(commands) > 3){
  margin = as.integer(commands[4])
}

# align matrix names
common = Reduce(intersect, list(rownames(data),rownames(survival)))
sprintf("%s samples", length(common))

data = data[common,,drop=F]
survival = survival[common,,drop=F]

# stop at low death rate
death_rate = sum(survival[,2])/dim(survival)[1]
if(death_rate < 0.05){
  warning("return by low death rate: ", death_rate)
  q()
  }

# split up survival and background
surv = Surv(survival[,1], survival[,2])

if(dim(survival)[2] > 2){
  B = survival[,3:dim(survival)[2], drop=F]
}else{
  B = survival[,c(), drop=F]
}

# build up regression data space
B = cbind(B, rep(0, dim(data)[1]))
B = as.data.frame(B)

N_B = ncol(B)
colnames(B)[N_B] = "pivot"

# iterate over features
features = colnames(data)
N = length(features)

result = NULL

step = round(max(N/NSTEP,1))

for (i in 1:N)
{
  # progress report
  if((i %% step) == 0){
    sprintf("%s", round(100 * i/N, 2))
  }
  
  fid = features[i]
  #if(!(tail(strsplit(fid, '@')[[1]], n=1) %in% c('AOAH', 'CR1L', 'COLQ', 'LY86', 'IFNG', 'ADAMTS7'))) next
  
  # part 1: overall regression
  arr = B[,N_B] = data[,i]
  
  arr_result = CoxPH_best_separation(B, surv, margin, NSTEP)
  
  if(sum(is.na(arr_result)) > 0){
    warning(paste0('Jump with failed continuous regression ', fid))
    next
  }
  
  if(is.null(result)){
    result = matrix(nrow = N, ncol=length(arr_result))
    colnames(result) = names(arr_result)
    rownames(result) = features
  }
  
  result[i,] = arr_result
}

if(is.null(result)){
  warning("stop by no meaning result")
  q()
}

mean.value = colMeans(data)
result = cbind(result, mean.value)

N = rep(length(common), dim(result)[1])
result = cbind(result, N)

writemat(result, output)

# temporary fix to add sample size
#result = readmat(output, F)
#N = rep(length(common), dim(result)[1])
#result = cbind(result, N)
#write.table(result, file=output, sep='\t', quote=F)
