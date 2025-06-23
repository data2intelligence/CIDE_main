suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(Hmisc))

penalty = 1e-3

readmat = function(mat, flag)
{
  mat = as.matrix(read.table(mat, sep='\t', header=T, check.names=F, quote=NULL))
  
  if(flag){
    symbol = mat[,"Symbol"]
    mat = mat[,1:(ncol(mat)-1)]
    class(mat) = "numeric"
    return (list("symbol" = symbol, "mat" = mat))
  }else{
    return (mat)
  }
}

writemat = function(result, output)
{
  colnames(result)[colnames(result) == "Pr(>|z|)"] = "p"
  result = result[!is.na(result[,"p"]),,drop=F]
  FDR = p.adjust(result[,"p"], method="fdr")
  result = cbind(result, FDR)
  write.table(result, file=output, sep='\t', quote=F)
}


CoxPH_best_separation = function(X, Y, margin, NSTEP)
{
  # part 1: continuous regression
  errflag = F
  
  # still keep warning results below
  X = as.matrix(X)
  
  coxph.fit = tryCatch(
    #coxph(Y~., data=X),
    coxph(Y~ridge(X, theta=penalty)),
    
    error = function(e){
      errflag <<- T
      }
    )
  
  #warning = function(w){
  #  errflag <<- T
  #}
  
  if(errflag) return (NA)
  
  n_r = nrow(X)
  n_c = ncol(X)
  
  arr_result = summary(coxph.fit)$coef
  rownames(arr_result) = colnames(X)
  
  z = arr_result[, 'coef']/arr_result[, 'se(coef)']
  arr_result = cbind(arr_result, z)
  
  arr_result = arr_result[n_c,]
  if(is.na(arr_result["z"])) return (NA)
  
  # no need to find optimal threshold
  if(is.null(margin)) return (arr_result)
  
  # part 2: find the optimal threshold
  arr = X[, n_c]
  
  vthres_arr = sort(arr)[(margin+1):(n_r-margin)]
  
  # for really long array, we don't need that resolution, shrink down to fixed length
  if(length(vthres_arr) > NSTEP){
    vthres_arr = seq(vthres_arr[1], tail(vthres_arr, n=1), length.out=NSTEP)
  }
  
  # these are missing values, not NULL not existing values
  zscore_opt = thres_opt = NA
  
  for(vthres in vthres_arr)
  {
    X[, n_c] = as.numeric(arr >= vthres)
    
    errflag = F
    coxph.fit = tryCatch(
      #coxph(Y~., data=X),
      coxph(Y~ridge(X, theta=penalty)),
      error = function(e) errflag <<- T
      )
    # warning = function(w) errflag <<- T
    if(errflag) next
    
    #z = summary(coxph.fit)$coef[n_c, "z"]
    
    r_temp = summary(coxph.fit)$coef[n_c,]
    z = r_temp["coef"] / r_temp["se(coef)"]
    
    if(is.na(z)) next
    
    if (is.na(zscore_opt)){
      zscore_opt = z
      thres_opt= vthres
    
    }else if(arr_result['z'] > 0){
      if(z > zscore_opt){
        zscore_opt = z
        thres_opt = vthres
      }
      
    }else{ # arr_result['z'] <= 0
      if(z < zscore_opt){
        zscore_opt = z
        thres_opt = vthres
      }
    }
  }
  
  arr_result['thres.opt'] = thres_opt
  arr_result['z.opt'] = zscore_opt
  
  return (arr_result)
}

partial_correlation = function(X, Y, B)
{
   if(dim(B)[2] > 0){
      X = lm(X~B)$residuals
      Y = lm(Y~B)$residuals
   }
   
   result = apply(X, 2, cor.test, y=Y)
   
   result = cbind(
      unlist(lapply(result, function(x) x$estimate)),
      unlist(lapply(result, function(x) x$statistic)),
      unlist(lapply(result, function(x) x$p.value))
   )
   
   colnames(result) = c('r', 't', 'p')
   rownames(result) = colnames(X)
   return (result)
}
