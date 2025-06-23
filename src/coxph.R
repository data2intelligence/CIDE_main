#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(Hmisc))

commands = commandArgs(trailingOnly = T)

mat = read.table(commands[1], sep='\t', header=T, check.names=F)
N = dim(mat)[2]

output = commands[2]

# stop at low death rate
death_rate = sum(mat[,2])/dim(mat)[1]
if(death_rate < 0.1){
  warning("return by low death rate")
  q()
  }

# assume the first and second columns are survival days and censor status
surv = Surv(mat[,1], mat[,2])

factors = as.data.frame(mat[,3:N, drop=F])

coxph.fit = coxph(surv~., data=factors)
reg.summary = summary(coxph.fit)$coef
rownames(reg.summary) = gsub('`','',rownames(reg.summary))

coxph.zph = cox.zph(coxph.fit, global=F, transform='rank')$table
rownames(coxph.zph) = gsub('`','',rownames(coxph.zph))

common = Reduce(intersect, list(rownames(reg.summary),rownames(coxph.zph)))
reg.summary = reg.summary[common,,drop=F]
coxph.zph = coxph.zph[common,,drop=F]

result = cbind(reg.summary, coxph.zph)
write.table(result, output, sep='\t', quote=F)


if(length(commands) >= 3)
{
  display_ratio = 3
  
  pdf(paste(output, "pdf", sep='.'), family="Helvetica")
  par(mar=c(5,5,5,1))
  
  inx = as.numeric(commands[3])
  arr = sort(factors[,inx])
  
  len = length(arr)
  
  value_r = c(arr[3*len/4], arr[len/4])
  
  if(sum(arr==1) + sum(arr==0) == len){value_r = c(1,0)}
  
  adjust_value = colMeans(factors)
  adjust_value = rbind(adjust_value, adjust_value)
  adjust_value[,inx] = value_r
  
  ltypes = c(1,2)
  lcols = c("red", "blue")
  lnames = rownames(adjust_value) = c("High", "Low")
  
  z = reg.summary[inx,"z"]
  pvalue = reg.summary[inx,"Pr(>|z|)"]
  z = format.df(z, 3, numeric.dollar=F)
  pvalue = format.df(pvalue, 3, numeric.dollar=F)
  title = paste0("Z=", z, " p=", pvalue)
  
  coxph.survfit = survfit(coxph.fit, newdata=as.data.frame(adjust_value), conf.type="none")
  plot(coxph.survfit, lty=ltypes, col=lcols, lwd=display_ratio, cex.lab=display_ratio, cex.axis=display_ratio, cex.main=display_ratio, xlab=colnames(mat)[1], ylab="Fraction", main=title)
  
  legend("topright", lnames, lty=ltypes, col=lcols, text.col = lcols, cex=display_ratio, lwd=display_ratio)
}
