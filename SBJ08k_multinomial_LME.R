cwd = '/home/knight/PRJ_Error'
setwd(cwd)

library(mlogit)
#library(hdf5r)
library(R.matlab)
library(stringr)
library(ordinal)

data_file <- 'data/GRP/stats/EpnRPE_DifFB_mLME_st0t6_WL05_WS25_HGm_F25t121_zbtS_sm0_l1_wn50_hfa_chancoef.mat'
beta.chan <- readMat(data_file)$beta.chan

chancat.ix <- beta.chan[[13]]
chancat.labels <- beta.chan[[14]]
chan.labels <- beta.chan[[9]]

nroi1 <- length(chan.labels[[1]][[1]]) 
nroi2 <- length(chan.labels[[2]][[1]]) 
d <- data.frame(
                category = rep('nonsignificant', nroi1 + nroi2),
                region = c(rep(c('MPFC'), nroi1),rep(c('INS'),nroi2))
                )
d$subject <- c(str_split_fixed(chan.labels[[1]][[1]], ' ',2)[,1],
               str_split_fixed(chan.labels[[2]][[1]], ' ',2)[,1])

d$subject <- str_split_fixed(d$subject,'"',2)[,2]

for (r in 1:length(chancat.ix)){
  for (chl in 1:length(chancat.labels)){
    d$category[chancat.ix[[r]][[1]][[chl]][[1]] + nroi1*(r-1)] <- chancat.labels[[chl]][[1]]
  }
}
d$category <- factor(d$category, levels = c('uRPE','sRPE','pRPE','nRPE','nonsignificant'))
d$region <- as.factor(d$region)
d$subject <- as.factor(d$subject)
## Run the model
#d2 <- dfidx(d, shape = 'wide', varying = 2:3, choice='category')
#m1 <- mlogit(category ~ region | subject, data = d2, reflevel = 'uRPE')
m1 <- clmm2(category ~ 1, nominal = ~ region, random = subject, data = d, Hess = TRUE)
m0 <- clmm2(category ~ 1, random = subject, data = d, Hess = TRUE)


