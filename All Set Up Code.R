
require(matlab)
require(expm)
require(rugarch)


############# VXO

setwd("C:/Users/mafagan/Documents");

vxo = read.csv(file="CBOE_VXO_8616.csv")

vxo.r_t = NULL

for (i in 1:(nrow(vxo)-1)){
  vxo.diff = log(vxo$Adj.Close[i]/vxo$Adj.Close[i+1])
  vxo.r_t = c(vxo.r_t,vxo.diff)
}

vxort = vxo.r_t

vxort[nrow(vxo)] = 0

vxo$Log.Diff = vxort

vxo.r_t = rev(vxo.r_t)

vxo.ts = ts(vxo.r_t,frequency = 252,start = c(1986,1))

vxo_ordered = vxo[rev(rownames(vxo)),]


######### SP100


sp100 = read.csv(file="SP100_Full_8316.csv")

sp100.r_t = NULL

for (i in 1:(nrow(sp100)-1)){
  sp100.diff = log(sp100$Adj.Close[i]/sp100$Adj.Close[i+1])
  sp100.r_t = c(sp100.r_t,sp100.diff)
}

sp100rt = sp100.r_t

sp100rt[nrow(sp100)] = 0

sp100$Log.Diff = sp100rt

sp100.r_t = rev(sp100.r_t)


sp100.ts = ts(sp100.r_t,frequency = 252,start = c(1983,1))

sp100_ordered = sp100[rev(rownames(sp100)),]


sp100_vxo = sp100_ordered[(3*252+3):nrow(sp100),]

sp100_vxo.ts = ts(sp100_vxo$Log.Diff,frequency = 252,start = c(1986,1))

par(mfrow=c(2,1))
plot(sp100_vxo.ts,type="l", main = "SP100 vs VXO", ylab = "SP100")#,ylim = c(-.5,1.5))
plot(vxo.ts,type="l", ylab = "VXO")#,ylim = c(-.5,1.5))

plot(sp100_vxo.ts,type="l", main = "SP100 vs VXO", ylab = "SP100",ylim = c(-.5,1.5))
plot(vxo.ts,type="l", ylab = "VXO",ylim = c(-.5,1.5))

par(mfrow=c(1,1))
plot(vxo.ts,type="l", main = "S&P 100 v VXO", col="black",lwd=1, ylab = "Return (%)")
lines(sp100_vxo.ts,type="l",col="cyan",lwd=1)
legend("topright",legend = c("VXO","SP100"),col = c("black","cyan"),lty = c(1,1),lwd=c(2,2))



########### VIX


VIX = read.csv(file="VIX_Hist_Ps.csv")

Vix.r_t = NULL

for (i in 1:(nrow(VIX)-1)){
  vix.diff = log(VIX$Adj.Close[i]/VIX$Adj.Close[i+1])
  Vix.r_t = c(Vix.r_t,vix.diff)
}

VIXrt = Vix.r_t

VIXrt[nrow(VIX)] = 0

VIX$Log.Diff = VIXrt

vix.r_t = rev(Vix.r_t)

vix.ts = ts(vix.r_t,frequency=252,start = c(1990,1))

###### SP 500 #######

SP = read.csv(file="SPFull.csv")

head(SP)

r_t = NULL

for (i in 1:(nrow(SP)-1)){
  log.diff = log(SP$Adj.Close[i]/SP$Adj.Close[i+1])
  r_t = c(r_t,log.diff)
}

r_t[nrow(SP)]=0

SP$Log.Diff = r_t


SP_90_16 = SP[1:(252*26),]

sp500.r_t = SP_90_16$Log.Diff

sp500.r_t = rev(sp500.r_t)

sp500.ts = ts(sp500.r_t,frequency = 252,start=c(1990,1))

par(mfrow=c(2,1))
plot(sp500.ts,type="l", main = "SP500 vs VIX", ylab = "SP500")#,ylim = c(-.5,1.5))
plot(vix.ts,type="l", ylab = "VIX")#,ylim = c(-.5,1.5))

plot(sp500.ts,type="l", main = "SP500 vs VIX", ylab = "SP500",ylim = c(-.5,.5))
plot(vix.ts,type="l", ylab = "VIX",ylim = c(-.5,.5))

par(mfrow=c(1,1))
plot(vix.ts,type="l", main = "S&P 500 v VIX", col="black",lwd=1, ylab = "Return (%)")
lines(sp500.ts,type="l",col="green",lwd=1)
legend("topright",legend = c("VIX","SP500"),col = c("black","green"),lty = c(1,1),lwd=c(2,2))

par(mfrow=c(4,1))
plot(sp500.ts, type="l",main="S&P 500",ylab="Return (%)")
plot(sp100_vxo.ts, type="l",main="S&P 100",ylab="Return (%)")
plot(vix.ts, type="l",main="VIX",ylab="Return (%)")
plot(vxo.ts, type="l",main="VXO",ylab="Return (%)")


#### Function Set Up


###################

T_mat_template = function(kbar){
  
  A = matrix(nrow = (2^kbar),ncol = (2^kbar),data=0);
  for( i in 0:(2^kbar-1)){
    for( j in i:(2^kbar-1)-i){ 
      
      A[i+1,j+1] = bitwXor(i,j);
      
    }
  } 
  
  A
}



transition_mat = function(A, input, kbar){
  
  b = 5.238742524 #1.661558594 #input[1]
  gamma_kbar = 0.934560319 #0.707649404 #input[3]
  kbar=9
  
  gamma = matrix(data = 0, nrow=kbar, ncol =2)
  gamma[1] = 1-(1-gamma_kbar)^(1/b^(kbar-1))
  
  for( i in 2:kbar ){
    gamma[i,1] = 1-(1-gamma[1,1])^(b^(i-1));
  }
  
  gamma = gamma*.5
  gamma[,2] = gamma[,1]
  gamma[,1] = 1 - gamma[,1]
  kbar1 = kbar+1
  kbar2 = 2^kbar
  prob = matrix(data=1,nrow = kbar2)
  
  for( i in 0:(2^kbar-1) ){  # Works out probability associated with each XOR number
    for( m in 1:kbar ) {
      prob[i+1,1] = prob[i+1,1] * gamma[kbar1-m, (as.numeric(intToBits(i)[m])+1)];
    }
  }  
  
  for( i in 0:(2^(kbar-1)-1) ){   # Copies probabilities to the transition matrix
    for( j in i:(2^(kbar-1)-1) ){
      A[kbar2-i,j+1] = prob[kbar2-A[i+1,j+1],1]; # Copies each probability to the other 8 symmetrical locations
      A[kbar2-j,i+1] =  A[kbar2-i,j+1];
      A[j+1,kbar2-i] =  A[kbar2-i,j+1];
      A[i+1,kbar2-j] =  A[kbar2-i,j+1];    
      A[i+1,j+1] = prob[A[i+1,j+1]+1,1];
      A[j+1,i+1] = A[i+1,j+1];
      A[kbar2-j,kbar2-i] = A[i+1,j+1];
      A[kbar2-i,kbar2-j] = A[i+1,j+1];
    }
  }
  
  A
  
}

gofm = function(input,kbar){
  
  kbar2 = 2^kbar
  
  m0 = input[2]
  m1=2-m0;
  g_m1 = 0:(kbar2-1);
  
  for( i in 1:(kbar2) ){
    g=1;
    for( j in 0:(kbar-1) ){
      if( bitwAnd(g_m1[i], (2^j)) != 0 ){
        g=g*m1
      }else{
        g=g*m0
      }
    }
    g_m1[i]=g;
  }
  
  g_m=sqrt(g_m1);
  
  g_m
  
}

pi_mat_extract = function( input, data, kbar){
  sigma = as.numeric(input[4])/sqrt(252);
  
  k2 =2^kbar;
  
  
  #############
  
  
  A_template = T_mat_template(kbar)
  
  
  require(matlab)
  
  A = transition_mat(A_template,input,kbar);
  g_m = gofm(input,kbar);
  Time = length(data);                      
  pi_mat = matrix(data=0,nrow=(Time+1),ncol = k2);     
  LLs = matrix(data=0,nrow=1,ncol=Time);
  pi_mat[1,] = (1/k2)*rep(1,times=k2);
  
  #*----------------------------------------------------------------------*
  #*                        Likelihood algorithm                          *
  #*----------------------------------------------------------------------*
  pa = (2*pi)^-0.5;
  s = repmat(sigma*g_m,Time,1);
  w_t = repmat(data,1,k2)#[1:(length(data)/2)],1,k2);
  w_t = pa*exp( - 0.5*((w_t/s)^2))/s; 
  w_t = w_t + 1e-16;
  
  for( t in 2:(Time+1)){          
    piA = (pi_mat[(t-1),]%*%A);
    C = (w_t[(t-1),]*piA); ft = rowSums(C);
    if( ft == 0 || is.nan(ft) || is.infinite(ft)){                     #This stop div by zero if probs are too low
      pi_mat[t,1] = 1;   
    }else{
      pi_mat[t,] = C / ft; 
    }
    
    LLs[,(t-1)] = log((w_t[(t-1),]%*%as.vector(piA)));
  } 
  
  LL=-rowSums(LLs);
  
  
  output = list(pi_mat=pi_mat,A=A,g_m=g_m)
  
  output
  
  
}





######################## Monte Carlo LL Analysis


getwd()
# 
# setwd("C:/Users/mafagan/Documents/Complexity Science")
# 
# LLs = read.table(file="Monte_Carlo_LL_Results", as.is = T)

setwd("C:/Users/mafagan/Documents/Complexity Science")
LLs = read.table(file="MC MSM LL Mat.csv", as.is = T, sep=",")
LLs = LLs[,-which(LLs<0)]

LLs_counter = matrix(data=NA,nrow=8,ncol=ncol(LLs))
LL_count_vect = NULL

for( column in 1:ncol(LLs) ){
  
   if( c(LLs[,column] == max(LLs[,column]))[5] == T ){
     LL_column = 5
   }else{
    LL_column = which.max(LLs[,column])
 }
  LL_count_vect = c(LL_count_vect,LL_column)
  
  LLs_counter[,column] = c(LLs[,column] == max(LLs[,column]))
}

times.k.chosen = c(sum(LL_count_vect==1)/length(LL_count_vect),
                  sum(LL_count_vect==2)/length(LL_count_vect),
                  sum(LL_count_vect==3)/length(LL_count_vect),
                  sum(LL_count_vect==4)/length(LL_count_vect),
                  sum(LL_count_vect==5)/length(LL_count_vect),
                  sum(LL_count_vect==6)/length(LL_count_vect),
                  sum(LL_count_vect==7)/length(LL_count_vect),
                  sum(LL_count_vect==8)/length(LL_count_vect))


Mean_LLs = matrix(NA, nrow=8,ncol=1)
Mean_diffs = matrix(NA, nrow=8,ncol=1)
SE_LLs = matrix(NA,nrow=8,ncol=1)

for( rows in 1:8 ){
  Mean_LLs[rows,] = mean(as.matrix(LLs[rows,]))
}
for( rows in 1:8 ){
  Mean_diffs[rows,] =Mean_LLs[rows,] -  Mean_LLs[5,]
}
for( rows in 1:8 ){
  SE_LLs[rows,] = sd(LLs[rows,])/sqrt(ncol(LLs))
}

Mean_LLs
Mean_diffs
SE_LLs


SE_sim = matrix(NA,nrow = 8,ncol=1)

for( thisrow in 1:8){
  SE_sim[thisrow,] = sd(LL_Diffs[thisrow,])/sqrt(ncol(LLs))
}

monte.carlo.LL.analysis.data = rbind(t(Mean_LLs),t(SE_LLs),times.k.chosen,t(Mean_diffs),t(SE_sim))
colnames(monte.carlo.LL.analysis.data) = c("kbar = 1","2","3","4","5","6","7","8")
rownames(monte.carlo.LL.analysis.data) = c("Mean LLs","LL SEs","% of times kbar = Max LL",
                                           "LL kbar=5 - LL kbar=x","SE of LL Differences")

monte.carlo.LL.analysis.data



############### Monte Carlo Parameter Analysis

setwd("C:/Users/mafagan/Documents/Complexity Science")

Params = read.table(file="MC MSM Param Mat.csv",as.is=T, sep = ",")

Params[4,] = (Params[4,])/sqrt(252)

TrueVals = matrix(nrow=4)
TrueVals[1] = 3
TrueVals[2] = 1.5
TrueVals[3] = .95
TrueVals[4] = .25

Params.mean = matrix(nrow=4)
Params.sse = matrix(nrow=4)
Params.rmse = matrix(nrow=4)
resid.sq = matrix(nrow=4,ncol=ncol(Params))

for( p in 1:4 ){
  Params.mean[p] = mean(as.matrix(Params[p,]))
  Params.sse[p] = sd(as.matrix(Params[p,]))/(sqrt(400))
  resid = Params[p,] - TrueVals[p]
  residsq = resid^2
  Params.rmse[p] = sqrt(mean(residsq))
}

monte.carlo.param.analysis.data = rbind(t(TrueVals),t(Params.mean),t(Params.sse),t(Params.rmse))
colnames(monte.carlo.param.analysis.data) = c("b","gamma_kbar","m0","sigma")
rownames(monte.carlo.param.analysis.data) = c("True Values", "Mean Simulated Value", "SE of simmed values","RMSe of simmed values")

monte.carlo.param.analysis.data



##################### GARCH v MSM LL Analysis

###### GARCH Fitting

spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE, archm = FALSE, 
                                    archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE))

stdspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE, archm = FALSE, 
                                       archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
                     distribution.model = "std")


sp5_garch_norm = ugarchfit(spec,sp500.r_t)
sp5_garch_std = ugarchfit(stdspec,sp500.r_t)

vix_garch_norm = ugarchfit(spec,vix.r_t)
vix_garch_std = ugarchfit(stdspec,vix.r_t)

sp1_garch_norm = ugarchfit(spec,sp100.r_t)
sp1_garch_std = ugarchfit(stdspec,sp100.r_t)

vxo_garch_norm = ugarchfit(spec,vxo.r_t)
vxo_garch_std = ugarchfit(stdspec,vxo.r_t)


#### Pulling out LLs for Model Comparison

sp5.garch.norm.LLs = -sp5_garch_norm@fit$log.likelihoods
sp1.garch.norm.LLs = -sp1_garch_norm@fit$log.likelihoods
vxo.garch.norm.LLs = -vxo_garch_norm@fit$log.likelihoods
vix.garch.norm.LLs = -vix_garch_norm@fit$log.likelihoods

sp5.garch.std.LLs = -sp5_garch_std@fit$log.likelihoods
sp1.garch.std.LLs = -sp1_garch_std@fit$log.likelihoods
vxo.garch.std.LLs = -vxo_garch_std@fit$log.likelihoods
vix.garch.std.LLs = -vix_garch_std@fit$log.likelihoods


getwd()
setwd("C:/Users/mafagan/Documents/Complexity Science")

MSM.LLs.sp5 = read.table('sp5_FINAL_MSM_LLs.csv',skip = 1,sep=",")
MSM.LLs.sp1 = read.table('sp1_FINAL_MSM_LLs.csv',skip = 1,sep=",")
MSM.LLs.vix = read.table('vix_FINAL_MSM_LLs.csv',skip = 1,sep=",")
MSM.LLs.vxo = read.table('vxo_FINAL_MSM_LLs.csv',skip = 1,sep=",")

sp5.MSM.LLs = NULL
sp1.MSM.LLs = NULL
vix.MSM.LLs = NULL
vxo.MSM.LLs = NULL

for(n in 1:length(MSM.LLs.sp5)){
  sp5.MSM.LLs = c(sp5.MSM.LLs,MSM.LLs.sp5[,n])
}
for(n in 1:length(MSM.LLs.sp1)){
  sp1.MSM.LLs = c(sp1.MSM.LLs,MSM.LLs.sp1[,n])
}
for(n in 1:length(MSM.LLs.vxo)){
  vxo.MSM.LLs = c(vxo.MSM.LLs,MSM.LLs.vxo[,n])
}
for(n in 1:length(MSM.LLs.vix)){
  vix.MSM.LLs = c(vix.MSM.LLs,MSM.LLs.vix[,n])
}

####### Control Panel for In Sample Comparison
# For MSM v Normal-GARCH, q  = 3; index.garch.LLs = index.garch.norm.LLs
# For MSM v Student-GARCH, q = 4; index.garch.LLs = index.garch.std.LLs

p = 4 # Number of Free Params for MSM. Always = 4 (b, m0, gamma_kbar, sigma) 
q = 3 # Number of Free Params for Garch. Norm = 3 (omega, alpha, beta); Student = 4 (omega, alpha, beta, v)

sp5.garch.LLs = sp5.garch.norm.LLs
sp1.garch.LLs = sp1.garch.norm.LLs
vxo.garch.LLs = vxo.garch.norm.LLs
vix.garch.LLs = vix.garch.norm.LLs

####### End Control Panel

### Vuong Tests

# Total Log-Likelihood comparison, adjusted for parameters
lrat.sp5 = sum(sp5.MSM.LLs) - sum(sp5.garch.LLs) - ((p/2)*log(length(sp5.MSM.LLs)) - (q/2)*log(length(sp5.garch.LLs)))
# Pointwise Log-Likelihood differences, calculated to obtain variance
li.sp5 = sp5.MSM.LLs - sp5.garch.LLs
omega2 = var(li.sp5)
# Adjust Log-Likelihood difference by square root of sample size
sp5.Vuong.score = sqrt(length(sp5.MSM.LLs))*(lrat.sp5/length(sp5.MSM.LLs))

# Calculate the Test's P-value with Z test
sp5.Vuong.Test = pnorm(sp5.Vuong.score,mean = 0,sd = sqrt(omega2))

# Want prob of observing value greater than sample value, so subtract by 1 (lower.tail = F would have worked as well)
sp5.Vuong.pval = 1 - sp5.Vuong.Test

## Rinse and Repeat for all indices

lrat.sp1 = sum(sp1.MSM.LLs) - sum(sp1.garch.LLs) - ((p/2)*log(length(sp1.MSM.LLs)) - (q/2)*log(length(sp1.garch.LLs)))
li.sp1 = sp1.MSM.LLs - sp1.garch.LLs
omega2 = var(li.sp1)
sp1.Vuong.score = sqrt(length(sp1.MSM.LLs))*(lrat.sp1/length(sp1.MSM.LLs))

sp1.Vuong.Test = pnorm(sp1.Vuong.score,mean = 0,sd = sqrt(omega2))

sp1.Vuong.pval = 1 - sp1.Vuong.Test


lrat.vix = sum(vix.MSM.LLs) - sum(vix.garch.LLs) - ((p/2)*log(length(vix.MSM.LLs)) - (q/2)*log(length(vix.garch.LLs)))
li.vix = vix.MSM.LLs - vix.garch.LLs
omega2 = var(li.vix)
vix.Vuong.score = sqrt(length(vix.MSM.LLs))*(lrat.vix/length(vix.MSM.LLs))

vix.Vuong.Test = pnorm(vix.Vuong.score,mean = 0,sd = sqrt(omega2))

vix.Vuong.pval = 1 - vix.Vuong.Test


lrat.vxo = sum(vxo.MSM.LLs) - sum(vxo.garch.LLs) - ((p/2)*log(length(vxo.MSM.LLs)) - (q/2)*log(length(vxo.garch.LLs)))
li.vxo = vxo.MSM.LLs - vxo.garch.LLs
omega2 = var(li.vxo)
vxo.Vuong.score = lrat.vxo/(sqrt(length(vxo.MSM.LLs)))

vxo.Vuong.Test = pnorm(vxo.Vuong.score,mean = 0,sd = sqrt(omega2))

vxo.Vuong.pval = 1 - vxo.Vuong.Test

####### Clarke Test

# Calculate difference in pointwise LLs, adjusted with a parameter penalty
di.sp5 = ( sp5.MSM.LLs - ((p/(2*length(sp5.MSM.LLs)))*
                            log(length(sp5.MSM.LLs)) )  ) - 
  ( (sp5.garch.LLs) - ( (q/(2*length(sp5.garch.LLs))) *
                          log(length(sp5.garch.LLs)) ) )

di.sp1 = ( sp1.MSM.LLs - ((p/(2*length(sp1.MSM.LLs)))*
                            log(length(sp1.MSM.LLs)) )  ) - 
  ( (sp1.garch.LLs) - ( (q/(2*length(sp1.garch.LLs))) *
                          log(length(sp1.garch.LLs)) ) )

di.vix = ( vix.MSM.LLs - ((p/(2*length(vix.MSM.LLs)))*
                            log(length(vix.MSM.LLs)) )  ) - 
  ( (vix.garch.LLs) - ( (q/(2*length(vix.garch.LLs))) *
                          log(length(vix.garch.LLs)) ) )

di.vxo = ( vxo.MSM.LLs - ((p/(2*length(vxo.MSM.LLs)))*
                            log(length(vxo.MSM.LLs)) )  ) - 
  ( (vxo.garch.LLs) - ( (q/(2*length(vxo.garch.LLs))) *
                          log(length(vxo.garch.LLs)) ) )

# Indicator to count number of positive differences
Bsp5 = length(di.sp5[di.sp5>=0])
# Total number of trials from which positive differences are compared
nsp5 = length(di.sp5)

# Binomial test to determine if probability of observed binomial trials is near .5
# Null Hypothesis is that test is <= .5; Alternative Hypothesis is that the test is > .5
sp5.nonp.test = binom.test(Bsp5,nsp5,p=.5,alternative = "greater")


Bsp1 = length(di.sp1[di.sp1>=0])
nsp1 = length(di.sp1)

sp1.nonp.test = binom.test(Bsp1,nsp1,p=.5,alternative = "greater")


Bvxo = length(di.vxo[di.vxo>=0])
nvxo = length(di.vxo)

vxo.nonp.test = binom.test(Bvxo,nvxo,p=.5,alternative = "greater")


Bvix = length(di.vix[di.vix>=0])
nvix = length(di.vix)

vix.nonp.test = binom.test(Bvix,nvix,p=.5,alternative = "greater")

MSM.LLs = c(sum(sp5.MSM.LLs),sum(sp1.MSM.LLs),sum(vix.MSM.LLs),sum(vxo.MSM.LLs))
GARCH.LLs = c(sum(sp5.garch.LLs),sum(sp1.garch.LLs),sum(vix.garch.LLs),sum(vxo.garch.LLs))

Vuong.pvals = rbind(sp5.Vuong.pval,sp1.Vuong.pval,vix.Vuong.pval,vxo.Vuong.pval)
Nonparam.pvals = rbind(sp5.nonp.test$p.value,sp1.nonp.test$p.value,vix.nonp.test$p.value,vxo.nonp.test$p.value)

All.pvals.norm = cbind(Vuong.pvals,Nonparam.pvals)

# Calculate BIC as defined by rugarch package. Takes standard BIC and scales it by the sample size of data points

sp5.BIC = (-2*(sum(sp5.MSM.LLs)) + p*log(length(sp5.MSM.LLs)))/length(sp5.MSM.LLs)
sp1.BIC = (-2*(sum(sp1.MSM.LLs)) + p*log(length(sp1.MSM.LLs)))/length(sp1.MSM.LLs)
vix.BIC = (-2*(sum(vix.MSM.LLs)) + p*log(length(vix.MSM.LLs)))/length(vix.MSM.LLs)
vxo.BIC = (-2*(sum(vxo.MSM.LLs)) + p*log(length(vxo.MSM.LLs)))/length(vxo.MSM.LLs)

MSM.BICs = c(sp5.BIC, sp1.BIC, vix.BIC, vxo.BIC)

sp5.g.BIC = (-2*(sum(sp5.garch.LLs)) + q*log(length(sp5.garch.LLs)))/length(sp5.garch.LLs)
sp1.g.BIC = (-2*(sum(sp1.garch.LLs)) + q*log(length(sp1.garch.LLs)))/length(sp1.garch.LLs)
vix.g.BIC = (-2*(sum(vix.garch.LLs)) + q*log(length(vix.garch.LLs)))/length(vix.garch.LLs)
vxo.g.BIC = (-2*(sum(vxo.garch.LLs)) + q*log(length(vxo.garch.LLs)))/length(vxo.garch.LLs)

GARCH.BICs = c(sp5.g.BIC, sp1.g.BIC, vix.g.BIC, vxo.g.BIC)

All.model.data = cbind(MSM.LLs,GARCH.LLs,MSM.BICs,GARCH.BICs,All.pvals.norm)

rownames(All.model.data) = c("S&P 500","S&P 100","VIX","VXO")

colnames(All.model.data) = c("MSM LLs","GARCH LLs", "MSM BICs", "GARCH BICs",
                             "Vuong Test P value","Clarke Test P value")

All.model.data

#### For the curious, this provides a kurtosis check to corroborate my choice of Clarke's test

require(moments)

kurtosis(sp5.MSM.LLs)
kurtosis(sp5.garch.std.LLs)
kurtosis(sp5.garch.norm.LLs)






################### Out of Sample Forecast Analysis

### First, refit GARCH models on half of each data set

sp500_garch_norm = ugarchfit(spec,sp500.r_t,out.sample = length(sp500.r_t)/2)
sp500_garch_std = ugarchfit(stdspec,sp500.r_t,out.sample = length(sp500.r_t)/2)

vix_garch_norm = ugarchfit(spec,vix.r_t, out.sample = length(vix.r_t)/2)
vix_garch_std = ugarchfit(stdspec,vix.r_t, out.sample = length(vix.r_t)/2)

sp100_garch_norm = ugarchfit(spec,sp100.r_t, out.sample = length(sp100.r_t)/2)
sp100_garch_std = ugarchfit(stdspec,sp100.r_t, out.sample = length(sp100.r_t)/2)

vxo_garch_norm = ugarchfit(spec,vxo.r_t, out.sample = length(vxo.r_t)/2)
vxo_garch_std = ugarchfit(stdspec,vxo.r_t, out.sample = length(vxo.r_t)/2)

#### MSM Fitting

## Refit MSM models by half of data set as well. Parameters calculated in MATLAB

setwd("C:/Users/mafagan/Documents/Complexity Science")

allpredparams = as.matrix(read.table('MSM_Prediction_Params.csv', header=TRUE,sep=","))

sp500.params = allpredparams[1,]
sp100.params = allpredparams[2,]
vix.params = allpredparams[3,]
vxo.params = allpredparams[4,]

# Extract necessary data from MSM estimation. Essentially, we estimate the parameters on half of the data, then use these parameters
# to evaluate data at time t. Then use available info at time t to forecast t+k.

# the function pi_mat_extract runs the MSM on the full data set to get an estimate of the composition of likely states for each t

sp500.msm.pred = pi_mat_extract(input = sp500.params,data = sp500.r_t,kbar=7)

sp100.msm.pred = pi_mat_extract(input = sp100.params,data = sp100.r_t,kbar=9)

vix.msm.pred = pi_mat_extract(input = vix.params,data = vix.r_t,kbar=10)

vxo.msm.pred = pi_mat_extract(input = vxo.params,data = vxo.r_t,kbar=9)







##### Setting Desired N-ahead Forecasts

n.ahead.vect = c(1,5,10,20,50)

############ Forecasting evaluation

for( n in n.ahead.vect ){
  
  #### SP500 MSM Forecast
  
  datanums = seq((length(sp500.r_t)/2+n),length(sp500.r_t),1) # Sets indexer to pull desired data
  
  sp500.pred.data2 = (sp500.r_t[datanums])^2
  
  pi_mat.sp500 = sp500.msm.pred$pi_mat[-1,] # first value of pi_mat is initializer
  half_pi.sp500 = pi_mat.sp500[(datanums-n),] # takes the latter half of pi_mat and scales it down by forecast horizon
                                              # example: with 1 day prediction, pi_mat = 3276, first data point = 3277
  
  require(expm) # Facilitates the creation of the forward transition mat by allowing matrix exponentiation
  
  A.sp500 = sp500.msm.pred$A%^%n  # take A from pi_mat_extract and raise it to prediction interval power
                                  # This is where markov property comes in. Apply discrete-time Kolmogorov forward equation
                                  # Xn+k = Xn * A^k
  gm.sp500 = sp500.msm.pred$g_m   # Pull M volatility multipliers in the g_m vector
  sigma.sp500 = sp500.params[4]/sqrt(252) # Convert volatility to daily level
  
  gm_preds = half_pi.sp500%*%A.sp500%*%(gm.sp500)^2 # Law of total expectation and application of Markov properties.
                                                    # Expected value of squared future returns = 
                                                    # distribution at time t %*% A ^ k forward periods %*% squared possible values
                                                    # E[squared multipler] @ time k = Xt %*% A^k %*% g_m^2, where %*% is the dotproduct
  
  msm_preds_sq = gm_preds*(sigma.sp500)^2 # Take E[squared multiplier] and multiply it by E[squared sigma] to get E[squared return]
  
  sp500.ols = lm(sp500.pred.data2~msm_preds_sq) # Regress predicted squared returns on to empricial squared returns
  
  sp500.coeffs = summary(sp500.ols)$coefficients # Take out coefficients for regression
  
  ### SP500 Garch Forecast
  
  # Use internal rugarch forecasting function
  # for guidelines , use ?ugarchforecast
  
  sp500.forecast.norm = ugarchforecast(sp500_garch_norm, n.ahead = n, n.roll=(length(sp500.r_t)/2)-n);
  sp500.forecast.std = ugarchforecast(sp500_garch_std, n.ahead = n, n.roll = (length(sp500.r_t)/2)-n);
  
  sp500.normgarch.preds = (sp500.forecast.norm@forecast$sigmaFor^2)[n,] # Extract squared returns at forecast horizon
  sp500.normgarch.ols = lm(sp500.pred.data2~sp500.normgarch.preds)
  
  sp500.stdgarch.preds = (sp500.forecast.std@forecast$sigmaFor^2)[n,]
  sp500.stdgarch.ols = lm(sp500.pred.data2~sp500.stdgarch.preds)
  
  sp500.normgarch.coefs = summary(sp500.normgarch.ols)$coefficients
  sp500.stdgarch.coefs = summary(sp500.stdgarch.ols)$coefficients
  
  ### Collect MSE and Adjusted R-Squared Values for All Models
  
  sp500.msm.r2 = summary(sp500.ols)$adj.r.squared
  sp500.msm.mse = mean((summary(sp500.ols)$residuals)^2)
  
  sp500.stdgarch.r2 = summary(sp500.stdgarch.ols)$adj.r.squared
  sp500.stdgarch.mse = mean((summary(sp500.stdgarch.ols)$residuals)^2)
  
  sp500.normgarch.r2 = summary(sp500.normgarch.ols)$adj.r.squared
  sp500.normgarch.mse = mean((summary(sp500.normgarch.ols)$residuals)^2)
  
  ### Test Significance of H0: Alpha == 0; Beta == 1
  
  sp500.msm.alpha.ptest = summary(sp500.ols)$coefficients[1,4]
  sp500.msm.beta.ptest = summary(lm(sp500.pred.data2 ~msm_preds_sq, offset= 1.00*msm_preds_sq))$coefficients[2,4]
  # By offsetting the beta coefficient by 1 I can use the lm p-value calculation on the null hypothesis that the coefficient is = 1
  
  sp500.normgarch.alpha.ptest = summary(sp500.normgarch.ols)$coefficients[1,4]
  sp500.normgarch.beta.ptest = summary(lm(sp500.pred.data2 ~sp500.normgarch.preds, offset= 1.00*sp500.normgarch.preds))$coefficients[2,4]
  
  sp500.stdgarch.alpha.ptest = summary(sp500.stdgarch.ols)$coefficients[1,4]
  sp500.stdgarch.beta.ptest = summary(lm(sp500.pred.data2 ~sp500.stdgarch.preds, offset= 1.00*sp500.stdgarch.preds))$coefficients[2,4]
  
  # Now, rinse and repeat for all indices
  
  ##### SP100 MSM Forecast
  
  datanums = seq((length(sp100.r_t)/2-.5+n),length(sp100.r_t),1)
  
  sp100.pred.data2 = (sp100.r_t[datanums])^2
  
  pi_mat.sp100 = sp100.msm.pred$pi_mat[-1,]
  half_pi.sp100 = pi_mat.sp100[(datanums-n),]
  
  require(expm)
  
  A.sp100 = sp100.msm.pred$A%^%n
  gm.sp100 = sp100.msm.pred$g_m
  sigma.sp100 = sp100.params[4]/sqrt(252)
  
  gm_preds = half_pi.sp100%*%A.sp100%*%(gm.sp100)^2
  
  msm_preds_sq = gm_preds*(sigma.sp100)^2
  head(msm_preds_sq)
  head(sp100.pred.data2)
  
  sp100.ols = lm(sp100.pred.data2~msm_preds_sq)
  
  sp100.coeffs = summary(sp100.ols)$coefficients
  
  
  ### SP100 Garch Forecast
  
  sp100.forecast.norm = ugarchforecast(sp100_garch_norm, n.ahead = n, n.roll=(length(sp100.r_t)/2)-n+.5);
  sp100.forecast.std = ugarchforecast(sp100_garch_std, n.ahead = n, n.roll = (length(sp100.r_t)/2)-n+.5);
  
  sp100.normgarch.preds = (sp100.forecast.norm@forecast$sigmaFor^2)[n,]
  sp100.normgarch.ols = lm(sp100.pred.data2~sp100.normgarch.preds)
  
  sp100.stdgarch.preds = (sp100.forecast.std@forecast$sigmaFor^2)[n,]
  sp100.stdgarch.ols = lm(sp100.pred.data2~sp100.stdgarch.preds)
  
  sp100.normgarch.coefs = summary(sp100.normgarch.ols)$coefficients
  sp100.stdgarch.coefs = summary(sp100.stdgarch.ols)$coefficients
  
  
  sp100.msm.r2 = summary(sp100.ols)$adj.r.squared
  sp100.msm.mse = mean((summary(sp100.ols)$residuals)^2)
  
  sp100.stdgarch.r2 = summary(sp100.stdgarch.ols)$adj.r.squared
  sp100.stdgarch.mse = mean((summary(sp100.stdgarch.ols)$residuals)^2)
  
  sp100.normgarch.r2 = summary(sp100.normgarch.ols)$adj.r.squared
  sp100.normgarch.mse = mean((summary(sp100.normgarch.ols)$residuals)^2)
  
  ######
  
  sp100.msm.alpha.ptest = summary(sp100.ols)$coefficients[1,4]
  sp100.msm.beta.ptest = summary(lm(sp100.pred.data2 ~msm_preds_sq, offset= 1.00*msm_preds_sq))$coefficients[2,4]
  
  sp100.normgarch.alpha.ptest = summary(sp100.normgarch.ols)$coefficients[1,4]
  sp100.normgarch.beta.ptest = summary(lm(sp100.pred.data2 ~sp100.normgarch.preds, offset= 1.00*sp100.normgarch.preds))$coefficients[2,4]
  
  sp100.stdgarch.alpha.ptest = summary(sp100.stdgarch.ols)$coefficients[1,4]
  sp100.stdgarch.beta.ptest = summary(lm(sp100.pred.data2 ~sp100.stdgarch.preds, offset= 1.00*sp100.stdgarch.preds))$coefficients[2,4]
  
  
  ##### VIX MSM Forecast
  
  
  datanums = seq((length(vix.r_t)/2-.5+n),length(vix.r_t),1)
  
  vix.pred.data2 = (vix.r_t[datanums])^2
  
  pi_mat.vix = vix.msm.pred$pi_mat[-1,]
  half_pi.vix = pi_mat.vix[(datanums-n),]
  
  require(expm)
  
  A.vix = vix.msm.pred$A%^%n
  gm.vix = vix.msm.pred$g_m
  sigma.vix = vix.params[4]/sqrt(252)
  
  gm_preds = half_pi.vix%*%A.vix%*%(gm.vix)^2
  
  msm_preds_sq = gm_preds*(sigma.vix)^2
  head(msm_preds_sq)
  head(vix.pred.data2)
  
  
  vix.ols = lm(vix.pred.data2~msm_preds_sq)
  
  vix.coeffs = summary(vix.ols)$coefficients
  
  
  ### VIX Garch Forecast
  
  vix.forecast.norm = ugarchforecast(vix_garch_norm, n.ahead = n, n.roll=(length(vix.r_t)/2)-n+.5);
  vix.forecast.std = ugarchforecast(vix_garch_std, n.ahead = n, n.roll = (length(vix.r_t)/2)-n+.5);
  
  vix.normgarch.preds = (vix.forecast.norm@forecast$sigmaFor^2)[n,]
  vix.normgarch.ols = lm(vix.pred.data2~vix.normgarch.preds)
  
  vix.stdgarch.preds = (vix.forecast.std@forecast$sigmaFor^2)[n,]
  vix.stdgarch.ols = lm(vix.pred.data2~vix.stdgarch.preds)
  
  vix.normgarch.coefs = summary(vix.normgarch.ols)$coefficients
  vix.stdgarch.coefs = summary(vix.stdgarch.ols)$coefficients
  
  #####
  
  vix.msm.r2 = summary(vix.ols)$adj.r.squared
  vix.msm.mse = mean((summary(vix.ols)$residuals)^2)
  
  vix.stdgarch.r2 = summary(vix.stdgarch.ols)$adj.r.squared
  vix.stdgarch.mse = mean((summary(vix.stdgarch.ols)$residuals)^2)
  
  vix.normgarch.r2 = summary(vix.normgarch.ols)$adj.r.squared
  vix.normgarch.mse = mean((summary(vix.normgarch.ols)$residuals)^2)
  
  #####
  
  vix.msm.alpha.ptest = summary(vix.ols)$coefficients[1,4]
  vix.msm.beta.ptest = summary(lm(vix.pred.data2 ~msm_preds_sq, offset= 1.00*msm_preds_sq))$coefficients[2,4]
  
  vix.normgarch.alpha.ptest = summary(vix.normgarch.ols)$coefficients[1,4]
  vix.normgarch.beta.ptest = summary(lm(vix.pred.data2 ~vix.normgarch.preds, offset= 1.00*vix.normgarch.preds))$coefficients[2,4]
  
  vix.stdgarch.alpha.ptest = summary(vix.stdgarch.ols)$coefficients[1,4]
  vix.stdgarch.beta.ptest = summary(lm(vix.pred.data2 ~vix.stdgarch.preds, offset= 1.00*vix.stdgarch.preds))$coefficients[2,4]
  
  
  ##### VXO MSM Forecast
  
  
  datanums = seq((length(vxo.r_t)/2-.5+n),length(vxo.r_t),1)
  
  vxo.pred.data2 = (vxo.r_t[datanums])^2
  
  pi_mat.vxo = vxo.msm.pred$pi_mat[-1,]
  half_pi.vxo = pi_mat.vxo[(datanums-n),]
  
  require(expm)
  
  A.vxo = vxo.msm.pred$A%^%n
  gm.vxo = vxo.msm.pred$g_m
  sigma.vxo = vxo.params[4]/sqrt(252)
  
  gm_preds = half_pi.vxo%*%A.vxo%*%(gm.vxo)^2
  
  msm_preds_sq = gm_preds*(sigma.vxo)^2
  head(msm_preds_sq)
  head(vxo.pred.data2)
  
  
  vxo.ols = lm(vxo.pred.data2~msm_preds_sq)
  
  vxo.coeffs = summary(vxo.ols)$coefficients
  
  
  ### VXO Garch Forecast
  
  vxo.forecast.norm = ugarchforecast(vxo_garch_norm, n.ahead = n, n.roll=(length(vxo.r_t)/2)-n+.5);
  vxo.forecast.std = ugarchforecast(vxo_garch_std, n.ahead = n, n.roll = (length(vxo.r_t)/2)-n+.5);
  
  vxo.normgarch.preds = (vxo.forecast.norm@forecast$sigmaFor^2)[n,]
  vxo.normgarch.ols = lm(vxo.pred.data2~vxo.normgarch.preds)
  
  vxo.stdgarch.preds = (vxo.forecast.std@forecast$sigmaFor^2)[n,]
  vxo.stdgarch.ols = lm(vxo.pred.data2~vxo.stdgarch.preds)
  
  vxo.normgarch.coefs = summary(vxo.normgarch.ols)$coefficients
  vxo.stdgarch.coefs = summary(vxo.stdgarch.ols)$coefficients
  
  ##### 
  
  vxo.msm.r2 = summary(vxo.ols)$adj.r.squared
  vxo.msm.mse = mean((summary(vxo.ols)$residuals)^2)
  
  vxo.stdgarch.r2 = summary(vxo.stdgarch.ols)$adj.r.squared
  vxo.stdgarch.mse = mean((summary(vxo.stdgarch.ols)$residuals)^2)
  
  vxo.normgarch.r2 = summary(vxo.normgarch.ols)$adj.r.squared
  vxo.normgarch.mse = mean((summary(vxo.normgarch.ols)$residuals)^2)
  
  ####
  
  vxo.msm.alpha.ptest = summary(vxo.ols)$coefficients[1,4]
  vxo.msm.beta.ptest = summary(lm(vxo.pred.data2 ~msm_preds_sq, offset= 1.00*msm_preds_sq))$coefficients[2,4]
  
  vxo.normgarch.alpha.ptest = summary(vxo.normgarch.ols)$coefficients[1,4]
  vxo.normgarch.beta.ptest = summary(lm(vxo.pred.data2 ~vxo.normgarch.preds, offset= 1.00*vxo.normgarch.preds))$coefficients[2,4]
  
  vxo.stdgarch.alpha.ptest = summary(vxo.stdgarch.ols)$coefficients[1,4]
  vxo.stdgarch.beta.ptest = summary(lm(vxo.pred.data2 ~vxo.stdgarch.preds, offset= 1.00*vxo.stdgarch.preds))$coefficients[2,4]
  
  #### Adjust Coefficient P-Values to reflect H0 of Beta = 1
  
  sp500.coeffs[2,4] = sp500.msm.beta.ptest
  sp500.normgarch.coefs[2,4] = sp500.normgarch.beta.ptest
  sp500.stdgarch.coefs[2,4] = sp500.stdgarch.beta.ptest
  
  sp100.coeffs[2,4] = sp100.msm.beta.ptest
  sp100.normgarch.coefs[2,4] = sp100.normgarch.beta.ptest
  sp100.stdgarch.coefs[2,4] = sp100.stdgarch.beta.ptest
  
  vix.coeffs[2,4] = vix.msm.beta.ptest
  vix.normgarch.coefs[2,4] = vix.normgarch.beta.ptest
  vix.stdgarch.coefs[2,4] = vix.stdgarch.beta.ptest
  
  vxo.coeffs[2,4] = vxo.msm.beta.ptest
  vxo.normgarch.coefs[2,4] = vxo.normgarch.beta.ptest
  vxo.stdgarch.coefs[2,4] = vxo.stdgarch.beta.ptest
  
  
  if( n == 1 ){
    
    ### Store all regression models and coefficients
    
    msm.one.step.ols = list(sp500.ols,sp100.ols,vix.ols,vxo.ols)
    msm.one.step.coefs = rbind(sp500.coeffs,sp100.coeffs,vix.coeffs,vxo.coeffs)
    
    normgarch.one.ols = list(sp500.normgarch.ols,sp100.normgarch.ols,vix.normgarch.ols,vxo.normgarch.ols)
    normgarch.one.coefs = rbind(sp500.normgarch.coefs,sp100.normgarch.coefs,vix.normgarch.coefs,vxo.normgarch.coefs)
    
    stdgarch.one.ols = list(sp500.stdgarch.ols,sp100.stdgarch.ols,vix.stdgarch.ols,vxo.stdgarch.ols)
    stdgarch.one.coefs = rbind(sp500.stdgarch.coefs,sp100.stdgarch.coefs,vix.stdgarch.coefs,vxo.stdgarch.coefs)
    
    ### Store all R-Squared and MSE Values
    
    msm.one.step.r2 = rbind(sp500.msm.r2,sp100.msm.r2,vix.msm.r2,vxo.msm.r2)
    msm.one.step.mse = rbind(sp500.msm.mse,sp100.msm.mse,vix.msm.mse,vxo.msm.mse)
    
    normgarch.one.step.r2 = rbind(sp500.normgarch.r2,sp100.normgarch.r2,vix.normgarch.r2,vxo.normgarch.r2)
    normgarch.one.step.mse = rbind(sp500.normgarch.mse,sp100.normgarch.mse,vix.normgarch.mse,vxo.normgarch.mse)
    
    stdgarch.one.step.r2 = rbind(sp500.stdgarch.r2,sp100.stdgarch.r2,vix.stdgarch.r2,vxo.stdgarch.r2)
    stdgarch.one.step.mse = rbind(sp500.stdgarch.mse,sp100.stdgarch.mse,vix.stdgarch.mse,vxo.stdgarch.mse)
    
    #### Store all Alpha and Beta Tests
    
    msm.one.step.a.tests = rbind(sp500.msm.alpha.ptest,sp100.msm.alpha.ptest,
                                 vix.msm.alpha.ptest,vxo.msm.alpha.ptest)
    msm.one.step.b.tests = rbind(sp500.msm.beta.ptest,sp100.msm.beta.ptest,vix.msm.beta.ptest,vxo.msm.beta.ptest)
    
    normgarch.one.step.a.tests = rbind(sp500.normgarch.alpha.ptest,sp100.normgarch.alpha.ptest,vix.normgarch.alpha.ptest,
                                       vxo.normgarch.alpha.ptest)
    normgarch.one.step.b.tests = rbind(sp500.normgarch.beta.ptest,sp100.normgarch.beta.ptest,vix.normgarch.beta.ptest,
                                       vxo.normgarch.beta.ptest)
    
    stdgarch.one.step.a.tests = rbind(sp500.stdgarch.alpha.ptest,sp100.stdgarch.alpha.ptest,vix.stdgarch.alpha.ptest,
                                      vxo.stdgarch.alpha.ptest)
    stdgarch.one.step.b.tests = rbind(sp500.stdgarch.beta.ptest,sp100.stdgarch.beta.ptest,vix.stdgarch.beta.ptest,
                                      vxo.stdgarch.beta.ptest)
    
  }else if( n == 5 ){
    msm.five.step.ols = list(sp500.ols,sp100.ols,vix.ols,vxo.ols)
    msm.five.step.coefs = rbind(sp500.coeffs,sp100.coeffs,vix.coeffs,vxo.coeffs)
    
    normgarch.five.ols = list(sp500.normgarch.ols,sp100.normgarch.ols,vix.normgarch.ols,vxo.normgarch.ols)
    normgarch.five.coefs = rbind(sp500.normgarch.coefs,sp100.normgarch.coefs,vix.normgarch.coefs,vxo.normgarch.coefs)
    
    stdgarch.five.ols = list(sp500.stdgarch.ols,sp100.stdgarch.ols,vix.stdgarch.ols,vxo.stdgarch.ols)
    stdgarch.five.coefs = rbind(sp500.stdgarch.coefs,sp100.stdgarch.coefs,vix.stdgarch.coefs,vxo.stdgarch.coefs)
    
    msm.five.step.r2 = rbind(sp500.msm.r2,sp100.msm.r2,vix.msm.r2,vxo.msm.r2)
    msm.five.step.mse = rbind(sp500.msm.mse,sp100.msm.mse,vix.msm.mse,vxo.msm.mse)
    
    normgarch.five.step.r2 = rbind(sp500.normgarch.r2,sp100.normgarch.r2,vix.normgarch.r2,vxo.normgarch.r2)
    normgarch.five.step.mse = rbind(sp500.normgarch.mse,sp100.normgarch.mse,vix.normgarch.mse,vxo.normgarch.mse)
    
    stdgarch.five.step.r2 = rbind(sp500.stdgarch.r2,sp100.stdgarch.r2,vix.stdgarch.r2,vxo.stdgarch.r2)
    stdgarch.five.step.mse = rbind(sp500.stdgarch.mse,sp100.stdgarch.mse,vix.stdgarch.mse,vxo.stdgarch.mse)
    
    #### Store all Alpha and Beta Tests
    
    msm.five.step.a.tests = rbind(sp500.msm.alpha.ptest,sp100.msm.alpha.ptest,
                                  vix.msm.alpha.ptest,vxo.msm.alpha.ptest)
    msm.five.step.b.tests = rbind(sp500.msm.beta.ptest,sp100.msm.beta.ptest,vix.msm.beta.ptest,vxo.msm.beta.ptest)
    
    normgarch.five.step.a.tests = rbind(sp500.normgarch.alpha.ptest,sp100.normgarch.alpha.ptest,vix.normgarch.alpha.ptest,
                                        vxo.normgarch.alpha.ptest)
    normgarch.five.step.b.tests = rbind(sp500.normgarch.beta.ptest,sp100.normgarch.beta.ptest,vix.normgarch.beta.ptest,
                                        vxo.normgarch.beta.ptest)
    
    stdgarch.five.step.a.tests = rbind(sp500.stdgarch.alpha.ptest,sp100.stdgarch.alpha.ptest,vix.stdgarch.alpha.ptest,
                                       vxo.stdgarch.alpha.ptest)
    stdgarch.five.step.b.tests = rbind(sp500.stdgarch.beta.ptest,sp100.stdgarch.beta.ptest,vix.stdgarch.beta.ptest,
                                       vxo.stdgarch.beta.ptest)
    
  }else if( n == 10 ){
    msm.ten.step.ols = list(sp500.ols,sp100.ols,vix.ols,vxo.ols)
    msm.ten.step.coefs = rbind(sp500.coeffs,sp100.coeffs,vix.coeffs,vxo.coeffs)
    
    normgarch.ten.ols = list(sp500.normgarch.ols,sp100.normgarch.ols,vix.normgarch.ols,vxo.normgarch.ols)
    normgarch.ten.coefs = rbind(sp500.normgarch.coefs,sp100.normgarch.coefs,vix.normgarch.coefs,vxo.normgarch.coefs)
    
    stdgarch.ten.ols = list(sp500.stdgarch.ols,sp100.stdgarch.ols,vix.stdgarch.ols,vxo.stdgarch.ols)
    stdgarch.ten.coefs = rbind(sp500.stdgarch.coefs,sp100.stdgarch.coefs,vix.stdgarch.coefs,vxo.stdgarch.coefs)
    
    msm.ten.step.r2 = rbind(sp500.msm.r2,sp100.msm.r2,vix.msm.r2,vxo.msm.r2)
    msm.ten.step.mse = rbind(sp500.msm.mse,sp100.msm.mse,vix.msm.mse,vxo.msm.mse)
    
    normgarch.ten.step.r2 = rbind(sp500.normgarch.r2,sp100.normgarch.r2,vix.normgarch.r2,vxo.normgarch.r2)
    normgarch.ten.step.mse = rbind(sp500.normgarch.mse,sp100.normgarch.mse,vix.normgarch.mse,vxo.normgarch.mse)
    
    stdgarch.ten.step.r2 = rbind(sp500.stdgarch.r2,sp100.stdgarch.r2,vix.stdgarch.r2,vxo.stdgarch.r2)
    stdgarch.ten.step.mse = rbind(sp500.stdgarch.mse,sp100.stdgarch.mse,vix.stdgarch.mse,vxo.stdgarch.mse)
    
    #### Store all Alpha and Beta Tests
    
    msm.ten.step.a.tests = rbind(sp500.msm.alpha.ptest,sp100.msm.alpha.ptest,
                                 vix.msm.alpha.ptest,vxo.msm.alpha.ptest)
    msm.ten.step.b.tests = rbind(sp500.msm.beta.ptest,sp100.msm.beta.ptest,vix.msm.beta.ptest,vxo.msm.beta.ptest)
    
    normgarch.ten.step.a.tests = rbind(sp500.normgarch.alpha.ptest,sp100.normgarch.alpha.ptest,vix.normgarch.alpha.ptest,
                                       vxo.normgarch.alpha.ptest)
    normgarch.ten.step.b.tests = rbind(sp500.normgarch.beta.ptest,sp100.normgarch.beta.ptest,vix.normgarch.beta.ptest,
                                       vxo.normgarch.beta.ptest)
    
    stdgarch.ten.step.a.tests = rbind(sp500.stdgarch.alpha.ptest,sp100.stdgarch.alpha.ptest,vix.stdgarch.alpha.ptest,
                                      vxo.stdgarch.alpha.ptest)
    stdgarch.ten.step.b.tests = rbind(sp500.stdgarch.beta.ptest,sp100.stdgarch.beta.ptest,vix.stdgarch.beta.ptest,
                                      vxo.stdgarch.beta.ptest)
    
  }else if( n == 20){
    msm.twenty.step.ols = list(sp500.ols,sp100.ols,vix.ols,vxo.ols)
    msm.twenty.step.coefs = rbind(sp500.coeffs,sp100.coeffs,vix.coeffs,vxo.coeffs)
    
    normgarch.twenty.ols = list(sp500.normgarch.ols,sp100.normgarch.ols,vix.normgarch.ols,vxo.normgarch.ols)
    normgarch.twenty.coefs = rbind(sp500.normgarch.coefs,sp100.normgarch.coefs,vix.normgarch.coefs,vxo.normgarch.coefs)
    
    stdgarch.twenty.ols = list(sp500.stdgarch.ols,sp100.stdgarch.ols,vix.stdgarch.ols,vxo.stdgarch.ols)
    stdgarch.twenty.coefs = rbind(sp500.stdgarch.coefs,sp100.stdgarch.coefs,vix.stdgarch.coefs,vxo.stdgarch.coefs)
    
    msm.twenty.step.r2 = rbind(sp500.msm.r2,sp100.msm.r2,vix.msm.r2,vxo.msm.r2)
    msm.twenty.step.mse = rbind(sp500.msm.mse,sp100.msm.mse,vix.msm.mse,vxo.msm.mse)
    
    normgarch.twenty.step.r2 = rbind(sp500.normgarch.r2,sp100.normgarch.r2,vix.normgarch.r2,vxo.normgarch.r2)
    normgarch.twenty.step.mse = rbind(sp500.normgarch.mse,sp100.normgarch.mse,vix.normgarch.mse,vxo.normgarch.mse)
    
    stdgarch.twenty.step.r2 = rbind(sp500.stdgarch.r2,sp100.stdgarch.r2,vix.stdgarch.r2,vxo.stdgarch.r2)
    stdgarch.twenty.step.mse = rbind(sp500.stdgarch.mse,sp100.stdgarch.mse,vix.stdgarch.mse,vxo.stdgarch.mse)
    
    #### Store all Alpha and Beta Tests
    
    msm.twenty.step.a.tests = rbind(sp500.msm.alpha.ptest,sp100.msm.alpha.ptest,
                                    vix.msm.alpha.ptest,vxo.msm.alpha.ptest)
    msm.twenty.step.b.tests = rbind(sp500.msm.beta.ptest,sp100.msm.beta.ptest,vix.msm.beta.ptest,vxo.msm.beta.ptest)
    
    normgarch.twenty.step.a.tests = rbind(sp500.normgarch.alpha.ptest,sp100.normgarch.alpha.ptest,vix.normgarch.alpha.ptest,
                                          vxo.normgarch.alpha.ptest)
    normgarch.twenty.step.b.tests = rbind(sp500.normgarch.beta.ptest,sp100.normgarch.beta.ptest,vix.normgarch.beta.ptest,
                                          vxo.normgarch.beta.ptest)
    
    stdgarch.twenty.step.a.tests = rbind(sp500.stdgarch.alpha.ptest,sp100.stdgarch.alpha.ptest,vix.stdgarch.alpha.ptest,
                                         vxo.stdgarch.alpha.ptest)
    stdgarch.twenty.step.b.tests = rbind(sp500.stdgarch.beta.ptest,sp100.stdgarch.beta.ptest,vix.stdgarch.beta.ptest,
                                         vxo.stdgarch.beta.ptest)
    
  }else{
    msm.fifty.step.ols = list(sp500.ols,sp100.ols,vix.ols,vxo.ols)
    msm.fifty.step.coefs = rbind(sp500.coeffs,sp100.coeffs,vix.coeffs,vxo.coeffs)
    
    normgarch.fifty.ols = list(sp500.normgarch.ols,sp100.normgarch.ols,vix.normgarch.ols,vxo.normgarch.ols)
    normgarch.fifty.coefs = rbind(sp500.normgarch.coefs,sp100.normgarch.coefs,vix.normgarch.coefs,vxo.normgarch.coefs)
    
    stdgarch.fifty.ols = list(sp500.stdgarch.ols,sp100.stdgarch.ols,vix.stdgarch.ols,vxo.stdgarch.ols)
    stdgarch.fifty.coefs = rbind(sp500.stdgarch.coefs,sp100.stdgarch.coefs,vix.stdgarch.coefs,vxo.stdgarch.coefs)
    
    msm.fifty.step.r2 = rbind(sp500.msm.r2,sp100.msm.r2,vix.msm.r2,vxo.msm.r2)
    msm.fifty.step.mse = rbind(sp500.msm.mse,sp100.msm.mse,vix.msm.mse,vxo.msm.mse)
    
    normgarch.fifty.step.r2 = rbind(sp500.normgarch.r2,sp100.normgarch.r2,vix.normgarch.r2,vxo.normgarch.r2)
    normgarch.fifty.step.mse = rbind(sp500.normgarch.mse,sp100.normgarch.mse,vix.normgarch.mse,vxo.normgarch.mse)
    
    stdgarch.fifty.step.r2 = rbind(sp500.stdgarch.r2,sp100.stdgarch.r2,vix.stdgarch.r2,vxo.stdgarch.r2)
    stdgarch.fifty.step.mse = rbind(sp500.stdgarch.mse,sp100.stdgarch.mse,vix.stdgarch.mse,vxo.stdgarch.mse)
    
    #### Store all Alpha and Beta Tests
    
    msm.fifty.step.a.tests = rbind(sp500.msm.alpha.ptest,sp100.msm.alpha.ptest,
                                   vix.msm.alpha.ptest,vxo.msm.alpha.ptest)
    msm.fifty.step.b.tests = rbind(sp500.msm.beta.ptest,sp100.msm.beta.ptest,vix.msm.beta.ptest,vxo.msm.beta.ptest)
    
    normgarch.fifty.step.a.tests = rbind(sp500.normgarch.alpha.ptest,sp100.normgarch.alpha.ptest,vix.normgarch.alpha.ptest,
                                         vxo.normgarch.alpha.ptest)
    normgarch.fifty.step.b.tests = rbind(sp500.normgarch.beta.ptest,sp100.normgarch.beta.ptest,vix.normgarch.beta.ptest,
                                         vxo.normgarch.beta.ptest)
    
    stdgarch.fifty.step.a.tests = rbind(sp500.stdgarch.alpha.ptest,sp100.stdgarch.alpha.ptest,vix.stdgarch.alpha.ptest,
                                        vxo.stdgarch.alpha.ptest)
    stdgarch.fifty.step.b.tests = rbind(sp500.stdgarch.beta.ptest,sp100.stdgarch.beta.ptest,vix.stdgarch.beta.ptest,
                                        vxo.stdgarch.beta.ptest)
    
  }
  
  
}

all.msm.coefs=cbind(msm.one.step.coefs,msm.five.step.coefs,msm.ten.step.coefs,msm.twenty.step.coefs,msm.fifty.step.coefs)
colnames(all.msm.coefs) = c("1-step estimate","1-step Std. Error","1-step T-value", "1-step P-value",
                            "5-step estimate","5-step Std. Error","5-step T-value", "5-step P-value",
                            "10-step estimate","10-step Std. Error","10-step T-value", "10-step P-value",
                            "20-step estimate","20-step Std. Error","20-step T-value", "20-step P-value",
                            "50-step estimate","50-step Std. Error","50-step T-value", "50-step P-value")
rownames(all.msm.coefs) = c("S&P 500 Alpha","S&P 500 Beta","S&P 100 Alpha","S&P 100 Beta",
                            "VIX Alpha","VIX Beta","VXO Alpha","VXO Beta")

all.normgarch.coefs = cbind(normgarch.one.coefs,normgarch.five.coefs,normgarch.ten.coefs,
                            normgarch.twenty.coefs,normgarch.fifty.coefs)
colnames(all.normgarch.coefs) = c("1-step estimate","1-step Std. Error","1-step T-value", "1-step P-value",
                                  "5-step estimate","5-step Std. Error","5-step T-value", "5-step P-value",
                                  "10-step estimate","10-step Std. Error","10-step T-value", "10-step P-value",
                                  "20-step estimate","20-step Std. Error","20-step T-value", "20-step P-value",
                                  "50-step estimate","50-step Std. Error","50-step T-value", "50-step P-value")
rownames(all.normgarch.coefs) = c("S&P 500 Alpha","S&P 500 Beta","S&P 100 Alpha","S&P 100 Beta",
                            "VIX Alpha","VIX Beta","VXO Alpha","VXO Beta")

all.stdgarch.coefs = cbind(stdgarch.one.coefs,stdgarch.five.coefs,stdgarch.ten.coefs,
                           stdgarch.twenty.coefs,stdgarch.fifty.coefs)
colnames(all.stdgarch.coefs) = c("1-step estimate","1-step Std. Error","1-step T-value", "1-step P-value",
                                 "5-step estimate","5-step Std. Error","5-step T-value", "5-step P-value",
                                 "10-step estimate","10-step Std. Error","10-step T-value", "10-step P-value",
                                 "20-step estimate","20-step Std. Error","20-step T-value", "20-step P-value",
                                 "50-step estimate","50-step Std. Error","50-step T-value", "50-step P-value")
rownames(all.stdgarch.coefs) = c("S&P 500 Alpha","S&P 500 Beta","S&P 100 Alpha","S&P 100 Beta",
                            "VIX Alpha","VIX Beta","VXO Alpha","VXO Beta")

all.msm.r2=cbind(msm.one.step.r2,msm.five.step.r2,msm.ten.step.r2,msm.twenty.step.r2,msm.fifty.step.r2)
colnames(all.msm.r2) = c("1-Step","5-Step","10-Step","20-Step","50-Step")
rownames(all.msm.r2) = c("SP500 MSM R^2","SP100 MSM R^2","VIX MSM R^2","VXO MSM R^2")

all.msm.mse=cbind(msm.one.step.mse,msm.five.step.mse,msm.ten.step.mse,msm.twenty.step.mse,msm.fifty.step.mse)
colnames(all.msm.mse) = c("1-Step","5-Step","10-Step","20-Step","50-Step")
rownames(all.msm.mse) = c("SP500 MSM MSE","SP100 MSM MSE","VIX MSM MSE","VXO MSM MSE")

all.normgarch.r2 = cbind(normgarch.one.step.r2,normgarch.five.step.r2,normgarch.ten.step.r2,
                         normgarch.twenty.step.r2,normgarch.fifty.step.r2)
colnames(normgarch.r2s) = c("1-Step","5-Step","10-Step","20-Step","50-Step")
rownames(normgarch.r2s) = c("SP500 Norm-GARCH R^2","SP100 Norm-GARCH R^2","VIX Norm-GARCH R^2","VXO Norm-GARCH R^2")

all.stdgarch.r2 = cbind(stdgarch.one.step.r2,stdgarch.five.step.r2,stdgarch.ten.step.r2,
                        stdgarch.twenty.step.r2,stdgarch.fifty.step.r2)
colnames(all.stdgarch.r2) = c("1-Step","5-Step","10-Step","20-Step","50-Step")
rownames(all.stdgarch.r2) = c("SP500 Std-GARCH R^2","SP100 Std-GARCH R^2","VIX Std-GARCH R^2","VXO Std-GARCH R^2")

all.normgarch.mse = cbind(normgarch.one.step.mse,normgarch.five.step.mse,normgarch.ten.step.mse,
                          normgarch.twenty.step.mse,normgarch.fifty.step.mse)
colnames(all.normgarch.mse) = c("1-Step","5-Step","10-Step","20-Step","50-Step")
rownames(all.normgarch.mse) = c("SP500 Norm-GARCH MSE","SP100 Norm-GARCH MSE","VIX Norm-GARCH MSE","VXO Norm-GARCH MSE")

all.stdgarch.mse = cbind(stdgarch.one.step.mse,stdgarch.five.step.mse,stdgarch.ten.step.mse,
                         stdgarch.twenty.step.mse,stdgarch.fifty.step.mse)
colnames(all.stdgarch.mse) = c("1-Step","5-Step","10-Step","20-Step","50-Step")
rownames(all.stdgarch.mse) = c("SP500 Std-GARCH MSE","SP100 Std-GARCH MSE","VIX Std-GARCH MSE","VXO Std-GARCH MSE")


all.prediction.stats = list(msm = all.msm.coefs, msm.r2 = all.msm.r2, msm.mse = all.msm.mse,
                             normgarch = all.normgarch.coefs, normgarch.r2 = all.normgarch.r2, normgarch.mse = all.normgarch.mse,
                             stdgarch = all.stdgarch.coefs, stdgarch.r2 = all.stdgarch.r2, stdgarch.mse = all.stdgarch.mse)
all.prediction.stats



###########


mse.tests.one = NULL
mse.tests.five = NULL
mse.tests.ten = NULL
mse.tests.twenty = NULL
mse.tests.fifty = NULL


for( i in 1:4 ){
  
  msm.resids = msm.one.step.ols[[i]]$residuals^2
  normgarch.resids = normgarch.one.ols[[i]]$residuals^2
  stdgarch.resids = stdgarch.one.ols[[i]]$residuals^2
  
  normgarch.stat = msm.one.step.mse[i] - normgarch.one.step.mse[i]
  norm.s=(var(normgarch.resids)+var(msm.resids))/length(msm.resids)
  
  stdgarch.stat = msm.one.step.mse[i] - stdgarch.one.step.mse[i]
  std.s = (var(stdgarch.resids)+var(msm.resids))/length(msm.resids)
  
  normgarch.score = normgarch.stat/sqrt(std.s)
  stdgarch.score = stdgarch.stat/sqrt(std.s)
  
  msm.v.normgarch = pnorm(normgarch.score,lower.tail = F)
  msm.v.stdgarch = pnorm(stdgarch.score, lower.tail = F)
  
  norm.std = c(msm.v.normgarch,msm.v.stdgarch)
  mse.tests.one = cbind(mse.tests.one,norm.std)
  
  
  msm.resids = msm.five.step.ols[[i]]$residuals^2
  normgarch.resids = normgarch.five.ols[[i]]$residuals^2
  stdgarch.resids = stdgarch.five.ols[[i]]$residuals^2
  
  normgarch.stat = msm.five.step.mse[i] - normgarch.five.step.mse[i]
  norm.s=(var(normgarch.resids)+var(msm.resids))/length(msm.resids)
  
  stdgarch.stat = msm.five.step.mse[i] - stdgarch.five.step.mse[i]
  std.s = (var(stdgarch.resids)+var(msm.resids))/length(msm.resids)
  
  normgarch.score = normgarch.stat/sqrt(std.s)
  stdgarch.score = stdgarch.stat/sqrt(std.s)
  
  msm.v.normgarch = pnorm(normgarch.score, lower.tail = F)
  msm.v.stdgarch = pnorm(stdgarch.score, lower.tail = F)
  
  norm.std = c(msm.v.normgarch,msm.v.stdgarch)
  mse.tests.five = cbind(mse.tests.five,norm.std)
  
  
  msm.resids = msm.ten.step.ols[[i]]$residuals^2
  normgarch.resids = normgarch.ten.ols[[i]]$residuals^2
  stdgarch.resids = stdgarch.ten.ols[[i]]$residuals^2
  
  normgarch.stat = msm.ten.step.mse[i] - normgarch.ten.step.mse[i]
  norm.s=(var(normgarch.resids)+var(msm.resids))/length(msm.resids)
  
  stdgarch.stat = msm.ten.step.mse[i] - stdgarch.ten.step.mse[i]
  std.s = (var(stdgarch.resids)+var(msm.resids))/length(msm.resids)
  
  normgarch.score = normgarch.stat/sqrt(std.s)
  stdgarch.score = stdgarch.stat/sqrt(std.s)
  
  msm.v.normgarch = pnorm(normgarch.score, lower.tail = F)
  msm.v.stdgarch = pnorm(stdgarch.score, lower.tail = F)
  
  norm.std = c(msm.v.normgarch,msm.v.stdgarch)
  mse.tests.ten = cbind(mse.tests.ten,norm.std)
  
  
  msm.resids = msm.twenty.step.ols[[i]]$residuals^2
  normgarch.resids = normgarch.twenty.ols[[i]]$residuals^2
  stdgarch.resids = stdgarch.twenty.ols[[i]]$residuals^2
  
  normgarch.stat = msm.twenty.step.mse[i] - normgarch.twenty.step.mse[i]
  norm.s=(var(normgarch.resids)+var(msm.resids))/length(msm.resids)
  
  stdgarch.stat = msm.twenty.step.mse[i] - stdgarch.twenty.step.mse[i]
  std.s = (var(stdgarch.resids)+var(msm.resids))/length(msm.resids)
  
  normgarch.score = normgarch.stat/sqrt(std.s)
  stdgarch.score = stdgarch.stat/sqrt(std.s)
  
  msm.v.normgarch = pnorm(normgarch.score, lower.tail = F)
  msm.v.stdgarch = pnorm(stdgarch.score, lower.tail = F)
  
  norm.std = c(msm.v.normgarch,msm.v.stdgarch)
  mse.tests.twenty = cbind(mse.tests.twenty,norm.std)
  
  
  msm.resids = msm.fifty.step.ols[[i]]$residuals^2
  normgarch.resids = normgarch.fifty.ols[[i]]$residuals^2
  stdgarch.resids = stdgarch.fifty.ols[[i]]$residuals^2
  
  normgarch.stat = msm.fifty.step.mse[i] - normgarch.fifty.step.mse[i]
  norm.s=(var(normgarch.resids)+var(msm.resids))/length(msm.resids)
  
  stdgarch.stat = msm.fifty.step.mse[i] - stdgarch.fifty.step.mse[i]
  std.s = (var(stdgarch.resids)+var(msm.resids))/length(msm.resids)
  
  normgarch.score = normgarch.stat/sqrt(std.s)
  stdgarch.score = stdgarch.stat/sqrt(std.s)
  
  msm.v.normgarch = pnorm(normgarch.score, lower.tail = F)
  msm.v.stdgarch = pnorm(stdgarch.score, lower.tail = F)
  
  norm.std = c(msm.v.normgarch,msm.v.stdgarch)
  mse.tests.fifty = cbind(mse.tests.fifty,norm.std)
  
}

all.mse.tests = rbind(mse.tests.one,mse.tests.five,mse.tests.ten,mse.tests.twenty,mse.tests.fifty)
colnames(all.mse.tests) = c("S&P 500","S&P 100", "VIX", "VXO")
rownames(all.mse.tests) = c("MSM v Norm-GARCH 1 step","MSM v Std-GARCH 1 step",
                       "MSM v Norm-GARCH 5 step","MSM v Std-GARCH 5 step",
                       "MSM v Norm-GARCH 10 step","MSM v Std-GARCH 10 step",
                       "MSM v Norm-GARCH 20 step","MSM v Std-GARCH 20 step",
                       "MSM v Norm-GARCH 50 step","MSM v Std-GARCH 50 step")
all.mse.tests