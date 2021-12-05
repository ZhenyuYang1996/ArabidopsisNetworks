# Title     : TODO
# Objective : TODO
# Created by: Admin
# Created on: 2021/5/30
get_miu = function(miu_par)
{
  #u1 <- abs(miu_par[1])/(1+abs(miu_par[2])*exp(-abs(miu_par[3])*c(1:8)))
  #u2 <- abs(miu_par[4])/(1+abs(miu_par[5])*exp(-abs(miu_par[6])*c(1:8)))
  u <- miu_par[1]/(1+miu_par[2]*exp(-miu_par[3]*c(1:8)))
  return (u);
}
#=========================================================
get_sad_continuous = function(sad_par)
{
  n <- 8
  phi<- sad_par[1]
  v2 <- sad_par[2]
  tmp <- (1-phi^2)
  sigma <- array(1, dim=c(n,n))
  for(i in 1:n)
  {
    sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp
    sigma[i:n,i] <- sigma[i,i:n]
  }
  sigma <- sigma * abs(v2)
  return(sigma);
}
#=========================================================
H0_hypothesis = function(snp_data_i)
{
  H0_pheno = pheno_data[which(snp_data_i != 9), ]
  sample_m = colMeans(H0_pheno)
  LS_init_par = c(35.56889,55.140007,0.8)
  LS = optim(LS_init_par, fn = function (par, pheno){y <- get_miu(par);return(sum((pheno-y)^2))}, pheno = sample_m, method = "L-BFGS-B",lower = LS_init_par*0.01,upper = LS_init_par*50)
  H0_init_par = c(LS$par, c(0.6,0.6))
  H0_par = optim(par = H0_init_par, fn = function(par,pheno){return(-sum(dmvnorm(pheno, get_miu(par[1:3]), get_sad_continuous(par[4:5]), log = TRUE)))}, pheno = H0_pheno, method = "BFGS", control=list(maxit=200000))
  return(c(-H0_par$value, H0_par$par))
}


H1_hypothesis = function(snp_data_i,H1_init_par)
{
  H1_pheno_0 = pheno_data[which(snp_data_i == 0), ]
  H1_pheno_1 = pheno_data[which(snp_data_i == 1), ]
  H1_par = optim(par = H1_init_par, fn = function(par,pheno_0, pheno_1){SAD_M = get_sad_continuous(par[7:8]); return(    -sum(dmvnorm(pheno_0, get_miu(par[1:3]), SAD_M, log = TRUE)) - sum(dmvnorm(pheno_1, get_miu(par[4:6]), SAD_M, log = TRUE))   )},  pheno_0 = H1_pheno_0, pheno_1 = H1_pheno_1, method = "BFGS", control=list(maxit=200000))
  return(c(-H1_par$value, H1_par$par))
}

get_LR = function(i)
{
  snp_i = snp_data[i,]
  H0_estimated = H0_hypothesis(snp_data_i = snp_i);
  init_par = c(rep(H0_estimated[2:4],2),H0_estimated[5:6])
  H1_estimated = H1_hypothesis(snp_i, init_par)
  LR = 2*(H1_estimated[1] - H0_estimated[1])
  return(c(LR,H1_estimated))
}

