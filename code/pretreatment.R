# Title     : TODO
# Objective : TODO
# Created by: Admin
# Created on: 2021/5/7

cross_merge_data = function (data_1, data_2)
{
  nrows = dim(data_1)[1]
  ncols = dim(data_1)[2]
  HH_LL = matrix(0, nrow = nrows, ncol = 2*ncols)
  for (i in 1:nrows)
  {
    for (j in 1:ncols)
    {
      k = (j-1)*2 + 1
      HH_LL[i,k] = data_1[i,j]
      HH_LL[i,k+1] = data_2[i,j]
    }
  }
  return(HH_LL)
}
###根据求得的参数计算遗传标准差
get_VG <- function(FunMap_par,marker_data,t)
{
  diff_vg <- c()

  for (a in 1:dim(marker_data)[1]) {

    AA <- which(as.numeric(marker_data[a,])==1)
    aa <- which(as.numeric(marker_data[a,])==0)
    Aa <- which(as.numeric(marker_data[a,])==2)
    all <- which(as.numeric(marker_data[a,])!=9)

    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)

    p1 <- (NAA*2+NAa)/((NAA+NAa+Naa)*2) #A基因频率
    p0 <- (Naa*2+NAa)/((NAA+NAa+Naa)*2) #a基因频率

    mean_AA <- logistic_fun(as.numeric(FunMap_par[a,3:5]),t)
    mean_aa <- logistic_fun(as.numeric(FunMap_par[a,6:8]),t)
    AE <- (mean_AA - mean_aa)/2

    if(NAa==0){ Vg <- 2*p1*p0*(AE^2)  } else{

      mean_Aa <- logistic_fun(FunMap_par[a,][35:44],t)

      DE <- mean_Aa - (mean_AA + mean_aa)/2
      Vg <- 2*p1*p0*((AE + (p1 - p0)*DE)^2) + 4*p1*p1*p0*p0*DE*DE

    }
    diff_vg <- rbind(diff_vg,Vg)
    #cat(a,"finished","\n")
  }

  #colnames(diff_vg) <- seq(min(t),max(t),length.out = length(t))
  return(sqrt(diff_vg))
}

check_abnormal_par = function (all_par)
{
  n_par = dim(all_par)[2]
  illegal_list = list()
  for (i in 1:n_par)
  {
    mean_i  = mean(abs(all_par[,i]))
    sd_i   = sd(abs(all_par[,i]))
    illegal_list[i] = list(which(all_par[,i]>(mean_i+20*sd_i)))
  }
  return(illegal_list)
}

logistic_fun = function (par, times)
{
  a = abs(par[1])
  b = abs(par[2])
  r = abs(par[3])
  res <- a/(1+b*exp(-r*times))
  return(res)
}

out_sub_cluster_data = function (belong_data, ori_eff_data, out_file_name)
{
  all_cluster_num = length(table(belong_data[,2]))
  for (i in 1:all_cluster_num)
  {
    file_name = paste(out_file_name,"-",i,".txt",sep = "")
    cluster_i_index = which( belong_data[,2] == (i-1) )
    out_data = ori_eff_data[cluster_i_index, ]
    write.table(out_data, file = file_name, row.names = FALSE, col.names = FALSE)
  }
  return(TRUE)
}
out_ori_data_belong <- function (header, belong_data, ori_data, out_file_name)
{
  if(!file.exists(out_file_name))
  {
    out_data <- rep(0,nrow(ori_data))
    all_cluster_num = length(table(belong_data[,2]))
    for (i in 1:all_cluster_num)
    {
      feature = paste(all_cluster_num,"-",i,sep = "")
      cluster_i_index = which( belong_data[,2] == (i-1) )
      out_data[cluster_i_index] = feature
    }
    write.table(out_data, file = out_file_name, row.names = FALSE, col.names = FALSE)
  }else {
    out_data  = unlist(read.table(out_file_name))
    all_cluster_num = length(table(belong_data[,2]))
    for (i in 1:all_cluster_num)
    {
      feature = paste(all_cluster_num,"-",i,sep = "")
      cluster_i_index = which( belong_data[,2] == (i-1) )
      out_data[cluster_i_index] = feature
    }
    write.table(out_data, file = out_file_name, row.names = FALSE, col.names = FALSE)
  }
  return(TRUE)
}
out_lasso_data = function (belong_data, ori_eff_data, out_file_name)
{
  all_cluster_num = length(table(belong_data[,2]))
  for (i in 1:all_cluster_num)
  {
    file_name = paste(out_file_name,"_",i,".txt",sep = "")
    cluster_i_index = which( belong_data[,2] == (i-1) )
    out_data = ori_eff_data[cluster_i_index, ]
    write.table(out_data, file = file_name, row.names = FALSE, col.names = FALSE)
  }
  return(TRUE)
}

smooth.optim <- function(times,para,y,nt=seq(1,6,length=150))
{
  allpar <- c()
  smooth.d <- c()
  dsmooth.d <- c()
  for(i in 1:dim(y)[1])
  {
    L <- optim(para,smL,DS1=y[i,],times=times,method="BFGS")
    allpar <- rbind(allpar,c(L$par,L$value))
    smooth.d <- rbind(smooth.d,Legendre.model(nt,L$par))
    dsmooth.d <- rbind(dsmooth.d,dLegendre.model(nt,L$par))
  }

  list(allpar=allpar,smooth.d=smooth.d,dsmooth.d=dsmooth.d)
}

eff_smooth <- function(i)
{
  para = rep(.1,6)
  L <- optim(par = para,smL,DS1=data[i,],times=times,method="BFGS")
  smooth.d <- Legendre.model(nt,L$par)
  return(smooth.d)
}

smL <- function(times,para,DS1){

  sum((DS1-Legendre.model(t=times,mu=para))^2)
}

dLegendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1]*0 + 1*mu[2];
  if (np.order>=2)
    L <- L + 0.5 * (6 * ti)* mu[3] ;
  if (np.order>=3)
    L <- L +0.5 * (15 * ti ^ 2 - 3)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)* mu[5];
  if (np.order>=5)
    L <- L + 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)*mu[6];
  if (np.order>=6)
    L <- L + (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *ti)* mu[7];
  if (np.order>=7)
    L <- L + (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)* mu[8];
  return(L);
}

Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}

