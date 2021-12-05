# Title     : TODO
# Objective : TODO
# Created by: Admin
# Created on: 2021/6/8

LMall <- function(NX,nt,nstep=30,order){

  stp <- (max(nt)-min(nt))/nstep
  res <- c()
  for(j in 1:nstep){

    tg1 <- Legendre.model11((j-1)*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg2 <- Legendre.model11(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg3 <- Legendre.model11(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg4 <- Legendre.model11(j*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tmp1 <- rbind(tg1,tg2,tg3,tg4)
    res <- rbind(res,tmp1)
  }
  res
}

fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){

  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL)
  sum((NG[,self]-(rowSums(odes)+NG[1,self]))^2)
}

ode.sovle.ind <- function(NG,fitpar,nconnect,nt,order,nstep,LL){

  stp <- (max(nt)-min(nt))/nstep
  index <- which(nconnect==1)

  ind.par <- matrix(fitpar[1:(length(index)*(order-1))],ncol=order-1,byrow=T)
  allrep  <- matrix(rep(0,length(index)),nrow=1)
  nn <- 1
  for(j in 1:nstep)
  {
    tmp_1 <- rowSums(t(apply(ind.par,1,"*",LL[nn,])))
    tmp_1[1] <- abs(tmp_1[1])
    tg1 <- tmp_1*NG[j,index]

    tmp_2 <- rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))
    tmp_2[1] <- abs(tmp_2[1])
    tg2 <- tmp_2*NG[j,index]

    tmp_3 <- rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))
    tmp_3[1] <- abs(tmp_3[1])
    tg3 <- tmp_3*NG[j,index]

    tmp_4 <- rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))
    tmp_4[1] <- abs(tmp_4[1])
    tg4 <- tmp_4*NG[j,index]

    #tmp_1 <- (t(apply(ind.par,1,"*",LL[nn,])))*NG[j,index]
    #tmp_1[1,] <- abs(tmp_1[1,])
    #tg1 <- (rowSums(tmp_1))
    #tg2 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))*NG[j,index])
    #tg3 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))*NG[j,index])
    #tg4 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))*NG[j,index])
    tmp <- allrep[j,] + stp*(tg1+2*tg2+2*tg3+tg4)/6
    allrep <- rbind(allrep,tmp)
    nn <- nn + 4
  }
  #allrep[,1] <-  abs(allrep[,1])
  return(allrep)
}

network_eff <- function(network_para,network_connect,effect,order,times,nstep){

  nt1 <- min(times)
  nt2 <- max(times)
  LL <- LMall(NX=1,nt=seq(nt1,nt2,(nt2-nt1)/nstep),nstep=nstep,order=order)

  nx <- dim(effect)[2]
  all_list <- c()
  for (y.c in 1:nx)
  {
    indexx    <- which(network_connect[y.c,]==1)
    effect_s  <- cbind(effect[,y.c], effect[,-y.c])
    nconnect  <- as.numeric(c(1, network_connect[y.c,-y.c]))
    fitpar    <- as.numeric(network_para[y.c,-1])
    B <- ode.sovle.ind(NG=(effect_s), fitpar = fitpar, nconnect = nconnect, nt = times, order = order, nstep = nstep, LL = LL)
    colnames(B) <- c(colnames(effect)[y.c], colnames(effect)[indexx])
    all_list    <- c(all_list, list(B))
  }
  return(all_list)
}

interType <- function(con,alle,sme){

  diag(con) <- 0
  nn <- dim(con)[1]
  connfest <- matrix(0,nrow=nn,ncol=nn)
  indp <- c()
  inter <- list()
  for(i in 1:nn){
    al <- alle[[i]]
    index <- which(as.numeric(colnames(al))==i)
    if(is.matrix(al)){
      lindp <- al[,index]
      linter <- al[,-index]
      indp <- cbind(indp,lindp)
      inter[[i]] <- linter
      rcor <- cor(as.numeric(sme[i,]),linter)
    }else{
      indp <- cbind(indp,al)
      inter[[i]] <- 0
      rcor <- 0
    }
    connfest[i,which(con[i,]==1)] <- as.numeric(rcor)
  }
  return(list(connfest=connfest,connect=con,indp=indp,inter=inter))
}

get_network_info <- function (connfest, nodes_info)
{
  after <- c()
  for (i in 1:dim(connfest)[1]){
    dep <- i
    ind <- which(connfest[i,]!=0)
    effect <- connfest[i,ind]
    one <- c()
    for (j in 1:length(ind)) {
      if(effect[j] >= 0){
        type <- '+'
      }else{
        type <- '-'
      }
      #one <- rbind(one,c(paste("M",ind[j]),paste("M", dep),abs(effect[j]),type))
      one <- rbind(one,c(nodes_info[ind[j]],nodes_info[dep],abs(effect[j]),type))
    }
    after <- rbind(after,one)
  }
  return(after)
}

Legendre.model11 <- function(t, np.order,tmin = NULL, tmax = NULL)
{
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- 1;
  if (np.order >= 2)
    L[2] <- 0.5 * (6 * ti)
  if (np.order >= 3)
    L[3] <- 0.5 * (15 * ti ^ 2 - 3)
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * 4 * ti ^ 3 - 60 * ti)
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
      ti)
  if (np.order >= 7)
    L[7] <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
      ti ^ 2 - 35)
  return(L);
}
