


dat_F <- function(marker,pheno,pi1=1:8,pi2=9:16){

  ind <- pheno[,1]
  n_ind <- length(ind)
  npheno <- pheno[,-1]
  ind_t <- colnames(marker)[-(1:3)]
  nmarker <- marker[,-c(1:3)]
  ind_t1 <- c()
  for(i in 1:length(ind_t)){
    ind_t1 <- c(ind_t1,as.numeric(strsplit(ind_t[i],"X")[[1]][2]))
  }
  nm_i <- c()
  for(i in 1:n_ind){
    nm_i <- c(nm_i,which(ind[i]==ind_t1))
  }
  nmarker1 <- nmarker[,nm_i]
  pheno1 <- npheno[,pi1];pheno2 <- npheno[,pi2]

  return(list(Marker=nmarker1,H=pheno1,D=pheno2))
}

varf <- function(p1,p2,marker){

  nm <- dim(marker)[1]
  eff_h <- c()
  eff_d <- c()
  for(i in 1:nm){
    sm <- marker[i,]
    if(any(sm=="--")){
      sm[which(sm=="--")] <- NA
    }
    ni <- names(table(as.character(as.matrix(sm))))
    dat1 <- c();dat2 <- c();
    for(j in 1:length(ni)){
      dat1 <- rbind(dat1,colMeans(p1[which(sm==ni[j]),]))
      dat2 <- rbind(dat2,colMeans(p2[which(sm==ni[j]),]))
    }
    h <- apply(dat1,2,var)
    d <- apply(dat2,2,var)
    eff_h <- cbind(eff_h,h)
    eff_d <- cbind(eff_d,d)
  }

 return(list(eh=eff_h,ed=eff_d))
}


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
  allrep <- matrix(rep(0,length(index)),nrow=1)
  nn <- 1
  for(j in 1:nstep){
    tg1 <- (rowSums(t(apply(ind.par,1,"*",LL[nn,])))*NG[j,index])
    tg2 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))*NG[j,index])
    tg3 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))*NG[j,index])
    tg4 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))*NG[j,index])
    tmp <- allrep[j,] +stp*(tg1+2*tg2+2*tg3+tg4)/6
    allrep <- rbind(allrep,tmp)
    nn <- nn + 4
  }
  allrep
}

optim.parallel <- function(connect,effect,n.cores,proc,order,times,nstep){

  diag(connect) <- 1
  nt1 <- min(times)
  nt2 <- max(times)

  LL <- LMall(NX=1,nt=seq(nt1,nt2,(nt2-nt1)/nstep),nstep=nstep,order=order)

  nx <- dim(effect)[2]

  grp <- floor(nx/n.cores)
  grp.i <- c()
  if(n.cores==1){
    grp.i <- c(grp.i,rep(1,nx))
  }else{
    for(ii in 1:n.cores){
      if(ii==n.cores){
        grp.i <- c(grp.i,rep(ii,nx-grp*(ii-1)))
      }else{
        grp.i <- c(grp.i,rep(ii,grp))
      }
    }
  }

  grp.ii <- unique(grp.i)

  res.list <- mclapply(grp.ii, function(i)
  {
    y.c <- 	which(grp.i==i)
    A <- sapply(y.c, proc, connect=connect,effect=effect,LL=LL,nstep=nstep,order=order,times=times);
    return (list(A));
  }, mc.cores=n.cores )

  res1 <- do.call("c", res.list)
  #res2 <- parallel.data.optim(res1,connect,times)
  return(res1)
}

parallel.data.optim <- function(rd,nm,ntt){

  nrd <- matrix(rd,nrow=length(ntt))
  nn <- dim(nm)[1]
  ki <- 0
  allist <- list()
  for(i in 1:nn){
    iii <- (which(nm[i,]==1))
    iiil <- length(iii)
    tmp.d <- nrd[,(ki+1):(ki+iiil)]
    if(is.matrix(tmp.d)){
      colnames(tmp.d) <- iii
    }else{
      names(tmp.d) <- iii
    }

    allist[[i]] <- tmp.d
    ki <- ki + iiil
  }

  return(allist)
}


ode.optim <- function(y.c,connect,effect,LL,nstep,order,times){
  res_2 <- rep(0,(order-1) * dim(connect)[1])
  indexx <- which(connect[y.c,]==1)
  para <- rep(0.0001,length(indexx)*(order-1))
  res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
               LL=LL,method="BFGS",control=list(maxit=10,trace=T))
  cat("Gene=",y.c," ",res$value,"\n")
  #A <- ode.sovle.ind(NG=(effect),res$par,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,LL=LL)
  res_2[1:length(res$par)] <- res$par
  return(c(res$value, res_2))
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

