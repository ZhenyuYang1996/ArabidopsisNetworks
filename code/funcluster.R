library("wavelets")
library("mvtnorm")

set.seed(1)

plot.wv <- function(Theta, Data, type=c("all.means", "single.cluster"), j)
{
  r <- log(ncol(Data$Y)/ncol(Data$Z), 2)
  W <- Theta$W
  J <- ncol(W)

  iw <- function(i, r)
  {
    wnew <- W[,i]
    j <- 1
    f <- function(i, ws) c(ws[i]/sqrt(2), ws[i]/sqrt(2))
    while(j <= r) {
      wnew <- as.vector(sapply(1:length(wnew), f, wnew))
      j <- j+1
    }
    return(wnew)
  }

  means.plot <- function(Theta, Data)
  {
    MU <- sapply(1:J, iw, r=r)
    clrs <- rainbow(J)

    plot(Data$tm, MU[,1], type="n", ylim=c(min(MU), max(MU)), 			xlab = "Time (min)",
         ylab="Normalized Expression")
    for(i in 1:J) points(Data$tm, MU[,i], type="l", lwd=3,
                         col=clrs[i])
  }

  cluster.plot <- function(Theta, Data, j)
  {
    MU <- iw(j,r)

    ind <- which(Theta$P[,j] > .9)
    cldat <- Data$Y[ind,]

    plot(Data$tm, cldat[1,], ylim=c(min(cldat), max(cldat)), type="n",
         xlab = "Time (min)",
         ylab = "Normalized Expression",
         main = paste("Cluster", j))
    for(i in 1L:nrow(cldat)) {
      points(Data$tm, cldat[i,], pch = 19, col="#0000FF30")
    }
    lines(Data$tm, MU, col="black")
  }

  type <- match.arg(type)
  switch(type,
         all.means=means.plot(Theta, Data),
         single.cluster=cluster.plot(Theta, Data, j)
  )
}


update.p <- function(Theta, Data, MP)
{
  ngene <- MP$ngene
  POld <- Theta$P
  sigsq <- Theta$Var$sigsq
  W <- Theta$W
  Z <- Data$Z
  J <- MP$J
  omega <- Theta$omega

  Sig <- build.sig(Theta, Data, MP)

  oneterm <- function(i,j)
  {
    f <- function(i,j) {
      omega[j] * dmvnorm(Z[i,], W[,j], Sig)
    }
    mapply(f, i, j)
  }

  PRaw <- outer(1:ngene, 1:J, oneterm)
  denom <- rowSums(PRaw)
  PNew <- PRaw/denom
  PNew[is.nan(PNew)] <- POld[is.nan(PNew)]
  PNew
}

update.omega <- function(Theta) colMeans(Theta$P)

filter.y <- function(Data, MP){
  r <- MP$r
  ngene <- MP$ngene
  Y <- Data$Y
  matrix(dwt(t(Y), "haar", r , fast = TRUE)@V[[paste("V", r, sep="")]], nrow=ngene, byrow=TRUE)
}

update.W <- function(Theta, Data, MP) {
  Z <- Data$Z
  P <- Theta$P
  (t(Z) %*% P) %*% diag(1/colSums(P))
}
#Result should be a K by J matrix, where K each colum is the
#new estimated mean for group j

build.sig <- function(Theta, Data, MP, rho) {
  if(missing(rho)) rho <- Theta$Var$rho
  M <- MP$M
  sigsq <- Theta$Var$sigsq

  sigsq*toeplitz(rho^(0:(M-1)))
}

build.siginv <- function(Theta, Data, MP, rho) {
  M <- MP$M
  sigsq <- Theta$Var$sigsq
  if(missing(rho)) rho <- Theta$Var$rho

  tmp <- toeplitz(c(1+rho^2, -rho, rep(0, M-2)))
  tmp[1,1] <- tmp[M,M] <- 1
  tmp/(sigsq*(1-rho^2))
}

update.sigsq <- function(Theta, Data, MP) {
  ngene <- MP$ngene
  J <- MP$J
  M <- MP$M
  Z <- Data$Z
  W <- Theta$W
  P <- Theta$P

  SigInvNorm <- build.siginv(Theta, Data, MP)*Theta$Var$sigsq
  oneterm <- function(i,j) {
    f <- function(i, j) {
      P[i,j]*(Z[i,] - W[,j]) %*%
        SigInvNorm %*% (Z[i,] - W[,j])
    }
    mapply(f,i,j)
  }

  sum(outer(1:ngene, 1:J, oneterm))/(ngene*M)
}

update.rho <- function(Theta, Data, MP) {
  P <- Theta$P
  omega <- Theta$omega
  ngene <- MP$ngene
  rho <- Theta$Var$rho
  sigsq <- Theta$Var$sigsq
  J <- MP$J
  M <- MP$M
  Z <- Data$Z
  W <- Theta$W

  build.d1mat <- function(rho) {
    a1 <- (2*rho/(1-rho^2)^2)*(1+rho^2)+2*rho/(1-rho^2)
    a2 <- (-2*rho^2/(1-rho^2)^2) - 1/(1-rho^2)
    tmp <- toeplitz(c(a1, a2, rep(0, M-2)))
    tmp[1,1] <- tmp[M,M] <- 2*rho/(1-rho^2)^2
    tmp
  }

  build.d2mat <- function(rho) {
    m1 <- toeplitz(c(1+rho^2, -rho, rep(0, M-2)))
    m1[1,1] <- m1[M,M] <- 1
    m2 <- toeplitz(c(2*rho, -1, rep(0, M-2)))
    m2[1,1] <- m2[M,M] <- 0
    m3 <- diag(c(0, rep(2, M-2), 0))

    ((2*(1-rho^2)+8*rho^2)/(1-rho^2)^3)*m1 +
      ((4*rho)/(1-rho^2)^2)*m2 + (1/(1-rho^2))*m3
  }

  d1 <- function(rho) {
    d1mat <- build.d1mat(rho)
    oneterm <- function(i,j) {
      f <- function(i,j) {
        P[i,j]*((M-1)*rho/(1-rho^2) -
                  (1/(2*sigsq)) *	(Z[i,] - W[,j]) %*% d1mat %*%
                  (Z[i,] - W[,j]))
      }
      mapply(f, i, j)
    }

    sum(outer(1:ngene, 1:J, oneterm))
  }

  d2 <- function(rho) {
    d2mat <- build.d2mat(rho)
    oneterm <- function(i,j) {
      f <- function(i,j) {
        P[i,j]*((M-1)*(1+rho^2)/(1-rho^2)^2 -
                  (1/(2*sigsq)) *	(Z[i,] - W[,j]) %*% d2mat %*%
                  (Z[i,] - W[,j]))
      }
      mapply(f, i, j)
    }

    sum(outer(1:ngene, 1:J, oneterm))
  }

  rho.new <- rho - d1(rho)/d2(rho)
  if(abs(rho.new) >= 1) sign(rho.new)*.5
  else rho.new
}

continue.em <- function(Theta.old, Theta, epsConv)
{
  delta <- max(abs(
    unlist(Theta.old[-match("ll", names(Theta.old))]) -
      unlist(Theta[-match("ll", names(Theta))])))
  cat("delta =", delta, "   epsConv =", epsConv, "\n")


  if ( delta < epsConv ) return(FALSE)
  else return(TRUE)
}

BiglogL <- function(Theta, Data, MP)
{
  P <- Theta$P
  omega <- Theta$omega
  J <- MP$J
  K <- MP$K
  Z <- Data$Z
  W <- Theta$W
  M <- MP$M

  ngene <- MP$ngene
  Sig <- build.sig(Theta, Data, MP)
  SigInv <- build.siginv(Theta, Data, MP)

  oneterm <- function(i, j) {
    f <- function(i,j)
      P[i,j]*(log(omega[j]) - M*log(2*pi)/2 - log(det(Sig))/2 - (Z[i,] - W[,j]) %*% SigInv %*% (Z[i,] - W[,j])/2)
    mapply(f, i, j)
  }
  tmp <- outer(1:ngene, 1:J, oneterm)
  tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
  sum(tmp)
}
#-----------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

runEM <- function(Theta, Data, MP, CP)
{
  max.iter <- CP$max.iter
  eps.conv <- CP$eps.conv

  CONT <- TRUE
  iter <- 0

  while(CONT)
  {
    iter <- iter + 1

    Theta.old <- Theta
    Theta <- ee(Theta, Data, MP)
    Theta <- em(Theta, Data, MP, CP, iter)

    Theta$ll <- BiglogL(Theta, Data, MP)

    if(CP$print.updates) {
      cat("Theta:\n")
      print(Theta$omega)
      print(Theta$Mu)
      print(Theta$Var)
    }

    CONT <- continue.em(Theta.old, Theta, eps.conv)
    if(iter >= max.iter) CONT <- FALSE

    if(CP$print.updates) {
      cat("Log Likelihood:", Theta$ll, "\n")
      cat("iteration = ", iter, "\n\n\n")
    }
  }

  return(list(Theta = Theta, Data = Data))
}

ee <- function(Theta, Data, MP)
{
  Theta$P <- update.p(Theta, Data, MP)
  return(Theta)
}

em <- function(Theta, Data, MP, CP, iter)
{
  Theta$omega <- update.omega(Theta)
  Theta$W <- update.W(Theta, Data, MP)
  Theta$Var$sigsq <- update.sigsq(Theta,Data, MP)

  if( iter > CP$ar.skip )
    Theta$Var$rho <- update.rho(Theta, Data, MP)

  return(Theta)
}
calculate_groups <- function(tmp_data, mean_matrix, p_matrix)
{
  col_dim <- dim(tmp_data)[2]
  belong_matrix <- matrix(data = NA, nrow = dim(tmp_data)[1], ncol = 2)
  belong_matrix[,1] <- apply(p_matrix, 1, FUN = function(x){return(which.max(x)-1)})
  belong_matrix[,2] <- apply(tmp_data, 1, FUN = function(x){return(which.min(dist(rbind(x,t(mean_matrix)))[1:col_dim])-1)})
  return(belong_matrix)
}

init.n.run <- function(tmp_data,nt, J=9, r=0, init_rho, init_sigsq, out_file_name)
{
  #数据预处理
  #tmp <- read.table(datafile, skip=2, row.names=1, na.string = "null")
  tmp <- tmp_data
  tmp <- as.matrix(tmp)
  colnames(tmp) <- NULL
  keep.rows <- apply(tmp, 1, function(rw) !any(is.na(rw)))  #找到有缺失数据的标记
  dat <- tmp[keep.rows,]   #去除缺失数据
  ngene <- nrow(dat)       #去除缺失后的标记的总个数
  times <- nt
  if(length(times) %% 2^r != 0) stop("Dimension Reduction too great")
  #整合算法需要的数据
  Data <- list(Y = dat, Z = NULL, tm = times)
  MP <- list(J=J, r=r, ngene=ngene, M=length(times)/(2^r))
  CP <- list(max.iter=100, eps.conv=.2, ar.skip=3, print.updates=TRUE)  #控制信息

  if(r==0){
    Data$Z <- dat
  }
  else {
    Data$Z <- filter.y(Data, MP)
  }
  #给定初始值
  om.init <- c(rep(1/(J), J))
  P.init <- runif(ngene*J)
  dim(P.init) <- c(ngene, J)
  for(i in 1:ngene) P.init[i,] <- P.init[i,]/sum(P.init[i,])

  W <- matrix(rnorm(MP$M*MP$J, mean=sqrt(2)^r), nrow = MP$M, ncol=J)
  W[,J] <- sqrt(2)^r

  rho <- init_rho
  Theta <- list(
    omega = om.init,
    P = P.init,
    W=W,
    Var = list(	rho= init_rho,
                sigsq = init_sigsq),
    ll = -Inf
  )
  #runEM(Theta, Data, MP, CP)
  mod <- runEM(Theta, Data, MP, CP)
  belong_matrix <- calculate_groups(tmp_data = dat, mean_matrix = mod$Theta$W, p_matrix = mod$Theta$P)
  write.table(belong_matrix, file = paste0(out_file_name, "_", J, "_", "Belong_matrix.txt"), row.names = FALSE, col.names = FALSE)
  write.table(mod$Theta$W,   file = paste0(out_file_name, "_", J, "_", "Means_matrix.txt"),  row.names = FALSE, col.names = FALSE)
  write.table(mod$Theta$ll,  file = paste0(out_file_name, "_", "log_likehood.txt"),  row.names = FALSE, col.names = FALSE)
}
