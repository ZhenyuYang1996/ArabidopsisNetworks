# Title     : Arabidopsis epistatic networks
# Objective : TODO
# Created by: YangZhenyu
setwd("")
rm(list = ls())
#setwd("/home/zyyang/Project/arabidopsis_remote")
#setwd("F:/yang/Project/arabidopsis_remote")
library(parallel)
library(mvtnorm)
library(glmnet)
source("code/sup.R")
source("code/funmap.R")
source("code/funCluster.R")
source("code/lasso_code.R")
source("code/network_sup.R")
source("code/pretreatment.R")
#=====================Setting global Parameters======================
times  <- 1:8
#functional cluster part
cl_num      <- 10
cluster_dim <- 8
#sub functional cluster part
cluster_J       <- 1
sub_cl_num      <- 10
sub_cluster_dim <- 30
nt=seq(times[1], times[length(times)],length = sub_cluster_dim)
#network part
sub_cluster_J  <- 2
legendre_order <- 6
#=====================data======================
pheno_data_location <- "input/Plant_height_data.txt"
snp_data_location   <- "input/snp_data.txt"
snp_info_location   <- "input/snp_info.txt"
funmap_par_location     <- "data/funmap_all_par.txt"
gene_var_file_location  <- "data/gene_var.txt"
out_file_location       <- "data/gene_var_C"
cluster_belong_file     <- paste0(out_file_location, "_", cl_num, "_", "Belong_matrix.txt")
cluster_J_snp_file      <- paste0("data/cluster_",cluster_J,"_snp_info.txt")
cluster_J_eff_file      <- paste0("data/cluster_",cluster_J,"_var_data.txt")
sub_cluster_file        <- paste0("data/sub_cluster_",cluster_J,"_var_C")
sub_cluster_belong_file <- paste0("data/sub_cluster_",cluster_J,"_var_C_",sub_cl_num,"_Belong_matrix.txt")
J_sub_J_gene_eff_file   <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_gene_eff.txt")
J_sub_J_snp_info_file   <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_snp_info.txt")
J_sub_J_lasso_file      <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_connect_mat.txt")
J_sub_J_net_para_file   <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_network_para.txt")
odee_file               <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_odee_res.Rdata")
net_res_file            <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_net_res.Rdata")
net_plot_file           <- paste0("data/sub_cluster_",cluster_J,"_",sub_cluster_J,"_net_plot_file.txt")
#——————————————————————————————Functional mapping————————————————————————————————
#==================Calculate LR value by parallel====================
pheno_data <- read.table(pheno_data_location)
snp_data   <- read.table(snp_data_location)
core.number <- detectCores() #
cl          <- makeCluster(getOption("cl.cores", core.number))
clusterExport(cl,c("dmvnorm", "get_miu","get_sad_continuous","H0_hypothesis" ,"H1_hypothesis","get_LR", "pheno_data", "snp_data"))
all_par     <- parSapply(cl,1:dim(snp_data)[1],function(c){return(get_LR(c))})
stopCluster(cl)
write.table(t(all_par), file = funmap_par_location, row.names = FALSE, col.names = FALSE)
#=====================Calculating genetic variance====================
funmap_par <- read.table(file = funmap_par_location)
gene_eff   <- get_VG(FunMap_par = as.matrix(funmap_par), marker_data = as.matrix(snp_data), t = times)
write.table(gene_eff, file = gene_var_file_location, row.names = FALSE, col.names = FALSE)
#————————————————————————————Functional cluster——————————————————————————————
#==================Run the functional cluster on the Linux/Windows===================
#R version
init.n.run(tmp_data = read.table(file = gene_var_file_location), nt = times, J = cl_num, r = 0, init_rho = 0.1, init_sigsq = 0.1, out_file_name = out_file_location)
#----------------------
#C++ version on Windows
#windows_command <- paste0("code/modfiy_cluster.exe ", cl_num, " ", cl_num, " ", cluster_dim, " ", gene_var_file_location, " ", out_file_location)
#system(command = windows_command)
#----------------------
#C++ version on linux
#linux_command   <- paste0("code/icc-Funcluster ", K_f, " ", K_t, " ", cluster_dim," ", cluster_J_eff_file, " ", sub_cluster_file)
#system(command = linux_command)
#====================================================================
snp_data      <- read.table(file = snp_data_location)
snp_info      <- read.table(file = snp_info_location)
funmap_par    <- read.table(file = funmap_par_location)
belong_data   <- read.table(file = cluster_belong_file)
cluster_J_snp_index <- which( belong_data[,2] == (cluster_J-1) )
if(length(cluster_J_snp_index) < sub_cl_num)
{
  write("Error: The number of sub-clusters is more than the number of samples", stderr())
}
cluster_J_gene_eff <- get_VG(FunMap_par = funmap_par[cluster_J_snp_index,],
                             marker_data = snp_data[cluster_J_snp_index,],
                             t = nt)
cluster_J_snp_info      <- matrix(0, nrow = length(cluster_J_snp_index), ncol = 3)
cluster_J_snp_info[,1]  <- rownames(snp_info[cluster_J_snp_index,])
cluster_J_snp_info[,2]  <- cluster_J_snp_index
for (i in 1:length(cluster_J_snp_index))
{
  cluster_J_snp_info[i,3] <- substring(snp_info[cluster_J_snp_index[i],5],
                                       gregexpr(":A", snp_info[cluster_J_snp_index[i],5])[[1]][1]+1)
}
write.table(cluster_J_snp_info, file = cluster_J_snp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cluster_J_gene_eff, file = cluster_J_eff_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
#————————————————————————————sub Functional cluster——————————————————————————————
#R version
init.n.run(tmp_data = read.table(file = cluster_J_eff_file), nt = nt, J = sub_cl_num, r = 0, init_rho = 0.6, init_sigsq = 0.8, out_file_name = sub_cluster_file)
#----------------------
#C++ version on Windows
#windows_command <- paste0("code/sub_cluster.exe ", sub_cl_num, " ", sub_cl_num, " ", sub_cluster_dim, " ", cluster_J_eff_file, " ", sub_cluster_file)
#system(command = windows_command)
#----------------------
#C++ version on linux
#linux_command   <- paste0("code/icc-Funcluster ", K_f, " ", K_t, " ", cluster_dim," ", cluster_J_eff_file, " ", sub_cluster_file)
#system(command = linux_command)
#
cluster_J_belong_data <- read.table(sub_cluster_belong_file)
cluster_J_snp_info    <- read.table(cluster_J_snp_file)
cluster_J_gene_eff <- read.table(file = cluster_J_eff_file)
cluster_J_sub_J_gene_eff  <- cluster_J_gene_eff[which(cluster_J_belong_data[,2] == (sub_cluster_J-1)), ]
cluster_J_sub_J_snp_info  <- cluster_J_snp_info[which(cluster_J_belong_data[,2] == (sub_cluster_J-1)), ]
write.table(cluster_J_sub_J_gene_eff, file = J_sub_J_gene_eff_file, row.names = FALSE, col.names = FALSE)
write.table(cluster_J_sub_J_snp_info, file = J_sub_J_snp_info_file, row.names = FALSE, col.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------------------
ori_data <- read.table(file = J_sub_J_gene_eff_file)
core.number <- detectCores()
cl <- makeCluster(getOption("cl.cores", core.number))
clusterExport(cl,c("ori_data", "parallel_lasso_new", "cv.glmnet", "glmnet"))
connect_matrix <- parSapply(cl,1:dim(ori_data)[1],function(c){return(parallel_lasso_new(c))})
stopCluster(cl)
write.table(t(connect_matrix), file = J_sub_J_lasso_file, row.names = FALSE, col.names = FALSE)
#---------------------network--------------------
#R version
DBH.odee <- optim.parallel(connect=read.table(file = J_sub_J_lasso_file),effect=t(read.table(J_sub_J_gene_eff_file)),
                           n.cores=1,proc=ode.optim,order=legendre_order,times=nt,nstep=sub_cluster_dim-1)
write.table(t(DBH.odee[[1]]), file = J_sub_J_net_para_file)
#----------------------
#C++ version on Windows
#network_command = paste0("code/network_yang.exe ", legendre_order, " ", sub_cluster_dim, " ", max(times), " ", J_sub_J_gene_eff_file, " ", J_sub_J_lasso_file, " ", J_sub_J_net_para_file)
#system(command = network_command)
#----------------------
#C++ version on linux
#network_command = paste0("code/icc-network ", 6, " ", length(nt), " ", max(times), " ", J_sub_J_gene_eff_file, " ", J_sub_J_lasso_file, " ", J_sub_J_net_para_file)
#system(command = network_command)
#——————————————————————————————————————————————————————————————————————————
network_para              <- read.table(J_sub_J_net_para_file)
cluster_J_sub_J_gene_eff  <- read.table(J_sub_J_gene_eff_file)
connect_matrix            <- read.table(J_sub_J_lasso_file)
#
rownames(cluster_J_sub_J_gene_eff) <- 1:dim(cluster_J_sub_J_gene_eff)[1]
net_odee          <- network_eff(network_para = network_para, network_connect = connect_matrix, effect = t(cluster_J_sub_J_gene_eff), order = 6, times = 1:8, nstep = 29)
#
net_stat          <- interType(con = connect_matrix, alle = net_odee, sme = cluster_J_sub_J_gene_eff)
save(net_odee, file = odee_file, version = 2)
save(net_stat, file = net_res_file, version = 2)
#
load(file = odee_file)
load(file = net_res_file)
cluster_J_sub_J_snp_info <- read.table(file = J_sub_J_snp_info_file)
net_plot_data <- get_network_info(connfest = net_stat[[1]], nodes_info = cluster_J_sub_J_snp_info[,3])
write.table(net_plot_data, file = net_plot_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
