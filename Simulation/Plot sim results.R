
rm(list=ls()) 

# set which results to use
result_gibbs <- readRDS("results_37.rds")
result_approx <- readRDS("results_85.rds")

pdf(file = "/Users/sjoerd/Desktop/Simulation/Figures/comb37.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

library(ggplot2)

##################################################################################################
############################################## Plot ##############################################
##################################################################################################

FPR_gibbs <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
TPR_gibbs <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
FPR_approx <- matrix(NA, nrow = nrow(result_approx[[1]]$lambdas), ncol = 25)
TPR_approx <- matrix(NA, nrow = nrow(result_approx[[1]]$lambdas), ncol = 25)
FPR_jgl <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
TPR_jgl <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
FPR_glasso <- matrix(NA, nrow = length(unique(result_gibbs[[1]]$lambdas[,1])), ncol = 25)
TPR_glasso <- matrix(NA, nrow = length(unique(result_gibbs[[1]]$lambdas[,1])), ncol = 25)
for(i in 1:25){
  FPR_gibbs[,i] <- result_gibbs[[i]]$FPR
  TPR_gibbs[,i] <- result_gibbs[[i]]$TPR
  FPR_approx[,i] <- result_approx[[i]]$FPR
  TPR_approx[,i] <- result_approx[[i]]$TPR
  FPR_jgl[,i] <- result_gibbs[[i]]$FPR_jgl
  TPR_jgl[,i] <- result_gibbs[[i]]$TPR_jgl
  FPR_glasso[,i] <- result_gibbs[[i]]$FPR_glasso
  TPR_glasso[,i] <- result_gibbs[[i]]$TPR_glasso
}
fpr_avg_gibbs <- rowMeans(FPR_gibbs)
tpr_avg_gibbs <- rowMeans(TPR_gibbs)
fpr_avg_approx <- rowMeans(FPR_approx)
tpr_avg_approx <- rowMeans(TPR_approx)
fpr_avg_jgl <- rowMeans(FPR_jgl)
tpr_avg_jgl <- rowMeans(TPR_jgl)
fpr_avg_glasso <- rowMeans(FPR_glasso)
tpr_avg_glasso <- rowMeans(TPR_glasso)

roc_data <- cbind(fpr_avg_gibbs, tpr_avg_gibbs, result_gibbs[[1]]$lambdas[,2])
roc_data <- as.data.frame(roc_data)
colnames(roc_data) <- c("FPR", "TPR", "l2")

fpr_avg_approx <- c(fpr_avg_approx, 0, 0, 0)
tpr_avg_approx <- c(tpr_avg_approx, 0, 0, 0)
lambdas2 <- c(result_approx[[1]]$lambdas[,2], 0, 0.1, 1.0)
roc_data_approx <- cbind(fpr_avg_approx, tpr_avg_approx, lambdas2)
roc_data_approx <- as.data.frame(roc_data_approx)
colnames(roc_data_approx) <- c("FPR", "TPR", "l2")
roc_data_approx[which(roc_data_approx$l2 == 0),"l2"] <- 2.01
roc_data_approx[which(roc_data_approx$l2 == 0.1),"l2"] <- 2.11
roc_data_approx[which(roc_data_approx$l2 == 1),"l2"] <- 3.01

fpr_avg_jgl <- c(fpr_avg_jgl, 0, 0, 0)
tpr_avg_jgl <- c(tpr_avg_jgl, 0, 0, 0)
lambdas2 <- c(result_gibbs[[1]]$lambdas[,2], 0, 0.1, 1.0)
roc_data_jgl <- cbind(fpr_avg_jgl, tpr_avg_jgl, lambdas2)
roc_data_jgl <- as.data.frame(roc_data_jgl)
colnames(roc_data_jgl) <- c("FPR", "TPR", "l2")
roc_data_jgl[which(roc_data_jgl$l2 == 0),"l2"] <- 0.01
roc_data_jgl[which(roc_data_jgl$l2 == 0.1),"l2"] <- 0.11
roc_data_jgl[which(roc_data_jgl$l2 == 1),"l2"] <- 1.01


fpr_avg_glasso <- c(fpr_avg_glasso, 1)
tpr_avg_glasso <- c(tpr_avg_glasso, 1)
roc_data_glasso <- cbind(fpr_avg_glasso, tpr_avg_glasso, rep(2, (length(fpr_avg_glasso))))
roc_data_glasso <- as.data.frame(roc_data_glasso)
colnames(roc_data_glasso) <- c("FPR", "TPR", "l2")

roc_data <- rbind(roc_data, roc_data_approx, roc_data_jgl, roc_data_glasso)
roc_data$l2 <- as.factor(roc_data$l2)
levels(roc_data$l2)
levels(roc_data$l2) <- c("Gibbs 0", "FGL 0", "Gibbs 0.1", "FGL 0.1", "Gibbs 1", "FGL 1", "GLASSO", "Approx 0", "Approx 0.1", "Approx1")
colnames(roc_data) <- c("FPR", "TPR", "Method")

ggplot(data=roc_data,
       aes(x=FPR, y=TPR, colour=Method, linetype = Method)) +
  geom_line() +
  #geom_abline(intercept = 0, slope = 1)+
  #annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1)+
  theme_classic()+
  scale_linetype_manual(values=c("solid", "solid", "dotdash", "dotdash", "dotted", "dotted", "solid", "solid", "dotdash", "dotted"))+
  scale_color_manual(values=c("#900C3F", "#2471A3", "#900C3F", "#2471A3", "#900C3F", "#2471A3", "#229954", "#FFC300", "#FFC300", "#FFC300"))+
  theme(legend.position = "none")

dev.off()


##################################################################################################
############################################## AUC ############################################### 
##################################################################################################

rm(list=ls()) 
setwd("/Users/sjoerd/Desktop/Simulation/")
result_gibbs <- readRDS("results_37.rds")
result_approx <- readRDS("results_85.rds")

FPR_gibbs <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
TPR_gibbs <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
FPR_approx <- matrix(NA, nrow = nrow(result_approx[[1]]$lambdas), ncol = 25)
TPR_approx <- matrix(NA, nrow = nrow(result_approx[[1]]$lambdas), ncol = 25)
FPR_jgl <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
TPR_jgl <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
FPR_glasso <- matrix(NA, nrow = length(unique(result_gibbs[[1]]$lambdas[,1])), ncol = 25)
TPR_glasso <- matrix(NA, nrow = length(unique(result_gibbs[[1]]$lambdas[,1])), ncol = 25)
for(i in 1:25){
  FPR_gibbs[,i] <- result_gibbs[[i]]$FPR
  TPR_gibbs[,i] <- result_gibbs[[i]]$TPR
  FPR_approx[,i] <- result_approx[[i]]$FPR
  TPR_approx[,i] <- result_approx[[i]]$TPR
  FPR_jgl[,i] <- result_gibbs[[i]]$FPR_jgl
  TPR_jgl[,i] <- result_gibbs[[i]]$TPR_jgl
  FPR_glasso[,i] <- result_gibbs[[i]]$FPR_glasso
  TPR_glasso[,i] <- result_gibbs[[i]]$TPR_glasso
}
fpr_avg_gibbs <- rowMeans(FPR_gibbs)
tpr_avg_gibbs <- rowMeans(TPR_gibbs)
fpr_avg_approx <- rowMeans(FPR_approx)
tpr_avg_approx <- rowMeans(TPR_approx)
fpr_avg_jgl <- rowMeans(FPR_jgl)
tpr_avg_jgl <- rowMeans(TPR_jgl)
fpr_avg_glasso <- rowMeans(FPR_glasso)
tpr_avg_glasso <- rowMeans(TPR_glasso)

# to correct for partial ROC curves, we need to add these values
# they do not change the AUC of proper curves
fpr_avg_approx <- c(fpr_avg_approx, 0, 0, 0)
tpr_avg_approx <- c(tpr_avg_approx, 0, 0, 0)
fpr_avg_jgl <- c(fpr_avg_jgl, 0, 0, 0)
tpr_avg_jgl <- c(tpr_avg_jgl, 0, 0, 0)
fpr_avg_glasso <- c(fpr_avg_glasso, 1)
tpr_avg_glasso <- c(tpr_avg_glasso, 1)


library(DescTools)

# results 1:21 are lambda2 = 0
# results 22:42 are lambda2 = 0.1
# results 43:63 are lambda2 = 1
auc_gibbs_0 <- AUC(x = fpr_avg_gibbs[1:21], y = tpr_avg_gibbs[1:21])
auc_gibbs_01 <- AUC(x = fpr_avg_gibbs[22:42], y = tpr_avg_gibbs[22:42])
auc_gibbs_1 <- AUC(x = fpr_avg_gibbs[43:63], y = tpr_avg_gibbs[43:63])
auc_gibbs_avg <- mean(c(auc_gibbs_0, auc_gibbs_01, auc_gibbs_1))

auc_approx_0 <- AUC(x = fpr_avg_approx[1:21], y = tpr_avg_approx[1:21])
auc_approx_01 <- AUC(x = fpr_avg_approx[22:42], y = tpr_avg_approx[22:42])
auc_approx_1 <- AUC(x = fpr_avg_approx[43:63], y = tpr_avg_approx[43:63])
auc_approx_avg <- mean(c(auc_approx_0, auc_approx_01, auc_approx_1))

auc_jgl_0 <- AUC(x = fpr_avg_jgl[1:21], y = tpr_avg_jgl[1:21])
auc_jgl_01 <- AUC(x = fpr_avg_jgl[22:42], y = tpr_avg_jgl[22:42])
auc_jgl_1 <- AUC(x = fpr_avg_jgl[43:63], y = tpr_avg_jgl[43:63])
auc_jgl_avg <- mean(c(auc_jgl_0, auc_jgl_01, auc_jgl_1))

auc_glasso <- AUC(x = fpr_avg_glasso, y = tpr_avg_glasso)

auc_gibbs_bc <- max(auc_gibbs_0, auc_gibbs_01, auc_gibbs_1)
auc_approx_bc <- max(auc_approx_0, auc_approx_01, auc_approx_1)
auc_jgl_bc <- max(auc_jgl_0, auc_jgl_01, auc_jgl_1)

##################################################################################################
############################################## FL ################################################  
##################################################################################################

fl_mean_gibbs0 <- vector("list", 25)
fl_mean_gibbs01 <- vector("list", 25)
fl_mean_gibbs1 <- vector("list", 25)
fl_mean_approx0 <- vector("list", 25)
fl_mean_approx01 <- vector("list", 25)
fl_mean_approx1 <- vector("list", 25)
fl_jgl_mean0 <- vector("list", 25)
fl_jgl_mean01 <- vector("list", 25)
fl_jgl_mean1 <- vector("list", 25)
fl_glasso_mean <- vector("list", 25)

for(i in 1:25){
  # we should remove every 1st number for fairness
  # as the numbers are too high to count
  fl_mean_gibbs0[[i]] <- result_gibbs[[i]]$FL[1:21][-1]
  fl_mean_gibbs01[[i]] <- result_gibbs[[i]]$FL[22:42][-1]
  fl_mean_gibbs1[[i]] <- result_gibbs[[i]]$FL[43:63][-1]
  fl_mean_approx0[[i]] <- result_approx[[i]]$FL[1:21][-1]
  fl_mean_approx01[[i]] <- result_approx[[i]]$FL[22:42][-1]
  fl_mean_approx1[[i]] <- result_approx[[i]]$FL[43:63][-1]
  fl_jgl_mean0[[i]] <- result_approx[[i]]$FL_jgl[1:21][-1]
  fl_jgl_mean01[[i]] <- result_approx[[i]]$FL_jgl[22:42][-1]
  fl_jgl_mean1[[i]] <- result_approx[[i]]$FL_jgl[43:63][-1]
  fl_glasso_mean[[i]] <- result_approx[[i]]$FL_glasso[-1]
}

fl_gibbs0 <- mean(unlist(fl_mean_gibbs0))
fl_gibbs01 <- mean(unlist(fl_mean_gibbs01))
fl_gibbs1 <- mean(unlist(fl_mean_gibbs1))
fl_gibbs_avg <- mean(c(fl_gibbs0, fl_gibbs01, fl_gibbs1))

fl_approx0 <- mean(unlist(fl_mean_approx0))
fl_approx01 <- mean(unlist(fl_mean_approx01))
fl_approx1 <- mean(unlist(fl_mean_approx1))
fl_approx_avg <- mean(c(fl_approx0, fl_approx01, fl_approx1))

fl_jgl_0 <- mean(unlist(fl_jgl_mean0))
fl_jgl_01 <- mean(unlist(fl_jgl_mean01))
fl_jgl_1 <- mean(unlist(fl_jgl_mean1))
fl_jgl_avg <- mean(c(fl_jgl_0, fl_jgl_01, fl_jgl_1))

fl_glasso <- mean(unlist(fl_glasso_mean))

fl_gibbs_bc <- min(fl_gibbs0, fl_gibbs01, fl_gibbs1)
fl_approx_bc <- min(fl_approx0, fl_approx01, fl_approx1)
fl_jgl_bc <- min(fl_jgl_0, fl_jgl_01, fl_jgl_1)


##################################################################################################
############################################## EL ################################################  
##################################################################################################

el_mean_approx0 <- vector("list", 25)
el_mean_approx01 <- vector("list", 25)
el_mean_approx1 <- vector("list", 25)
el_mean_gibbs0 <- vector("list", 25)
el_mean_gibbs01 <- vector("list", 25)
el_mean_gibbs1 <- vector("list", 25)
el_jgl_mean0 <- vector("list", 25)
el_jgl_mean01 <- vector("list", 25)
el_jgl_mean1 <- vector("list", 25)
el_glasso_mean <- vector("list", 25)

for(i in 1:25){
  # we should remove every 1st number for fairness
  # as the numbers are too high to count
  el_mean_approx0[[i]] <- result_approx[[i]]$EL[1:21][-1]
  el_mean_approx01[[i]] <- result_approx[[i]]$EL[22:42][-1]
  el_mean_approx1[[i]] <- result_approx[[i]]$EL[43:63][-1]
  el_mean_gibbs0[[i]] <- result_gibbs[[i]]$EL[1:21][-1]
  el_mean_gibbs01[[i]] <- result_gibbs[[i]]$EL[22:42][-1]
  el_mean_gibbs1[[i]] <- result_gibbs[[i]]$EL[43:63][-1]
  el_jgl_mean0[[i]] <- result_approx[[i]]$EL_jgl[1:21][-1]
  el_jgl_mean01[[i]] <- result_approx[[i]]$EL_jgl[22:42][-1]
  el_jgl_mean1[[i]] <- result_approx[[i]]$EL_jgl[43:63][-1]
  el_glasso_mean[[i]] <- result_approx[[i]]$EL_glasso[-1]
}

el_gibbs0 <- mean(unlist(el_mean_gibbs0))
el_gibbs01 <- mean(unlist(el_mean_gibbs01))
el_gibbs1 <- mean(unlist(el_mean_gibbs1))
el_gibbs_avg <- mean(c(el_gibbs0, el_gibbs01, el_gibbs1))

el_approx0 <- mean(unlist(el_mean_approx0))
el_approx01 <- mean(unlist(el_mean_approx01))
el_approx1 <- mean(unlist(el_mean_approx1))
el_approx_avg <- mean(c(el_approx0, el_approx01, el_approx1))

el_jgl0 <- mean(unlist(el_jgl_mean0))
el_jgl01 <- mean(unlist(el_jgl_mean01))
el_jgl1 <- mean(unlist(el_jgl_mean1))
el_jgl_avg <- mean(c(el_jgl0, el_jgl01, el_jgl1))

el_glasso <- mean(unlist(el_glasso_mean))

el_gibbs_bc <- min(el_gibbs0, el_gibbs01, el_gibbs1)
el_approx_bc <- min(el_approx0, el_approx01, el_approx1)
el_jgl_bc <- min(el_jgl0, el_jgl01, el_jgl1)

##################################################################################################
############################################# Make table #########################################  
##################################################################################################

library(xtable)
temp_tab <- round(c(auc_gibbs_avg, auc_approx_avg, 0, 0, fl_gibbs_avg, fl_approx_avg, el_gibbs_avg,
        el_approx_avg, auc_gibbs_bc, auc_approx_avg, 0, 0, fl_gibbs_bc, fl_approx_bc,
        el_gibbs_bc, el_approx_bc),3)

a <- as.character(c(auc_gibbs_avg, "/", auc_approx_avg, 0, 0, fl_gibbs_avg, fl_approx_avg, el_gibbs_avg,
  el_approx_avg, auc_gibbs_bc, auc_approx_avg, 0, 0, fl_gibbs_bc, fl_approx_bc,
  el_gibbs_bc, el_approx_bc))

temp_tab <- as.data.frame(matrix(temp_tab, nrow = 1))
xtable(temp_tab, digits = 2)


##################################################################################################
############################################## Load many #########################################  
##################################################################################################

gibbsAdress <- data.frame(rep(NA, 48))
gibbsAdress[,1] <- c("results_49.rds",
                     "results_50.rds",
                     "results_51.rds",
                     "results_52.rds",
                     "results_53.rds",
                     "results_54.rds",
                     "results_55.rds",
                     "results_56.rds",
                     "results_57.rds",
                     "results_58.rds",
                     "results_59.rds",
                     "results_60.rds",
                     "results_61.rds",
                     "results_62.rds",
                     "results_63.rds",
                     "results_64.rds",
                     "results_65.rds",
                     "results_66.rds",
                     "results_67.rds",
                     "results_68.rds",
                     "results_69.rds",
                     "results_70.rds",
                     "results_71.rds",
                     "results_72.rds",
                     "results_73.rds",
                     "results_74.rds",
                     "results_75.rds",
                     "results_76.rds",
                     "results_77.rds",
                     "results_78.rds",
                     "results_79.rds",
                     "results_80.rds",
                     "results_81.rds",
                     "results_82.rds",
                     "results_83.rds",
                     "results_84.rds",
                     "results_85.rds",
                     "results_86.rds",
                     "results_87.rds",
                     "results_88.rds",
                     "results_89.rds",
                     "results_90.rds",
                     "results_91.rds",
                     "results_92.rds",
                     "results_93.rds",
                     "results_94.rds",
                     "results_95.rds",
                     "results_96.rds")
approxAdress <- data.frame(rep(NA, 48))
approxAdress[,1] <- c("results_1.rds",
                     "results_2.rds",
                     "results_3.rds",
                     "results_4.rds",
                     "results_5.rds",
                     "results_6.rds",
                     "results_7.rds",
                     "results_8.rds",
                     "results_9.rds",
                     "results_10.rds",
                     "results_11.rds",
                     "results_12.rds",
                     "results_13.rds",
                     "results_14.rds",
                     "results_15.rds",
                     "results_16.rds",
                     "results_17.rds",
                     "results_18.rds",
                     "results_19.rds",
                     "results_20.rds",
                     "results_21.rds",
                     "results_22.rds",
                     "results_23.rds",
                     "results_24.rds",
                     "results_25.rds",
                     "results_26.rds",
                     "results_27.rds",
                     "results_28.rds",
                     "results_29.rds",
                     "results_30.rds",
                     "results_31.rds",
                     "results_32.rds",
                     "results_33.rds",
                     "results_34.rds",
                     "results_35.rds",
                     "results_36.rds",
                     "results_37.rds",
                     "results_38.rds",
                     "results_39.rds",
                     "results_40.rds",
                     "results_41.rds",
                     "results_42.rds",
                     "results_43.rds",
                     "results_44.rds",
                     "results_45.rds",
                     "results_46.rds",
                     "results_47.rds",
                     "results_48.rds")

#ncomb <- 48
ncomb <- 16
library(DescTools)

res_df1 <- matrix(NA, 16, 12)
res_df2 <- matrix(NA, 16, 9)

idx <- 33:48

for(comb in 1:ncomb){
  result_gibbs <- readRDS(gibbsAdress[idx[comb],1])
  result_approx <- readRDS(approxAdress[idx[comb],1])
  
  FPR_gibbs <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
  TPR_gibbs <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
  FPR_approx <- matrix(NA, nrow = nrow(result_approx[[1]]$lambdas), ncol = 25)
  TPR_approx <- matrix(NA, nrow = nrow(result_approx[[1]]$lambdas), ncol = 25)
  FPR_jgl <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
  TPR_jgl <- matrix(NA, nrow = nrow(result_gibbs[[1]]$lambdas), ncol = 25)
  FPR_glasso <- matrix(NA, nrow = length(unique(result_gibbs[[1]]$lambdas[,1])), ncol = 25)
  TPR_glasso <- matrix(NA, nrow = length(unique(result_gibbs[[1]]$lambdas[,1])), ncol = 25)
  for(i in 1:25){
    FPR_gibbs[,i] <- result_gibbs[[i]]$FPR
    TPR_gibbs[,i] <- result_gibbs[[i]]$TPR
    FPR_approx[,i] <- result_approx[[i]]$FPR
    TPR_approx[,i] <- result_approx[[i]]$TPR
    FPR_jgl[,i] <- result_gibbs[[i]]$FPR_jgl
    TPR_jgl[,i] <- result_gibbs[[i]]$TPR_jgl
    FPR_glasso[,i] <- result_gibbs[[i]]$FPR_glasso
    TPR_glasso[,i] <- result_gibbs[[i]]$TPR_glasso
  }
  fpr_avg_gibbs <- rowMeans(FPR_gibbs)
  tpr_avg_gibbs <- rowMeans(TPR_gibbs)
  fpr_avg_approx <- rowMeans(FPR_approx)
  tpr_avg_approx <- rowMeans(TPR_approx)
  fpr_avg_jgl <- rowMeans(FPR_jgl)
  tpr_avg_jgl <- rowMeans(TPR_jgl)
  fpr_avg_glasso <- rowMeans(FPR_glasso)
  tpr_avg_glasso <- rowMeans(TPR_glasso)
  
  # to correct for partial ROC curves, we need to add these values
  # they do not change the AUC of proper curves
  fpr_avg_approx <- c(fpr_avg_approx, 0, 0, 0)
  tpr_avg_approx <- c(tpr_avg_approx, 0, 0, 0)
  fpr_avg_jgl <- c(fpr_avg_jgl, 0, 0, 0)
  tpr_avg_jgl <- c(tpr_avg_jgl, 0, 0, 0)
  fpr_avg_glasso <- c(fpr_avg_glasso, 1)
  tpr_avg_glasso <- c(tpr_avg_glasso, 1)

  auc_gibbs_0 <- AUC(x = fpr_avg_gibbs[1:21], y = tpr_avg_gibbs[1:21])
  auc_gibbs_01 <- AUC(x = fpr_avg_gibbs[22:42], y = tpr_avg_gibbs[22:42])
  auc_gibbs_1 <- AUC(x = fpr_avg_gibbs[43:63], y = tpr_avg_gibbs[43:63])
  auc_gibbs_avg <- mean(c(auc_gibbs_0, auc_gibbs_01, auc_gibbs_1))
  
  auc_approx_0 <- AUC(x = fpr_avg_approx[1:21], y = tpr_avg_approx[1:21])
  auc_approx_01 <- AUC(x = fpr_avg_approx[22:42], y = tpr_avg_approx[22:42])
  auc_approx_1 <- AUC(x = fpr_avg_approx[43:63], y = tpr_avg_approx[43:63])
  auc_approx_avg <- mean(c(auc_approx_0, auc_approx_01, auc_approx_1))
  
  auc_jgl_0 <- AUC(x = fpr_avg_jgl[1:21], y = tpr_avg_jgl[1:21])
  auc_jgl_01 <- AUC(x = fpr_avg_jgl[22:42], y = tpr_avg_jgl[22:42])
  auc_jgl_1 <- AUC(x = fpr_avg_jgl[43:63], y = tpr_avg_jgl[43:63])
  auc_jgl_avg <- mean(c(auc_jgl_0, auc_jgl_01, auc_jgl_1))
  
  auc_glasso <- AUC(x = fpr_avg_glasso, y = tpr_avg_glasso)
  
  auc_gibbs_bc <- max(c(auc_gibbs_0, auc_gibbs_01, auc_gibbs_1))
  auc_approx_bc <- max(c(auc_approx_0, auc_approx_01, auc_approx_1))
  auc_jgl_bc <- max(c(auc_jgl_0, auc_jgl_01, auc_jgl_1))
  
  fl_mean_gibbs0 <- vector("list", 25)
  fl_mean_gibbs01 <- vector("list", 25)
  fl_mean_gibbs1 <- vector("list", 25)
  fl_mean_approx0 <- vector("list", 25)
  fl_mean_approx01 <- vector("list", 25)
  fl_mean_approx1 <- vector("list", 25)
  fl_jgl_mean0 <- vector("list", 25)
  fl_jgl_mean01 <- vector("list", 25)
  fl_jgl_mean1 <- vector("list", 25)
  fl_glasso_mean <- vector("list", 25)
  
  for(i in 1:25){
    # we should remove every 1st number for fairness
    # as the numbers are too high to count
    fl_mean_gibbs0[[i]] <- result_gibbs[[i]]$FL[1:21][-1]
    fl_mean_gibbs01[[i]] <- result_gibbs[[i]]$FL[22:42][-1]
    fl_mean_gibbs1[[i]] <- result_gibbs[[i]]$FL[43:63][-1]
    fl_mean_approx0[[i]] <- result_approx[[i]]$FL[1:21][-1]
    fl_mean_approx01[[i]] <- result_approx[[i]]$FL[22:42][-1]
    fl_mean_approx1[[i]] <- result_approx[[i]]$FL[43:63][-1]
    fl_jgl_mean0[[i]] <- result_approx[[i]]$FL_jgl[1:21][-1]
    fl_jgl_mean01[[i]] <- result_approx[[i]]$FL_jgl[22:42][-1]
    fl_jgl_mean1[[i]] <- result_approx[[i]]$FL_jgl[43:63][-1]
    fl_glasso_mean[[i]] <- result_approx[[i]]$FL_glasso[-1]
  }
  
  fl_gibbs0 <- mean(unlist(fl_mean_gibbs0))
  fl_gibbs01 <- mean(unlist(fl_mean_gibbs01))
  fl_gibbs1 <- mean(unlist(fl_mean_gibbs1))
  fl_gibbs_avg <- mean(c(fl_gibbs0, fl_gibbs01, fl_gibbs1))
  
  fl_approx0 <- mean(unlist(fl_mean_approx0))
  fl_approx01 <- mean(unlist(fl_mean_approx01))
  fl_approx1 <- mean(unlist(fl_mean_approx1))
  fl_approx_avg <- mean(c(fl_approx0, fl_approx01, fl_approx1))
  
  fl_jgl_0 <- mean(unlist(fl_jgl_mean0))
  fl_jgl_01 <- mean(unlist(fl_jgl_mean01))
  fl_jgl_1 <- mean(unlist(fl_jgl_mean1))
  fl_jgl_avg <- mean(c(fl_jgl_0, fl_jgl_01, fl_jgl_1))
  
  fl_glasso <- mean(unlist(fl_glasso_mean))
  
  fl_gibbs_bc <- min(c(fl_gibbs0, fl_gibbs01, fl_gibbs1))
  fl_approx_bc <- min(c(fl_approx0, fl_approx01, fl_approx1))
  fl_jgl_bc <- min(c(fl_jgl_0, fl_jgl_01, fl_jgl_1))
  
  el_mean_approx0 <- vector("list", 25)
  el_mean_approx01 <- vector("list", 25)
  el_mean_approx1 <- vector("list", 25)
  el_mean_gibbs0 <- vector("list", 25)
  el_mean_gibbs01 <- vector("list", 25)
  el_mean_gibbs1 <- vector("list", 25)
  el_jgl_mean0 <- vector("list", 25)
  el_jgl_mean01 <- vector("list", 25)
  el_jgl_mean1 <- vector("list", 25)
  el_glasso_mean <- vector("list", 25)
  
  for(i in 1:25){
    # we should remove every 1st number for fairness
    # as the numbers are too high to count
    el_mean_approx0[[i]] <- result_approx[[i]]$EL[1:21][-1]
    el_mean_approx01[[i]] <- result_approx[[i]]$EL[22:42][-1]
    el_mean_approx1[[i]] <- result_approx[[i]]$EL[43:63][-1]
    el_mean_gibbs0[[i]] <- result_gibbs[[i]]$EL[1:21][-1]
    el_mean_gibbs01[[i]] <- result_gibbs[[i]]$EL[22:42][-1]
    el_mean_gibbs1[[i]] <- result_gibbs[[i]]$EL[43:63][-1]
    el_jgl_mean0[[i]] <- result_approx[[i]]$EL_jgl[1:21][-1]
    el_jgl_mean01[[i]] <- result_approx[[i]]$EL_jgl[22:42][-1]
    el_jgl_mean1[[i]] <- result_approx[[i]]$EL_jgl[43:63][-1]
    el_glasso_mean[[i]] <- result_approx[[i]]$EL_glasso[-1]
  }
  
  el_gibbs0 <- mean(unlist(el_mean_gibbs0))
  el_gibbs01 <- mean(unlist(el_mean_gibbs01))
  el_gibbs1 <- mean(unlist(el_mean_gibbs1))
  el_gibbs_avg <- mean(c(el_gibbs0, el_gibbs01, el_gibbs1))
  
  el_approx0 <- mean(unlist(el_mean_approx0))
  el_approx01 <- mean(unlist(el_mean_approx01))
  el_approx1 <- mean(unlist(el_mean_approx1))
  el_approx_avg <- mean(c(el_approx0, el_approx01, el_approx1))
  
  el_jgl0 <- mean(unlist(el_jgl_mean0))
  el_jgl01 <- mean(unlist(el_jgl_mean01))
  el_jgl1 <- mean(unlist(el_jgl_mean1))
  el_jgl_avg <- mean(c(el_jgl0, el_jgl01, el_jgl1))
  
  el_glasso <- mean(unlist(el_glasso_mean))
  
  el_gibbs_bc <- min(c(el_gibbs0, el_gibbs01, el_gibbs1))
  el_approx_bc <- min(c(el_approx0, el_approx01, el_approx1))
  el_jgl_bc <- min(c(el_jgl0, el_jgl01, el_jgl1))

  res_df1[comb,] <- as.numeric((round(c(auc_gibbs_avg, auc_approx_avg, fl_gibbs_avg, fl_approx_avg, el_gibbs_avg,
          el_approx_avg, auc_gibbs_bc, auc_approx_bc, fl_gibbs_bc, fl_approx_bc,
          el_gibbs_bc, el_approx_bc),3)))
  
  res_df2[comb,] <- c(round(c(auc_jgl_avg, auc_glasso, fl_jgl_avg, fl_glasso, el_jgl_avg,
                             el_glasso, auc_jgl_bc, fl_jgl_bc,
                             el_jgl_bc),3))
}



library(xtable)

xtable(res_df1, digits = 2)
xtable(res_df2, digits = 2)



