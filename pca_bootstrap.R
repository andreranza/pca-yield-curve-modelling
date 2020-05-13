# PACKAGES
rm(list = ls())
library(xts)
library(ggplot2)
library(lubridate)
library(grid)
library(stringr)
library(RColorBrewer)
cols = brewer.pal(n = 3,name = "Set1")
library(highcharter)

# DATA LOADING
Ts <- read.table("/Users/Andrea/Desktop/thesis/data/USTREASURY_YIELD_UPDATED.csv")
Ts <- xts(Ts,ymd(rownames(Ts)))
Ts <- window(Ts, start="2006-02-09",end="2020-01-29")

Mat_labs <- colnames(Ts) # maturities
PCs_labs <- rep("PC",dim(Ts)[2]) # PCs
for(i in 1:dim(Ts)[2]){
  PCs_labs[i] <- str_c(PCs_labs[i], as.character(i), sep="")
}

# INTEREST RATES CHANGES IN BASIS POINTS
X <- apply(Ts,2,diff)*100
dates <- index(Ts)[-1]
X <- xts(X,dates)

# PCA CORRELATION MATRIX
CorX <- cor(X)
W <- eigen(CorX)$vectors
Z <- X%*%W*(-1)
Z <- xts(Z,dates) 
colnames(Z) <- PCs_labs

# GRAPH FIRST 3 PCS
highchart(type = "stock") %>% 
  hc_title(text = "First three PCs - Corr. Matrix - Absolute rates changes in Bps") %>% 
  hc_add_series(Z[, PCs_labs[1]], name = PCs_labs[1]) %>% 
  hc_add_series(Z[, PCs_labs[2]], name = PCs_labs[2]) %>%   
  hc_add_series(Z[, PCs_labs[3]], name = PCs_labs[3]) %>%
  hc_add_theme(hc_theme_flat()) %>% 
  hc_navigator(enabled = FALSE) %>% 
  hc_scrollbar(enabled = FALSE) %>% 
  hc_exporting(enabled = TRUE) %>% 
  hc_legend(enabled = TRUE) 

# GRAPH EIGENVECTORS CORRELATION MATRIX
colnames(W) <- PCs_labs # prcomp only for labs
W <- as.data.frame(W)
W3 <- W[,c(1:3)]
W3 <- stack(list(PC1 = W3[,1], PC2 = W3[,2], PC3 = W3[,3]))
W3[,"Maturity"] <- Mat_labs 
colnames(W3) <- c("LoadingValue", "Eigenvectors", "Maturity")

ggplot(W3, aes(x = Maturity,
               y = LoadingValue, 
               colour = Eigenvectors, 
               group = Eigenvectors)) + 
  geom_line(size=1.5) +
  scale_x_discrete(limits = colnames(X)) +
  geom_point(size=2.5) +
  theme_light() +
  ylim(c(-0.5,0.75)) +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position="top",
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1"),
        strip.background = element_rect("#EDF0F1"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))
ggsave("/Users/Andrea/colour_correct_EIGENVECTORS.png",device = "png", dpi = 320, height=210, width=297, units="mm")
#ggsave("/Users/Andrea/Desktop/thesis/images/colour_correct.png", 
       #device = "png", dpi = 320, height=210, width=297, units="mm")
        

# PCA REPRESENTATION CORRELATION
Z3 <- Z[,1:3]
# Xhat contains in each column the vectors 
# of the Perc. int. rate changes at distinct mat.
R <- Z3%*%t(W[,1:3])
colnames(R) <- Mat_labs
Xhat <- xts(R,dates)
head(R)
R <- xts(R,dates)
highchart(type = "stock") %>% 
  hc_add_series(R[, Mat_labs[1]], name = Mat_labs[1]) %>% 
  hc_add_series(R[, Mat_labs[2]], name = Mat_labs[2]) %>%   
  hc_add_series(R[, Mat_labs[3]], name = Mat_labs[3]) %>%
  hc_add_series(R[, Mat_labs[4]], name = Mat_labs[4]) %>%
  hc_add_series(R[, Mat_labs[5]], name = Mat_labs[5]) %>%
  hc_add_series(R[, Mat_labs[6]], name = Mat_labs[6]) %>%
  hc_add_series(R[, Mat_labs[8]], name = Mat_labs[8]) %>%
  hc_add_series(R[, Mat_labs[9]], name = Mat_labs[9]) %>%
  hc_add_series(R[, Mat_labs[10]], name = Mat_labs[10]) %>%
  hc_add_series(R[, Mat_labs[11]], name = Mat_labs[11]) %>% 
  hc_add_theme(hc_theme_flat()) %>% 
  hc_navigator(enabled = FALSE) %>% 
  hc_scrollbar(enabled = FALSE) %>% 
  hc_exporting(enabled = TRUE) %>% 
  hc_legend(enabled = TRUE) 

# ABSOLUTE DIFFERENCE BTW cor(X) and cor(R)
abs(cor(X)-cor(Z))

# PCA ON SIX CONSECUTIVE TIME WINDOWS
len_window <- dim(X)[1]/6 # 582
n_windows = dim(X)[1]/len_window
tot_days = dim(X)[1]

# Create Indexes to slice X into 6 equal windows 
points <- c()
for(window in 1:n_windows){
  points <- c(points, window*len_window)
}

# Populate a list with the 6 Time series
Xs <- vector(mode = "list", length = n_windows)
for(p in points){
  Xs[[p/len_window]] <- X[(p-len_window):p]
}

# Prcomp on each window
# Prcomp is equivalent to spectral decomposition but use SVD
# For each window we have:
# Absolute changes in basis points + EigenVals + EigenVecs + PCs
windows_pca <- lapply(Xs, FUN = prcomp, scale = FALSE, center = FALSE)

# Select eigenvectors from each ts of windows_pca 
selected_eig_vecs = vector(mode = "list", length = n_windows)
for(i in 1:length(selected_eig_vecs)){
  selected_eig_vecs[[i]] <- windows_pca[[i]][[2]][,c("PC1","PC2","PC3")]
}

# Select eigenvalues from each ts of windows_pca 
selected_eig_vals = vector(mode = "list", length = n_windows)
for(i in 1:length(selected_eig_vals)){
  selected_eig_vals[[i]] <- windows_pca[[i]][[1]]
}
sum_selected_eigneval <- lapply(selected_eig_vals,sum)

var_exp <- vector(mode = "list", length = n_windows)
for(i in 1:length(selected_eig_vecs)){
  var_exp[[i]] <- (selected_eig_vals[[i]]/sum_selected_eigneval[[i]])*100
}

var_exp_first_3 <- vector(mode="list",length = n_windows)
for(i in 1:length(var_exp_first_3)){
  var_exp_first_3[[i]] <- sum(var_exp[[i]][1:3])
}

# EIGENVECTORS GRAPHS ON EACH 
# We investigate the persistence of the fundamental structure:
# parallel shift, tilt, curvature
# Transform Eigenvectors into data.frame for ggplot2 
selected_eig_vecs_df <- lapply(selected_eig_vecs, 
                               FUN = as.data.frame)

# Create Labels 
mat <- row.names(selected_eig_vecs_df[[1]])
mat <- str_replace(mat,"MAT", "")
maturities <- lapply(selected_eig_vecs_df, row.names)
PCs <- lapply(selected_eig_vecs_df,colnames)

# Stack the PCs of the 6 windows for ggplot2 faceting 
stacked_dfs <- lapply(selected_eig_vecs_df, stack)
for(i in 1:length(selected_eig_vecs_df)){
  stacked_dfs[[i]] <- cbind(stacked_dfs[[i]],mat)
}

# Recover start-date and end-date for each window
startdate <- lapply(Xs,first)
endate <- lapply(Xs,last)

startdate <- lapply(startdate, index)
endate <- lapply(endate,index)

plots <- vector(mode="list", length = length(stacked_dfs))
for(i in 1:length(stacked_dfs)){
  plots[[i]] <- ggplot(stacked_dfs[[i]], aes(x = mat, 
                                             y = values, 
                                             group = ind, 
                                             colour = ind, 
                                             fill=ind)) + 
    geom_line(size=1) +
    theme_light() +
    scale_x_discrete(limits = mat) +
    labs(x = "Maturity", y = "Loadings", colour = "Eigenvectors:") +
    scale_y_continuous(limits=c(-1.5,1.5), breaks = seq(-5,5,0.25)) +
    theme(legend.position = c(0.5,0.85), 
          legend.title=element_text(size = 16),
          legend.text=element_text(size = 14),
          axis.text.x=element_text(size = rel(1.5)),
          axis.text.y=element_text(size=rel(1.5)),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          panel.grid.minor = element_blank(),
          plot.background = element_rect("#EDF0F1"),
          panel.background = element_rect(fill = "#EDF0F1"),
          legend.background = element_rect("#EDF0F1"),
          strip.background = element_rect("#EDF0F1"),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(str_c("Eigenvectors of window:","from",startdate[[i]],"to",endate[[i]],"\nWindow length: 582 Trading Days" , sep = " ")) +
    ggsave(str_c("PCA_window_colour_correct",as.character(i), ".png", sep=""), dpi = 320, width = 210, height = 210, units="mm", device = "png")
}


# BOOTSTRAP 
# Start from data TS: Absolute changes in basis points
# We create 10.000 bootstrap samples each made of 582 obs.
# Ratio 1/6 = 582/3492
# Each sample is stored into the list 'boot_samples'
reps=10000
boot_samples <- vector(mode="list", length = reps)
set.seed(1)
for(rep in 1:reps){
  boot_samples[[rep]] <- X[sample(1:dim(X)[1], 582, replace=TRUE),]
}

covariances <- lapply(boot_samples,cov)
E <- lapply(covariances, eigen, symmetric = TRUE)

# BOOTSTRAP EIGENVALUES
# Isolate the first eigenvalue for each bootstrap sample
# In total 10.000 lambda1
eigenvalues_first = vector(mode = "list", length = reps)
for(i in 1:length(E)){
  eigenvalues_first[[i]] <- E[[i]][[1]][1]     
}
eig_first <- unlist(eigenvalues_first)
summary(eig_first)
round(quantile(eig_first,c(0.025,0.975)),2)
round(sd(eig_first),2)

# Isolate the second eigenvalue for each bootstrap sample
# In total 10.000 lambda2
eigenvalues_second = vector(mode = "list", length = reps)
for(i in 1:length(E)){
  eigenvalues_second[[i]] <- E[[i]][[1]][2]      
}
eig_second <- unlist(eigenvalues_second)
summary(eig_second)
round(quantile(eig_second,c(0.025,0.975)),2)
round(sd(eig_second),2)

# Isolate the third eigenvalue for each bootstrap sample
# In total 10.000 lambda3
eigenvalues_third = vector(mode = "list", length = reps)
for(i in 1:length(E)){
  eigenvalues_third[[i]] <- E[[i]][[1]][3]      
}
eig_third <- unlist(eigenvalues_third)
summary(eig_third)
round(quantile(eig_third,c(0.025,0.975)),2)
round(sd(eig_third),2)

# We organize the 10.000*lambda1, 10.000*lambda2, 10.000*lambda3
# In order to have a stacked data.frame for ggplot2 faceting

# Create labels "first", "second", "third"
label_first <- rep("first",reps)
label_second <- rep("second",reps)
label_third <- rep("third",reps)

labels <- c(label_first,label_second,label_third)
eigenvalues <- c(eig_first,eig_second,eig_third)

eigenvalues_df <- data.frame(Lambda=as.numeric(eigenvalues), 
                             Label=as.character(labels))

# Densities
ggplot(eigenvalues_df, aes(x=Lambda, 
                           y=..density.., 
                           fill = Label)) +
  geom_histogram(bins = 250) +
  geom_density(size=0.5,colour = "grey60") +
  facet_grid(Label ~ .) +
  theme_light() +
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="top",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1"),
        strip.background = element_rect("#EDF0F1")) +
  ggtitle("Eigenvalues Bootstrap Estimates") 
#ggsave("/Users/Andrea/Desktop/thesis/images/colour_correct2.png", device = "png", dpi = 320, height=210, width=210, units="mm")

# Box plots
ggplot(eigenvalues_df, aes(x=Label, y=Lambda)) + 
  geom_boxplot(outlier.size=1.5,outlier.shape=21) + 
  theme_light() +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(size = rel(2)),
        axis.text.y=element_text(size=rel(1.25)),
        plot.title = element_text(hjust=0.5)) +
  scale_x_discrete(labels=c(expression(lambda[1]),
                            expression(lambda[2]),
                            expression(lambda[3]))) +
  ggtitle("Eigenvalues Bootstrap Estimates, Num of reps. = 10.000, Size = 582")

#### Bootstrapped Variance Explained by the first three eigenvalues ####
eigenvalues_all = vector(mode = "list", length = reps)
for(i in 1:length(E)){
  eigenvalues_all[[i]] <- E[[i]][[1]]      
}

eigenvalue_sum <- lapply(eigenvalues_all, sum)
eig_sum_all <- unlist(eigenvalue_sum)

# Bootstrapped First three Eigenvalues in a 3 column dataframe
eigenvalues_df_2 <- as.data.frame(cbind(eig_first, 
                                        eig_second, 
                                        eig_third), 
                                  stringsAsFactors = FALSE)
dim(eigenvalues_df_2)
colnames(eigenvalues_df_2) <- c("Lambda1","Lambda2","Lambda3")

# Sum of the first three eigenvalues on each bootstrap sample
sum_first_eig <- apply(eigenvalues_df_2, 1, sum)

# Ratio between first three eigenvals vs all eigenvals
explained <- data.frame(Var_Explained = sum_first_eig/eig_sum_all)
head(explained)
summary(explained)
round(quantile(explained[,"Var_Explained"],c(0.025,0.975)),3)
round(sd(explained[,"Var_Explained"]),3)

# Kernel estimated density of the variance explained 
# by the first three eigenvals 
ggplot(explained, aes(x=Var_Explained, y=..density..)) + 
  geom_histogram(bins = 50, colour = "grey60", fill = cols[2]) +
  geom_density() +
  theme_light() +
  xlab(expression((lambda[1]+lambda[2]+lambda[3]) / 
                    (lambda[1]+lambda[2]+...+lambda[11]))) +
  theme(axis.text.x = element_text(size = rel(1.25)), 
        axis.title.x = element_text(size = rel(1.25)),
        plot.title = element_text(hjust = 0.5),
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1"),
        strip.background = element_rect("#EDF0F1")) +
  ggtitle("Variance Explained by the first three PCs on 10.000 Bootstrap samples") +
ggsave("/Users/Andrea/Desktop/thesis/images/colour_correct3.png",device = "png", dpi = 320, height=210, width=210, units="mm")

# Bootstrap loadings of the first 3 eigenvectors
# Create a list that contains the first three eigenvectors for each bootstrap sample
# We take them from the original list that contained everything, named E
# Index traial
E[[1]][[2]][,c(1,2,3)]
loadings_eig = vector(mode = "list", length = reps)
length(loadings_eig)
for(i in 1:length(E)){
  loadings_eig[[i]] <- E[[i]][[2]][,c(1,2,3)]
}

# Now we create PC1 which is a matrix with:
# 11 columns corresponding to the numb. of loadings of the FIRST eigenvector
# 10.000 rows corresponding to the estimated FIRST eigenvector on the i-th 
# bootstrap sample, with i = 1,...,10.000
# In few, words we organize each bootstrap first eigenvector as it was transposed
w1 <- matrix(ncol = 11, nrow = reps)
# Same as above but for the SECOND eigenvector
w2 <- matrix(ncol = 11, nrow = reps)
# Same as above but for the THIRD eigenvector
w3 <- matrix(ncol = 11, nrow = reps)

# Populate PC1, PC2, PC3
for(i in 1:reps){
  w1[i,] <- unlist(loadings_eig[[i]][,1])
  w2[i,] <- unlist(loadings_eig[[i]][,2])
  w3[i,] <- unlist(loadings_eig[[i]][,3])
}

# Now we stack PC1, PC2, PC3 together by row for ggplot2
# The first columns contains the first loading on PC1, PC2, PC3,
# The second column contains the second loading on PC1, PC2, PC3
# and so on.
PCs_boot <- as.data.frame(rbind(w1,w2,w3))
dim(PCs_boot)
head(PCs_boot)

# We rename the 11 columns of PCs_boot
# in order to identify the loadings of the first three eigenvectors
colnames(PCs_boot) <- Mat_labs
# loadings_labels <- c()
# for(i in 1:11){
#   loadings_labels[i] <- str_c("a",as.character(i))
# }
# loadings_labels
# colnames(PCs_boot) <- loadings_labels

# Now delicate part.
# Stack PCs_boot
# The rows of PCs_boots are stacked in column one over the other
A <- as.data.frame(stack(PCs_boot))

# Associate correctly the labels PC1, PC2, PC3 for faceting
# In PCs_labs we 330.000 rows which corresponds to 11 
# groups of rows with each PC1*10k PC2*10k PC3*10k
PCs_labs_330.000 <- vector(mode="list", length = 11)
for(i in 1:11){
  PCs_labs_330.000[[i]] <- c(rep("w1",10000),
                             rep("w2",10000),
                             rep("w3",10000))
}
labsss <- unlist(PCs_labs_330.000)

# Associate the labs to each each loading of each eigenvector
A <- cbind(A,labsss)

# Facet graph of bootstrapped loadings of first 3 eigenvectors
# Estimated using Centered Percentage Changes 
# (It takes some time)
ggplot(A,aes(x=values, fill = labsss)) + 
  geom_histogram(bins=500) + 
  facet_grid(ind~labsss) +
  theme_light() +
  guides(fill=FALSE) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eeeeee",
                                        colour="grey"),
        strip.text = element_text(colour="black"),
        plot.title = element_text(hjust=0.5),
        plot.background = element_rect("#EDF0F1"),
        panel.background = element_rect(fill = "#EDF0F1"),
        legend.background = element_rect("#EDF0F1")) +
  labs(x="Value of the Loadings", y="Frequency")
ggsave("/Users/Andrea/Desktop/thesis/images/colour_correct4.png",device = "png", dpi = 320, height=297, width=210, units="mm")
