library(dplyr)
library(CEC)
library(reshape2)
library(TraMineR)
library(tidyr)
library(janitor)
library(caret)
library(car)
library(rgl)
library(mice)
library(car)
library(philentropy)
library(anocva)
library(dbscan)
library(MASS)
#library(SNFtool)


save.image("v_progress.RData")
load("v2.RData")

set.seed(499)

#Reads data into memory
data_LOW <- read.table(file = 'anonymized_ch_low.tsv', sep = ",", header = TRUE)
data_MID <- read.table(file = 'anonymized_ch_mid.tsv', sep = ",", header = TRUE)
data_LOW <- remove_empty(data_LOW, which = c("rows","cols"))
data_MID <- remove_empty(data_MID, which = c("rows","cols"))
customeractions <- t(read.csv("fromcustomer.csv", header = TRUE))
bankactions     <- t(read.csv("frombank.csv", header = TRUE))


#Calcualte PDFs of events for each customer -MID
t1       <- dplyr::select(data_LOW, "NEW_CLIENT_NO", "EVENT_DEFINITION_CODE", "EVENT_START_DATE") 
mlt1     <- melt(t1,id.vars =  names(t1)[c(1,3)])
freq_1   <- dcast(mlt1, NEW_CLIENT_NO  ~ value, value.var = c("NEW_CLIENT_NO"),fun.aggregate =length)
freq_low <- data.frame(freq_1[,-1], row.names = freq_1[,1])
dist_low <- diag(1/rowSums(freq_low)) %*% data.matrix(freq_low)

#Calcualte PDFs of events for each customer -LOW
t2       <- dplyr::select(data_MID, "NEW_CLIENT_NO", "EVENT_DEFINITION_CODE", "EVENT_START_DATE") 
mlt2     <- melt(t2,id.vars =  names(t2)[c(1,3)])
freq_2   <- dcast(mlt2, NEW_CLIENT_NO  ~ value, value.var = c("NEW_CLIENT_NO"),fun.aggregate =length)
freq_mid <- data.frame(freq_2[,-1], row.names = freq_2[,1])
dist_mid <- diag(1/rowSums(freq_mid)) %*% data.matrix(freq_mid)

#Calcualte seperate PDFs for bank and customer actions for each customer - MID
intersectbank      <- intersect(bankactions, names(freq_mid))
intersectcustomer  <- intersect(customeractions, names(freq_mid))
intersectbank2      <- intersect(bankactions, names(freq_low))
intersectcustomer2  <- intersect(customeractions, names(freq_low))

freq_mid_bank <- dplyr::select(freq_mid, !!intersectbank)
freq_mid_customer <- dplyr::select(freq_mid, !!intersectcustomer)

freq_low_bank <- dplyr::select(freq_low, !!intersectbank2)
freq_low_customer <- dplyr::select(freq_low, !!intersectcustomer2)

dist_mid_bank     <- diag(1/rowSums(freq_mid_bank)) %*% data.matrix(freq_mid_bank)
dist_mid_customer <- diag(1/rowSums(freq_mid_customer)) %*% data.matrix(freq_mid_customer)
dist_mid_customer[is.na(dist_mid_customer)] <- 0

dist_low_bank     <- diag(1/rowSums(freq_low_bank)) %*% data.matrix(freq_low_bank)
dist_low_customer <- diag(1/rowSums(freq_low_customer)) %*% data.matrix(freq_low_customer)
dist_low_customer[is.na(dist_low_customer)] <- 0
dist_low_bank[is.na(dist_low_bank)] <- 0

#MID + LOW --- didn'work
#freqnew <- freqnew[, !duplicated(t(f1reqnew))]
#bind seperate datasets
# new_names <- intersect(names(freq_low), names(freq_mid))
# freq <- bind_rows(freq_low[, new_names], freq_mid[, names(freq_mid)])
# dist <- diag(1/rowSums(freq)) %*% data.matrix(freq)




#MIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMID
#MIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMID
#MIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMID
#MIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMIDMID




#Calculate various distance metrics - MID

CHmid_abnormal <- as.dist(distance(freq_mid, method = "squared_chi"))

CHmid <- as.dist(distance(dist_mid, method = "squared_chi"))
JSmid <- as.dist(distance(dist_mid, method = "jensen-shannon"))
KJmid <- as.dist(distance(dist_mid, method = "kumar-johnson"))
JCmid <- as.dist(distance(dist_mid, method = "jaccard"))
SCmid <- as.dist(distance(dist_mid, method = "squared_chord"))
TJmid <- as.dist(distance(dist_mid, method =  "taneja"))
CBmid <- as.dist(distance(dist_mid, method = "canberra"))
LRmid <- as.dist(distance(dist_mid, method = "lorentzian"))


#distance metrics MID-bank only / customer only
CHmid_bank <- as.dist(distance(dist_mid_bank, method = "squared_chi"))
JSmid_bank <- as.dist(distance(dist_mid_bank, method = "jensen-shannon"))
JCmid_bank <- as.dist(distance(dist_mid_bank, method = "jaccard"))


CHmid_customer <- as.dist(distance(dist_mid_customer, method = "squared_chi"))
JSmid_customer <- as.dist(distance(dist_mid_customer, method = "jensen-shannon"))
JCmid_customer <- as.dist(distance(dist_mid_customer, method = "jaccard"))

CHmid_customer[CHmid_customer<=0] <- 10^(-8) 
JSmid_customer[JSmid_customer<=0] <- 10^(-8)
JCmid_customer[JCmid_customer<=0] <- 10^(-8)



#Perform Non-metric multidimensional scaling for visualisation -MID
NMDS_CH_mid_abnormal <-  isoMDS(CHmid_abnormal, k=3, maxit=20)

NMDS_CH_mid <- isoMDS(CHmid, k=3, maxit=20)
NMDS_JS_mid <- isoMDS(JSmid, k=3, maxit=20)
NMDS_KJ_mid <- isoMDS(KJmid, k=3, maxit=20)
NMDS_JC_mid <- isoMDS(JCmid, k=3, maxit=20)
NMDS_SC_mid <- isoMDS(SCmid, k=3, maxit=20)
NMDS_TJ_mid <- isoMDS(TJmid, k=3, maxit=20)
NMDS_CB_mid <- isoMDS(CBmid, k=3, maxit=20)
NMDS_LR_mid <- isoMDS(LRmid, k=3, maxit=20)

#NM MDS - Bank only / customer only
NMDS_CH_mid_bank <- isoMDS(CHmid_bank, k=3, maxit=20)
NMDS_JS_mid_bank <- isoMDS(JSmid_bank, k=3, maxit=20)
NMDS_JC_mid_bank <- isoMDS(JCmid_bank, k=3, maxit=20)

NMDS_CH_mid_customer <- isoMDS(CHmid_customer, k=3, maxit=20)
NMDS_JS_mid_customer <- isoMDS(JSmid_customer, k=3, maxit=20)
NMDS_JC_mid_customer <- isoMDS(JCmid_customer, k=3, maxit=20)



#HDBSCAN - spectral + hierarchial best of both worlds -MID
mid_CH_hdbscan_bank <- hdbscan(CHmid_bank, minPts = 20)
mid_JS_hdbscan_bank <- hdbscan(JSmid_bank, minPts = 10)
mid_JS_hdbscan <- hdbscan(JSmid, minPts = 10)

mid_CH_hdbscan_customer <- hdbscan(CHmid_customer, minPts = 20)
mid_JS_hdbscan_customer <- hdbscan(JSmid_customer, minPts = 10)
mid_JC_hdbscan_customer <- hdbscan(JCmid_customer, minPts = 20)

mid_CH_hdbscan_customer <- hdbscan(CHmid_customer, minPts = 20)
mid_JS_hdbscan_customer <- hdbscan(JSmid_customer, minPts = 20)



#plot - MID customer
c<-NMDS_JS_mid[["points"]]
car::scatter3d(c[,1], c[,2], c[,3],xlab="x", ylab="y", zlab="z",surface=FALSE, point.col="blue", sphere.size=1.5, axis.ticks=TRUE)


#plot - MID bank
b<-NMDS_JS_mid[["points"]]
car::scatter3d(b[,1], b[,2], b[,3],xlab="x", ylab="y", zlab="z",surface=FALSE, point.col="blue", sphere.size=1.5, axis.ticks=TRUE)

#plot - MID
a<-NMDS_JS_mid_customer[["points"]]
car::scatter3d(a[,1], a[,2], a[,3],xlab="y", ylab="z", zlab="x",
               surface=FALSE, point.col="blue", sphere.size=1.5, 
               axis.ticks=TRUE,
               groups=factor(mid_JS_hdbscan_customer[["cluster"]]))

play3d( spin3d( axis = c(0, 1, 0), rpm = 3), duration = 15 )

# Save like gif
movie3d(
        movie="JSmid", 
        spin3d( axis = c(0, 1, 0), rpm = 7),
        duration = 10, 
        type = "gif", 
        clean = TRUE
)

        

#Hierarchial Aglomative Clustering - MID
mid_CH_hclust_customer <- hclust(CHmid_customer)
plot(mid_CH_hclust, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = NULL, ylab = "Height")
cut_midCH_customer <- cutree(mid_CH_hclust_customer, k=3)
print(c(sum(cut_midCH==1),sum(cut_midCH==2),sum(cut_midCH==3)))

#Spectral Clustering - MID
mid_CH_spectral <- spectralClustering(as.matrix(CHmid), 3, type=3)
print(c(sum(mid_CH_spectral==1),sum(mid_CH_spectral==2),sum(mid_CH_spectral==3)))

table(mid_JS_spectral, mid_LR_spectral)

#Hierarchial clustering using Wasserstein Distance -------DW
#WH_hclust(as.data.frame(freq_mid))




#LOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOW
#LOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOW
#LOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOWLOW




#CH     <- as.dist(distance(dist,     method = "squared_chi"))
CHlow  <- as.dist(distance(dist_low, method = "squared_chi"))
JSlow  <- as.dist(distance(dist_low, method = "jensen-shannon"))
JClow  <- as.dist(distance(dist_low, method = "jaccard"))
SClow  <- as.dist(distance(dist_low, method = "squared_chord"))
TJlow  <- as.dist(distance(dist_low, method =  "taneja"))
#CBlow  <- as.dist(distance(dist_low, method = "canberra"))
#LRlow <- as.dist(distance(dist_low, method = "lorentzian"))
#KJlow <- as.dist(distance(dist_low, method = "kumar-johnson"))

#distance metrics LOW-bank only / customer only
CHlow_bank <- as.dist(distance(dist_low_bank, method = "squared_chi"))
JSlow_bank <- as.dist(distance(dist_low_bank, method = "jensen-shannon"))

CHlow_customer <- as.dist(distance(dist_low_customer, method = "squared_chi"))
JSlow_customer <- as.dist(distance(dist_low_customer, method = "jensen-shannon"))

CHlow_customer[CHlow_customer<=0] <- 10^(-8) 
JSlow_customer[JSlow_customer<=0] <- 10^(-8)

CHlow_bank[CHlow_bank<=0] <- 10^(-8) 
JSlow_bank[JSlow_bank<=0] <- 10^(-8)

CHlow[CHlow<=0] <- 10^(-8)
JSlow[JSlow<=0] <- 10^(-8)
JClow[JClow<=0] <- 10^(-8)
SClow[SClow<=0] <- 10^(-8)



if (any(CHlow==0)) {print("bad")} else {print("ok")}

sum(JSmid==0)

NMDS_CH_low <-  isoMDS(CHlow, k=3, maxit=10)
NMDS_JS_low <-  isoMDS(JSlow, k=3, maxit=10)
NMDS_JC_low <-  isoMDS(JClow, k=3, maxit=10)
NMDS_SC_low <-  isoMDS(SClow, k=3, maxit=10)

NMDS_CH_low_bank <-  isoMDS(CHlow_bank, k=3, maxit=10)
NMDS_JS_low_bank <-  isoMDS(JSlow_bank, k=3, maxit=10)
NMDS_CH_low_customer <-  isoMDS(CHlow_customer, k=3, maxit=10)
NMDS_JS_low_customer <-  isoMDS(JSlow_customer, k=3, maxit=10)

d<-NMDS_JS_customer[["points"]]
palette(rainbow(19))
car::scatter3d(d[,1], d[,2], d[,3],xlab="y", ylab="z", zlab="x",
               surface=FALSE, surface.col=1:19, sphere.size=3, 
               axis.ticks=TRUE,
               groups=factor(low_JS_hdbscan_customer[["cluster"]]))


low_JS_hdbscan <- hdbscan(JSlow,minPts = 30)
low_JS_hdbscan_bank <- hdbscan(JSlow_bank, minPts = 30)
low_JS_hdbscan_customer <- hdbscan(JSlow_customer, minPts=30)

#CH_clusters <- hdbscan(CH, minPts = 50)

#Hierarchial Aglomative Clustering
CH_hclust_low <- hclust(CHlow)
plot(CH_hclust_low, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram",
     sub = NULL, xlab = NULL, ylab = "Height")
cut_CH_low <- cutree(CH_hclust_low, k=3)
print(c(sum(cut_CH_low==1),sum(cut_CH_low==2),sum(cut_CH_low==3)))

# Spectral Clustering is too slow for this application - for low
#PCA meaningless - 100 dimensions for 90% POEV
#freq_1_pca <- PCA(freq_1[, 2:127], scale. = TRUE)