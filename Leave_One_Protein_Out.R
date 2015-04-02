#!/usr/bin/env Rscript
C <- read.csv("data/Compound.csv", header=TRUE)
P <- read.csv("data/Protein.csv", header=TRUE)
Y <- read.csv("data/Activity.csv", header=TRUE)
index <- read.csv("data/Protein_Index.csv", header=TRUE)

set.seed(400)
compound_corr <- cor(C)
dim(compound_corr)
compound_corr <- compound_corr[1:44, 1:44]
compound_highCorr <- findCorrelation(compound_corr, cutoff = .70)
length(compound_highCorr)
highCorrRemoveCompound <- C[, -compound_highCorr]
Activity <- Y$Activity
Compound <- cbind(Activity, highCorrRemoveCompound)
set.seed(500)
P = P[, -nearZeroVar(P)]
protein_corr <- cor(P)
dim(protein_corr)
protein_corr <- protein_corr[1:390, 1:390]
protein_highCorr <- findCorrelation(protein_corr, cutoff = .70)
length(protein_highCorr)
highCorrRemoveProtein <- P[, -protein_highCorr]
Protein <- cbind(Activity, highCorrRemoveProtein)
## Leave-One-Out-Protein for compound alone
set.seed(433330)
dat <- cbind(index, Compound)
subs <- unique(dat$Protein)
model_these <- vector(mode = "list", length= length(subs))
for(i in seq_along(subs))
  model_these[[i]] <- which(dat$Protein != subs[i])
names(model_these) <- paste0("Protein", subs)
LOPO_compound_alone <- train(x= dat[, 3:26],
                             y= dat[, 2],
                             method= "J48",
                             trControl=trainControl(method="cv",
                                                    index = model_these,
                                                    classProbs= TRUE))
pred <- predict(LOPO_compound_alone, dat)
actual_LOPO_compound_alone <- dat$Activity
matrix_LOPO_compound_alone <- table(pred, actual_LOPO_compound_alone)
###LOPO for Protein Alone
set.seed(44442)
dat <- cbind(index, Protein)
subs <- unique(dat$Protein)
model_these <- vector(mode = "list", length= length(subs))
for(i in seq_along(subs))
  model_these[[i]] <- which(dat$Protein != subs[i])
names(model_these) <- paste0("Protein", subs)
LOPO_protein_alone <- train(x= dat[, 3:128],
                             y= dat[, 2],
                             method= "J48",
                             trControl=trainControl(method="cv",
                                                    index = model_these,
                                                    classProbs= TRUE))
pred <- predict(LOPO_protein_alone, dat)
actual_LOPO_protein_alone <- dat$Activity
matrix_LOPO_protein_alone <- table(pred, actual_LOPO_protein_alone)
## LOPO for crossterms
set.seed(400023)
c <- Compound
p <- Protein
train_c <- c[-1]
train_p <- p[-1]
crossterms <- getCPI(train_c, train_p, type = "tensorprod")
crossterms <- data.frame(crossterms)
dat <- cbind(index, Activity, crossterms)
subs <- unique(dat$Protein)
model_these <- vector(mode = "list", length= length(subs))
for(i in seq_along(subs))
  model_these[[i]] <- which(dat$Protein != subs[i])
names(model_these) <- paste0("Protein", subs)
LOPO_crossterms <- train(x= dat[, 3:3026],
                             y= dat[, 2],
                             method= "J48",
                             trControl=trainControl(method="cv",
                                                    index = model_these,
                                                    classProbs= TRUE))
pred <- predict(LOPO_crossterms, dat)
actual_LOPO_crossterms <- dat$Activity
matrix_LOPO_crossterms <- table(pred, actual_LOPO_crossterms)
set.seed(43333)
####LOPO for protein_protein Crossterms
PxP_crossterms <- getCPI(train_p, train_p, type = "tensorprod")
TranposedIndexed_protein <- t(PxP_crossterms)
index1 <- which(duplicated(TranposedIndexed_protein))
removed_protein_train <- TranposedIndexed_protein[-index1, ]
protein_final_PxP_crossterms <- t(removed_protein_train)
protein_final_PxP_crossterms <- data.frame(protein_final_PxP_crossterms)
dat <- cbind(index, Activity, protein_final_PxP_crossterms)
subs <- unique(dat$Protein)
model_these <- vector(mode = "list", length= length(subs))
for(i in seq_along(subs))
  model_these[[i]] <- which(dat$Protein != subs[i])
names(model_these) <- paste0("Protein", subs)
LOPO_PxP_crossterms <- train(x= dat[, 3:8003],
                         y= dat[, 2],
                         method= "J48",
                         trControl=trainControl(method="cv",
                                                index = model_these,
                                                classProbs= TRUE))
pred <- predict(LOPO_PxP_crossterms, dat)
actual_LOPO_PxP_crossterms <- dat$Activity
matrix_LOPO_PxP_crossterms <- table(pred, actual_LOPO_PxP_crossterms)
set.seed(320)
### LOPO for compound-compound crossterms
CxC_crossterms <- getCPI(train_c, train_c, type = "tensorprod")
TranposedIndexed_compound <- t(CxC_crossterms)
index2 <- which(duplicated(TranposedIndexed_compound))
removed_compound_train <- TranposedIndexed_compound[-index2, ]
protein_final_CxC_crossterms <- t(removed_compound_train)
protein_final_CxC_crossterms <- data.frame(protein_final_CxC_crossterms)
dat <- cbind(index, Activity, protein_final_CxC_crossterms)
subs <- unique(dat$Protein)
model_these <- vector(mode = "list", length= length(subs))
for(i in seq_along(subs))
  model_these[[i]] <- which(dat$Protein != subs[i])
names(model_these) <- paste0("Protein", subs)
LOPO_CxC_crossterms <- train(x= dat[, 3:129],
                             y= dat[, 2],
                             method= "J48",
                             trControl=trainControl(method="cv",
                                                    index = model_these,
                                                    classProbs= TRUE))
pred <- predict(LOPO_CxC_crossterms, dat)
actual_LOPO_CxC_crossterms <- dat$Activity
matrix_LOPO_CxC_crossterms <- table(pred, actual_LOPO_CxC_crossterms)
## Statistical Assessment
set.seed(9999)
z <- data.frame(matrix_LOPO_compound_alone)
y <- data.frame(matrix_LOPO_protein_alone)
x <- data.frame(matrix_LOPO_crossterms)
w <- data.frame(matrix_LOPO_PxP_crossterms)
v <- data.frame(matrix_LOPO_CxC_crossterms)
results <- cbind(z$Freq, y$Freq, x$Freq, w$Freq, v$Freq)
data <- data.frame(results)
m = ncol(data)
ACC  <- matrix(nrow = m, ncol = 1)
SENS  <- matrix(nrow = m, ncol = 1)
SPEC  <-matrix(nrow = m, ncol = 1)
MCC <- matrix(nrow = m, ncol = 1)

for(i in 1:m){ 
  ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
  SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
  SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
  MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
  MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
  MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
  MCC4  =  sqrt(MCC2)*sqrt(MCC3)
  
  
  MCC[i,1]  = MCC1/MCC4
}

Model <- c("Compound",
           "Protein",
           "Compound-Protein",
           "Protein-Protein",
           "Compound-Compound")
Perf = data.frame(ACC,SENS,SPEC,MCC)
is.num <- sapply(Perf, is.numeric)
Perf[is.num] <- lapply(Perf[is.num], round, 2)
table <- data.frame(Model, Perf)
pdf("Performance_of_Leave_One_Protein_Out.pdf")
grid.table(table,gpar.coretext=gpar(fontsize=12),
           show.box = FALSE,
           show.csep = FALSE,
           core.just="center",
           equal.width=FALSE,
           show.rownames=FALSE,
           row.just = "right",
           show.vlines = FALSE,
           gpar.corefill = gpar(fill = "white", col = "white"),
           gpar.colfill = gpar(fill = "white", col = "white"))
dev.off()


