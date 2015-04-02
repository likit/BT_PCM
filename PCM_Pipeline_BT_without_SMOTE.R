#!/usr/bin/env Rscript
C <- read.csv("data/Compound.csv", header=TRUE)
P <- read.csv("data/Protein.csv", header=TRUE)
Y <- read.csv("data/Activity.csv", header=TRUE)
set.seed(400)
compound_corr <- cor(C)
dim(compound_corr)
compound_corr <- compound_corr[1:44, 1:44]
compound_highCorr <- findCorrelation(compound_corr, cutoff = .70)
length(compound_highCorr)
highCorrRemoveCompound <- C[, -compound_highCorr]
Activity <- Y$Activity
InputKennardCompound <- cbind(Activity, highCorrRemoveCompound)
set.seed(500)
P = P[, -nearZeroVar(P)]
protein_corr <- cor(P)
dim(protein_corr)
protein_corr <- protein_corr[1:390, 1:390]
protein_highCorr <- findCorrelation(protein_corr, cutoff = .70)
length(protein_highCorr)
highCorrRemoveProtein <- P[, -protein_highCorr]
InputKennardProtein <- cbind(Activity, highCorrRemoveProtein)
C_Inactive <- subset(InputKennardCompound, Activity == "Inactive")
C_Active <- subset(InputKennardCompound, Activity == "Active")
P_Inactive <- subset(InputKennardProtein, Activity == "Inactive")
P_Active <- subset(InputKennardProtein, Activity == "Active")
# Spliting the data set into training set and testing test using Kennard Stone Algorithm
Inactive <- C_Inactive
x <- data.frame(Inactive)
sel <- kenStone(x[-1], k = 195, metric = "mahal", pc=2)
trainInactiveCompound <- x[sel$model, ]
testInactiveCompound <- x[sel$test, ]
Active <- C_Active
x <- data.frame(Active)
sel <- kenStone(x[-1], k = 21, metric = "mahal", pc=2)
trainActiveCompound <- x[sel$model, ]
testActiveCompound <- x[sel$test, ]
Inactive <- P_Inactive
x <- data.frame(Inactive)
sel <- kenStone(x[-1], k = 195, metric = "mahal", pc=2)
trainInactiveProtein <- x[sel$model, ]
testInactiveProtein <- x[sel$test, ]
Active <- P_Active
x <- data.frame(Active)
sel <- kenStone(x[-1], k = 21, metric = "mahal", pc=2)
trainActiveProtein <- x[sel$model, ]
testActiveProtein <- x[sel$test, ]
# Preparing input data for Proteocheometrics modeling
set.seed(567)
C_Train <- rbind(trainInactiveCompound, trainActiveCompound)
C_Test <- rbind(testInactiveCompound, testActiveCompound)
P_Train <- rbind(trainInactiveProtein, trainActiveProtein)
P_Test <- rbind(testInactiveProtein, testActiveProtein)
c <- C_Train
p <- P_Train
train_c <- c[-1]
train_p <- p[-1]
crossTrain <- getCPI(train_c, train_p, type = "tensorprod")
c <- C_Test
p <- P_Test
test_c <- c[-1]
test_p <- p[-1]
crossTest <- getCPI(test_c, test_p, type = "tensorprod")
Activity_Train <- C_Train$Activity
Activity_Test <- C_Test$Activity
CxP_Train <- cbind(Activity_Train, data.frame(crossTrain))
CxP_Test <- cbind(Activity_Test, data.frame(crossTest))
protein_selfcrosstrain <- getCPI(train_p, train_p, type = "tensorprod")
TranposedIndexed_protein <- t(protein_selfcrosstrain)
index1 <- which(duplicated(TranposedIndexed_protein))
removed_protein_train <- TranposedIndexed_protein[-index1, ]
protein_finalselfcrosstrain <- t(removed_protein_train)
PxP_Train <- cbind(Activity_Train, data.frame(protein_finalselfcrosstrain))
protein_selfcrosstest <- getCPI(test_p, test_p, type = "tensorprod")
TransposedIndexed_protein <- t(protein_selfcrosstest)
index2 <- which(duplicated(TransposedIndexed_protein))
remove_protein_test <- TransposedIndexed_protein[-index2, ]
protein_finalselfcrosstest <- t(remove_protein_test)
PxP_Test <- cbind(Activity_Test, data.frame(protein_finalselfcrosstest))
compound_selfcrosstrain <- getCPI(train_c, train_c, type = "tensorprod")
TranposedIndexed_compound <- t(compound_selfcrosstrain)
index3 <- which(duplicated(TranposedIndexed_compound))
removed_compound_train <- TranposedIndexed_compound[-index3, ]
compound_finalselfcrosstrain <- t(removed_compound_train)
CxC_Train <- cbind(Activity_Train, data.frame(compound_finalselfcrosstrain))
compound_selfcrosstest <- getCPI(test_c, test_c, type = "tensorprod")
TransposedIndexed_compound <- t(compound_selfcrosstest)
index4 <- which(duplicated(TransposedIndexed_compound))
remove_compound_test <- TransposedIndexed_compound[-index4, ]
compound_finalselfcrosstest <- t(remove_compound_test)
CxC_Test <- cbind(Activity_Test, data.frame(compound_finalselfcrosstest))
##Bulding a Predictive Proteocheometrics Modelling with Caret Package
set.seed(4000)
## for Compound alone
fitControl <- trainControl(method= "repeatedcv", number = 10, repeats = 10)
compound_alone_model <- train(x = C_Train[,2:25],
                              y = C_Train[,1], method = "J48", trControl= fitControl)
pred <- predict(compound_alone_model, C_Train)
actual_10CV <- C_Train$Activity
matrix_C_alone_10CV <- table(pred, actual_10CV)
prediction <- predict(compound_alone_model, C_Test)
actual_test <- C_Test$Activity  
matrix_C_alone_test <- table(prediction, actual_test)
## for Protein alone
set.seed(5532)
protein_alone_model <- train(x = P_Train[,2:127],
                             y = P_Train[,1], method = "J48", trControl= fitControl)
pred <- predict(protein_alone_model, P_Train)
actual_10CV <- P_Train$Activity
matrix_P_alone_10CV <- table(pred, actual_10CV)
prediction <- predict(protein_alone_model, P_Test)
actual_test <- P_Test$Activity  
matrix_P_alone_test <- table(prediction, actual_test)
## for Compound_Portein Crossterms
set.seed(6667)
compound_protein_model <- train(x = CxP_Train[,2:3024],
                                y = CxP_Train[,1], method = "J48", trControl= fitControl)
pred <- predict(compound_protein_model, CxP_Train)
actual_10CV <- CxP_Train$Activity_Train
matrix_compound_protein_10CV <- table(pred, actual_10CV)
prediction <- predict(compound_protein_model, CxP_Test)
actual_test <- CxP_Test$Activity_Test 
matrix_compound_protein_test <- table(prediction, actual_test)
## for Protein_Protein Crossterms
set.seed(6667)
protein_protein_model <- train(x = PxP_Train[,2:7875],
                               y = PxP_Train[,1], method = "J48", trControl= fitControl)
pred <- predict(protein_protein_model, PxP_Train)
actual_10CV <- PxP_Train$Activity_Train
matrix_protein_protein_10CV <- table(pred, actual_10CV)
prediction <- predict(protein_protein_model, PxP_Test)
actual_test <- PxP_Test$Activity_Test 
matrix_protein_protein_test <- table(prediction, actual_test)
## For Compound_Compound Crossterms
set.seed(6667)
compound_compound_model <- train(x = CxC_Train[,2:276],
                                 y = CxC_Train[,1], method = "J48", trControl= fitControl)
pred <- predict(compound_compound_model, CxC_Train)
actual_10CV <- CxC_Train$Activity_Train
matrix_compound_compound_10CV <- table(pred, actual_10CV)
prediction <- predict(compound_compound_model, CxC_Test)
actual_test <- CxC_Test$Activity_Test 
matrix_compound_compound_test <- table(prediction, actual_test)
## Model Performance Evaluation
z <- data.frame(matrix_C_alone_10CV)
y <- data.frame(matrix_C_alone_test)
x <- data.frame(matrix_P_alone_10CV)
w <- data.frame(matrix_P_alone_test)
v <- data.frame(matrix_compound_protein_10CV)
u <- data.frame(matrix_compound_protein_test)
t <- data.frame(matrix_protein_protein_10CV)
s <- data.frame(matrix_protein_protein_test)
r <- data.frame(matrix_compound_compound_10CV)
q <- data.frame(matrix_compound_compound_test)
results <- cbind(z$Freq, y$Freq, x$Freq, w$Freq, v$Freq,
                 u$Freq, t$Freq, s$Freq, r$Freq, q$Freq)
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

Model <- c("10FoldCV Compound Alone",
           "External Compound Alone",
           "10FoldCV Protein Alone",
           "External Protein Alone",
           "10FoldCV Compound-Protein",
           "External Compound-Protein",
           "10FoldCV Protein-Protein",
           "External Protein-Protein",
           "10FoldCV Compound-Compound",
           "External Compound-Compound")
Perf = data.frame(ACC,SENS,SPEC,MCC)
is.num <- sapply(Perf, is.numeric)
Perf[is.num] <- lapply(Perf[is.num], round, 2)
table <- data.frame(Model, Perf)
pdf("Performance_Table_of_PCM_without_SMOTE.pdf")
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
