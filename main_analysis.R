setwd("C:/Users/aless/Documents/UniTn/SecondSemester/NetworkDataAnalysis/NitrogenStarvation/")
library('GEOquery')
library('FactoMineR')
library('factoextra')
library('plotly')
library('genefilter')
library('randomForest')
library("MASS")
library("pROC")
library(RColorBrewer)
library(glmnet)
library(e1071)
library('rScudo')
library("igraph")
library('caret')
library(shape)

gse<-getGEO('GSE110171')
str(gse)

# selecting "`tissue:ch1`", "`treatment:ch1`", "characteristics_ch1.2"
pheno<-gse$GSE110171_series_matrix.txt.gz@phenoData@data[, c(43,44,12)]
pheno
pheno$`treatment:ch1`<- ifelse(pheno$`treatment:ch1`=='9 mM nitrate', 'highN', 'lowN')
# selecting useful information for annotation
anno<-gse$GSE110171_series_matrix.txt.gz@featureData@data[, c(8:17)]
str(anno)
# extracting expression data
ex <- exprs(gse[[1]])
dim(ex)
# relaxing memory
rm(gse)

# removing zero variance genes
ex <- ex[which(apply(ex, 1, var) != 0), ]
dim(ex)

# filter genes by t-test
tt <- rowttests(ex, as.factor(pheno$`treatment:ch1`))
keepers <- which(p.adjust(tt$p.value, method = 'bonferroni')<0.1)
keepers <- which(tt$p.value<0.1)
filt_ex <- ex[keepers,]
dim(filt_ex)

# box plot
boxplot(ex, ylab='expression values', xaxt = 'n', main='Unnormalized data')
boxplot(scale(ex), ylab='expression values', xaxt = 'n', main='Normalized data')
boxplot(filt_ex)
boxplot(scale(filt_ex))

# pca
pca<-PCA(t(ex), graph = F)
summary(pca)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_ind(pca, geom = 'point', col.ind = as.factor(pheno$`treatment:ch1`), addEllipses = TRUE, legend.title = 'Condition', ellipse.level=0.69, mean.point = F, res=1000)

filt_pca <- PCA(t(filt_ex), graph = F)
fviz_pca_ind(filt_pca, geom = 'point', col.ind = as.factor(pheno$`treatment:ch1`), addEllipses = TRUE, legend.title = 'Condition', ellipse.level=0.69, mean.point = F, res=1000)



#k-means
##k selection by silhouette
fviz_nbclust(t(ex), kmeans, method = 'silhouette')
fviz_nbclust(t(ex), kmeans, method = 'wss')
kmeans_result<-kmeans(t(ex), centers = 4, nstart = 25)
fviz_cluster(kmeans_result, data = t(ex),
             palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot",
             labelsize = 0,
             pointsize = 2,
             res = 1000
)

#hierarchical
res <- hcut(t(ex), k = 4, stand = TRUE)
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

#Sampling train and test sets
rm(.Random.seed, envir=globalenv())

which_shoot_low <- pheno$`tissue:ch1`=='shoot' & pheno$`treatment:ch1`=='lowN'
which_shoot_high <- pheno$`tissue:ch1`=='shoot' & pheno$`treatment:ch1`=='highN'
which_root_low <- pheno$`tissue:ch1`=='root' & pheno$`treatment:ch1`=='lowN'
which_root_high <- pheno$`tissue:ch1`=='root' & pheno$`treatment:ch1`=='highN'

which_train_shoot_low <- sample(1:nrow(pheno[which_shoot_low, ]), 2*nrow(pheno[which_shoot_low, ])/3)
which_train_shoot_high <- sample(1:nrow(pheno[which_shoot_high, ]), 2*nrow(pheno[which_shoot_high, ])/3)
which_train_root_low <- sample(1:nrow(pheno[which_root_low, ]), 2*nrow(pheno[which_root_low, ])/3)
which_train_root_high <- sample(1:nrow(pheno[which_root_high, ]), 2*nrow(pheno[which_root_high, ])/3)

which_test_shoot_low <- setdiff(1:nrow(pheno[which_shoot_low, ]), which_train_shoot_low)
which_test_shoot_high <- setdiff(1:nrow(pheno[which_shoot_high, ]), which_train_shoot_high)
which_test_root_low <- setdiff(1:nrow(pheno[which_root_low, ]), which_train_root_low)
which_test_root_high <- setdiff(1:nrow(pheno[which_root_high, ]), which_train_root_high)

train <- cbind(ex[, which_root_high][, which_train_root_high], ex[, which_shoot_high][, which_train_shoot_high], ex[, which_root_low][, which_train_root_low], ex[, which_shoot_low][, which_train_shoot_low])
test <- cbind(ex[, which_root_high][, which_test_root_high], ex[, which_shoot_high][, which_test_shoot_high], ex[, which_root_low][, which_test_root_low], ex[, which_shoot_low][, which_test_shoot_low])
train_pheno <- rbind(pheno[which_root_high, ][which_train_root_high, ], pheno[which_shoot_high, ][which_train_shoot_high, ], pheno[which_root_low, ][which_train_root_low, ], pheno[which_shoot_low, ][which_train_shoot_low, ])
test_pheno <- rbind(pheno[which_root_high, ][which_test_root_high, ], pheno[which_shoot_high, ][which_test_shoot_high, ], pheno[which_root_low, ][which_test_root_low, ], pheno[which_shoot_low, ][which_test_shoot_low, ])

filt_train <- cbind(filt_ex[, which_root_high][, which_train_root_high], filt_ex[, which_shoot_high][, which_train_shoot_high], filt_ex[, which_root_low][, which_train_root_low], filt_ex[, which_shoot_low][, which_train_shoot_low])
filt_test <- cbind(filt_ex[, which_root_high][, which_test_root_high], filt_ex[, which_shoot_high][, which_test_shoot_high], filt_ex[, which_root_low][, which_test_root_low], filt_ex[, which_shoot_low][, which_test_shoot_low])

train_cond <- cbind(as.data.frame(t(train)), as.factor(train_pheno$`treatment:ch1`))
colnames(train_cond)[ncol(train_cond)] <- 'CONDITION'
test_cond <- cbind(as.data.frame(t(test)), as.factor(test_pheno$`treatment:ch1`))
colnames(test_cond)[ncol(test_cond)] <- 'CONDITION'

filt_train_cond <- cbind(as.data.frame(t(filt_train)), as.factor(train_pheno$`treatment:ch1`))
colnames(filt_train_cond)[ncol(filt_train_cond)] <- 'CONDITION'
filt_test_cond <- cbind(as.data.frame(t(filt_test)), as.factor(test_pheno$`treatment:ch1`))
colnames(filt_test_cond)[ncol(filt_test_cond)] <- 'CONDITION'

#RF
set.seed(42)
rf <- randomForest(x=t(train), y=as.factor(train_pheno$`treatment:ch1`), ntree=10000)
rf_pred <- predict(rf, t(test))
rf_confusion = table(rf_pred, test_pheno$`treatment:ch1`)

filt_rf <- randomForest(x=t(filt_train), y=as.factor(train_pheno$`treatment:ch1`), ntree=10000)
filt_rf_pred <- predict(filt_rf, t(filt_test))
filt_rf_confusion = table(filt_rf_pred, test_pheno$`treatment:ch1`)

#selection of most important genes in RF
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp, decreasing = T)
plot(c(1:nrow(ex)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16)
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25]
t <- is.element(rownames(ex),gn.25)
sig.ex <- ex[t,]

#Check differences with RF based on filtered genes
filt_imp = abs(filt_rf$importance[,])
filt_imp = names(filt_imp[order(filt_imp, decreasing = T)])
intersect(filt_imp[1:25], gn.25)

#heatmap of expression value most important genes in RF
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.ex) <- pheno$`treatment:ch1`
csc = rep(hmcol[50],72)
csc[pheno$`treatment:ch1`=='9 mM nitrate'] <- hmcol[200]
heatmap(sig.ex, scale="row", col=hmcol, ColSideColors=csc)

# lda
mod <- lda(CONDITION ~ ., data=filt_train_cond, prior = c(0.5,0.5))
plot(mod)
pred_lda <- predict(mod, filt_test_cond)
lda_genes <- coef(mod)[,1][order(coef(mod)[,1], decreasing = T)]
lda_confusion <- table(pred_lda$class, test_pheno$`treatment:ch1`)
roc_lda <- plot.roc(as.numeric(pred_lda$class),
                    as.numeric(as.factor(test_pheno$`treatment:ch1`)))
auc_lda <- auc(as.numeric(pred_lda$class),
               as.numeric(as.factor(test_pheno$`treatment:ch1`)))

#lasso
cfit=cv.glmnet(t(train),as.factor(train_pheno$`treatment:ch1`),standardize=FALSE,family="binomial")
plot(cfit)
pred_lasso <- predict(cfit,t(test), type="class", s= cfit$lambda.min)
lasso_confusion <- table(pred_lasso, test_pheno$`treatment:ch1`)
lasso_genes <- rownames(coef(cfit, s=cfit$lambda.min))[which(coef(cfit, s=cfit$lambda.min) != 0)][2:4]
roc_lasso <- plot.roc(as.numeric(as.factor(pred_lasso)),
                      as.numeric(as.factor(test_pheno$`treatment:ch1`)))
auc_lasso <- auc(as.numeric(as.factor(pred_lasso)),
                 as.numeric(as.factor(test_pheno$`treatment:ch1`)))

#svm
svmfit <- svm(CONDITION ~ ., data = filt_train_cond, kernel = "linear", cost = 10, scale = FALSE, type = "C")
svmpred <- predict(svmfit, filt_test_cond)
svmconf <- table(svmpred, test_pheno$`treatment:ch1`)

# rScudo
trainRes <- scudoTrain(train, groups = as.factor(train_pheno$`treatment:ch1`),
                       nTop = 75, nBottom = 75, alpha = 0.1, pAdj = 'bonferroni')
trainNet <- scudoNetwork(trainRes, N = 0.28)
scudoPlot(trainNet, vertex.label = NA, res=1000)
testRes <- scudoTest(trainRes, test, as.factor(test_pheno$`treatment:ch1`), 
                     nTop = 75,nBottom = 75)
testNet <- scudoNetwork(testRes, N = 0.30)
scudoPlot(testNet, vertex.label = NA)
write_graph(trainNet, 'rscudoTrain.graphml', 'graphml')

testClust <- igraph::cluster_spinglass(testNet, spins = 2)
plot(testClust, testNet, vertex.label = NA)
# feature selection was performed before, the internal way is buggy
classRes <- scudoClassify(filt_train, filt_test, N = 0.28,
                          nTop = 75, nBottom = 75,
                          trainGroups = as.factor(train_pheno$`treatment:ch1`), featureSel = F)
conf_rscudo <- caret::confusionMatrix(classRes$predicted, as.factor(test_pheno$`treatment:ch1`))
roc_rscudo <- plot.roc(as.numeric(classRes$predicted),
                       as.numeric(as.factor(test_pheno$`treatment:ch1`)))
auc_rscudo <- auc(as.numeric(classRes$predicted),
                  as.numeric(as.factor(test_pheno$`treatment:ch1`)))


# Model selection
ctrl <- trainControl(method =  "repeatedcv", number = 8, repeats = 10, 
                     classProbs = T, summaryFunction = twoClassSummary)

cl <- parallel::makePSOCKcluster(6)
doParallel::registerDoParallel(cl)

# best lda
fit.lda <- caret::train(CONDITION~., data=filt_train_cond, method="lda",
                        trControl=ctrl)
preds_lda <- predict(fit.lda, filt_test_cond)
conf_lda <-table(preds_lda, test_pheno$`treatment:ch1`)
auc_lda <- auc(as.numeric(as.factor(preds_lda)),
               as.numeric(as.factor(test_pheno$`treatment:ch1`)))

# best rf
fit.rf <- caret::train(CONDITION~., data=filt_train_cond, ntree=1000, method="rf",
                       trControl=ctrl)
preds_rf <- predict(fit.rf, filt_test_cond)
conf_rf <-table(preds_rf, test_pheno$`treatment:ch1`)
auc_rf <- auc(as.numeric(as.factor(preds_rf)),
              as.numeric(as.factor(test_pheno$`treatment:ch1`)))
imp <- (abs(fit.rf$finalModel$importance[,]))
o <- order(imp, decreasing = T)
plot(imp[o][1:200], ylab = 'Importance')
Arrows(170,1.05,188,0.15,lwd=2, arr.type="triangle", col = 'Purple')
imp <- gsub("`", "", names(imp[o]))
rf_genes <- imp[1:188]

# best lasso
fit.lasso <- train(CONDITION~., data=filt_train_cond,
                   method = "glmnet",
                   family = "binomial",
                   tuneGrid = expand.grid(alpha = 1,
                                          lambda = seq(0,1,by=0.05)),
                   trControl = ctrl)
preds_lasso <- predict(fit.lasso, filt_test_cond)
conf_lasso <-table(preds_lasso, test_pheno$`treatment:ch1`)
auc_lasso <- auc(as.numeric(as.factor(preds_lasso)),
                 as.numeric(as.factor(test_pheno$`treatment:ch1`)))
lasso_genes <- rownames(coef(fit.lasso$finalModel, 
                             s=fit.lasso$finalModel$lambdaOpt))[
                               which(coef(fit.lasso$finalModel, 
                                          s=fit.lasso$finalModel$lambdaOpt) != 0)]
lasso_genes <- gsub("`", "", lasso_genes[-1])

# best svm
fit.svm <- train(CONDITION~., data=filt_train_cond, method='svmLinear2', 
                 type='C', trControl=ctrl)
preds_svm <- predict(fit.svm, filt_test_cond)
conf_svm <-table(preds_svm, test_pheno$`treatment:ch1`)
auc_svm <- auc(as.numeric(as.factor(preds_svm)),
               as.numeric(as.factor(test_pheno$`treatment:ch1`)))

# best rScudo
model <- scudoModel(nTop = c(50, 75, 100), nBottom = c(50, 75, 100), N = c(0.20, 0.30, 0.40))
fit.rscudo <- train(CONDITION~., data=filt_train_cond, method = model, trControl=ctrl)

preds_rscudo <- scudoClassify(filt_train, filt_test, N = fit.rscudo$bestTune$N,
                              nTop = fit.rscudo$bestTune$nTop,
                              nBottom = fit.rscudo$bestTune$nBottom, 
                              trainGroups = as.factor(train_pheno$`treatment:ch1`), 
                              featureSel = F)$predicted
conf_rscudo <- table(preds_rscudo, test_pheno$`treatment:ch1`)
auc_rscudo <- auc(as.numeric(preds_rscudo),
                  as.numeric(as.factor(test_pheno$`treatment:ch1`)))
best_rscudo <- scudoTrain(filt_train, groups = as.factor(train_pheno$`treatment:ch1`),
                          nTop = fit.rscudo$bestTune['nTop'][[1]], 
                          nBottom = fit.rscudo$bestTune['nBottom'][[1]], 
                          featureSel = F)
rscudo_genes <- unique(c(consensusUpSignatures(best_rscudo)[,1], consensusUpSignatures(best_rscudo)[,2], 
                         consensusDownSignatures(best_rscudo)[,1], consensusUpSignatures(best_rscudo)[,2]))

rf_rscudo_genes <- union(rscudo_genes, rf_genes)
parallel::stopCluster(cl)

res<-resamples(list(LDA=fit.lda, LASSO=fit.lasso, SVM=fit.svm, RF=fit.rf, rScudo=fit.rscudo))
ggplot(res) + labs(y = "AUC") + ggtitle('Models comparison')

write.csv(rscudo_genes, file = 'rscudo_probes.txt', quote = F, row.names = F, col.names = F)
write.csv(rf_genes, file = 'rf_probes.txt', quote = F, row.names = F, col.names = F)
write.csv(lasso_genes, file = 'lasso_probes.txt', quote = F, row.names = F, col.names = F)
write.csv(rf_rscudo_genes, file = 'rf_rscudo_probes.txt', quote = F, row.names = F, col.names = F)
