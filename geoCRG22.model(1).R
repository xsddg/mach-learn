library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)      
inputFile="normalize.txt"     
geneFile="interGenes.txt"     
setwd("C:\\biowolf\\geoCRG\\22.model")  

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)

mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)


p_fun=function(object, newdata){
	predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="Control", 0, 1)

explainer_rf=explain(mod_rf, label = "RF",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_rf=model_performance(explainer_rf)
explainer_svm=explain(mod_svm, label = "SVM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_svm=model_performance(explainer_svm)
explainer_xgb=explain(mod_xgb, label = "XGB",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)
explainer_glm=explain(mod_glm, label = "GLM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_glm=model_performance(explainer_glm)

pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()

pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
legend('bottomright',
	   c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
	     paste0('SVM: ',sprintf("%.03f",roc2$auc)),
	     paste0('XGB: ',sprintf("%.03f",roc3$auc)),
	     paste0('GLM: ',sprintf("%.03f",roc4$auc))),
	   col=c("red","blue","green","yellow"), lwd=2, bty = 'n')
dev.off()

importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)

pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()

geneNum=5     
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)
