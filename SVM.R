library(e1071)
library(dplyr)



wgd_data=read.table("ks2ratio.txt",header=T,row.names = 1)
wgd_data$Newick_label=rownames(wgd_data)
x=cbind(wgd_data$ks_m,log(wgd_data$ratio))
colnames(x)=c("ks","log_ratio")
type=wgd_data$type
x=data.frame(x, type=as.factor(type))

svm.model <- svm(type~.,data=x,kernel="linear",cost=5,scale = FALSE)
summary(svm.model)
plot(svm.model, x, ks ~ log_ratio)

plot(x[,1:2],col=x[,3],pch=16)
points(x[svm.model$index,c(1,2)],col="blue",cex=2)
w=t(svm.model$coefs)%*%svm.model$SV
b=-svm.model$rho
abline(a=-b/w[1,2],b=-w[1,1]/w[1,2],col="red",lty=5)

pp=ggplot(wgd_data[which(wgd_data$type=="TRUE"),],aes(x=ks_m,y=ratio,label=Newick_label))+
  geom_point(mapping=aes(shape="0"),size=4,colour="red")+
  theme_classic()+
  geom_text(size=5)+
  stat_smooth(color = "blue",formula=y ~log(x),method = "lm",se=FALSE,lty=2)+
  geom_point(data=wgd_data[which(wgd_data$type!="TRUE"),],mapping=aes(x=ks_m,y=ratio,label=Newick_label),shape=3,size=2,colour="blue",show.legend=FALSE)

