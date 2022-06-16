library(e1071)
library(dplyr)
library(ggplot2)
#######################
lm_eqn <- function(df,x,y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*", p-value="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(as.numeric(summary(m)$coefficients[,4][2]), digits = 3)))
  as.character(as.expression(eq));
}
######################

wgd_data=read.csv("SVM_data.csv",header=T,row.names = 1)
wgd_data$Newick_label=rownames(wgd_data)
x=cbind(wgd_data$median,log(wgd_data$GDratio))
##The original GD median with scale is used by default and can be replaced with covariance support results
colnames(x)=c("ks","log_ratio")

if("type"  %in% colnames(wgd_data)){
type=wgd_data$type
x=data.frame(x, type=as.factor(type))

svm.model <- svm(type~.,data=x,kernel="linear",cost=5,scale = FALSE)
summary(svm.model)
plot(svm.model, x, ks ~ log_ratio)



w=t(svm.model$coefs)%*%svm.model$SV
b=-svm.model$rho
x$node=wgd_data$Newick_label

p=ggplot(x[which(x$type=="TRUE"),],aes(x=as.numeric(ks),y=as.numeric(log_ratio),label=node))+
  geom_point(size=4,colour="red")+
  theme_classic()+
  geom_text(size=5,vjust = "inward", hjust = "inward")+
  stat_smooth(color = "blue",formula=y ~x,method = "lm",lty=2)+
  labs(title ="", x="Ks",y="ln(GD_raito)",size=5)+
  geom_point(data=x[which(x$type!="TRUE"),],mapping=aes(x=as.numeric(ks),y=as.numeric(log_ratio),label=node),shape=3,size=3,colour="blue",show.legend=FALSE)
p1 <- p + 
  geom_text(x = 1.5, y = -0.5, label =lm_eqn(x[which(x$type=="TRUE"),],as.numeric(x[which(x$type=="TRUE"),"ks"]),as.numeric(x[which(x$type=="TRUE"),"log_ratio"])), 
            parse = TRUE,colour="blue",size=5)
p2 <- p1 +geom_abline(aes(intercept=-b/w[1,2],slope=-w[1,1]/w[1,2]),lty=2,colour="red",size=1) +
  geom_text(x =0.5, y = -2.5,label = paste("SVM: y = ",round(-w[1,1]/w[1,2],2),"x",round(-b/w[1,2],2)),colour="red",size=5)
print(p2)

}else {
  x=as.data.frame(cbind(wgd_data$median,log(wgd_data$GDratio)))
  colnames(x)=c("ks","log_ratio")
  x[,"type"]=x[,"ks"]
  for(i in 1:nrow(x)){
    if(x[i,"log_ratio"]>x[i,"ks"]*(-1.28)-1.3){
      x[i,"type"]="TRUE";
    } else {
      x[i,"type"]="FALSE";
    }
  }
  ###########No know WGD loci were given, using the classification planes obtained in the 68 plants data.
  a=-1.28
  b=-1.3
  x$node=wgd_data$Newick_label
  
  p=ggplot(x[which(x$type=="TRUE"),],aes(x=as.numeric(ks),y=as.numeric(log_ratio),label=node))+
    geom_point(size=4,colour="red")+
    theme_classic()+
    geom_text(size=5,vjust = "inward", hjust = "inward")+
    stat_smooth(color = "blue",formula=y ~x,method = "lm",lty=2)+
    labs(title ="", x="Ks",y="ln(GD_raito)",size=5)+
    geom_point(data=x[which(x$type=="FALSE"),],mapping=aes(x=as.numeric(ks),y=as.numeric(log_ratio),label=node),shape=3,size=3,colour="blue",show.legend=FALSE)
  p1 <- p + 
    geom_text(x = 1.5, y = -0.5, label =lm_eqn(x[which(x$type=="TRUE"),],as.numeric(x[which(x$type=="TRUE"),"ks"]),as.numeric(x[which(x$type=="TRUE"),"log_ratio"])), 
              parse = TRUE,colour="blue",size=5)
  p2 <- p1 +geom_abline(aes(intercept=b,slope=a),lty=2,colour="red",size=1) +
    geom_text(x =0.5, y = -2.5,label = paste("SVM: y = "," +",a,"x",b),colour="red",size=5)
  print(p2)
  
}
