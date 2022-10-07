options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
options(repos=structure(c(CRAN="http://cran.us.r-project.org")))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("e1071", quietly = TRUE)) BiocManager::install("e1071")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr")

library(e1071)
library(dplyr)
library(ggplot2)
lm_eqn <- function(df,x,y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*", p-value="~p, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3),
                        p = format(as.numeric(summary(m)$coefficients[,4][2]), digits = 3)))
  as.character(as.expression(eq));
}
                                                                  
a = list.files("sp_kaks_out/")                                                      
dir = paste("./sp_kaks_out/",a,sep="")                                     
n = length(dir)                                                                
merge.data = read.csv(file = dir[1],sep="\t",header = F)  
for (i in 2:n){
  new.data = read.csv(file = dir[i], sep="\t",header = F)
  merge.data = rbind(merge.data,new.data)
}
merge.data=merge.data[which(merge.data$V6<=2),]
kaks_table=aggregate(merge.data[,6],by=list(type=merge.data[,2]),mean)


ratio_data=read.table("Tree2GD_out/summarytable.txt",header = T)
wgd_data=merge(x = ratio_data,y = kaks_table,by.x = "Newick_label",by.y = "type")
wgd_data=wgd_data[grep("^phyto", wgd_data[,1]),]  
wgd_data$GDratio=as.numeric(sub("%", "", wgd_data$GDratio))/100
wgd_data[is.na(wgd_data)]<-0


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
  x=data.frame(cbind(wgd_data$x,log(wgd_data$GDratio)))
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
