options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
options(repos=structure(c(CRAN="http://cran.us.r-project.org")))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ggtree", quietly = TRUE)) BiocManager::install("ggtree")
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
if (!requireNamespace("ggridges", quietly = TRUE)) BiocManager::install("ggridges")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) BiocManager::install("RColorBrewer")
if (!requireNamespace("ggplotify", quietly = TRUE)) BiocManager::install("ggplotify")
library("ggtree")
library("RColorBrewer")
library("ggplot2")
library("ggridges")
library("ggplotify")

tree=read.tree("./Tree2GD_out/Phtree.nwk")
summary=read.table("./Tree2GD_out/summarytable.txt",header = T)
gdtype=read.table("./Tree2GD_out/GDtype_stat.txt",header = T)
dollop=read.csv("./dollop.out.scv",header = F)
kaks_file_list=list.files("./sp_kaks_out/")
all_kaks_result=data.frame()
kaks_list=c()
for(F in kaks_file_list){
  sp_name=strsplit(F,"_")[[1]][1]
  kaks_list=c(kaks_list,sp_name)
  tmp=read.table(paste('./sp_kaks_out/',F,sep = ""))
  colnames(tmp)=c("sp","level","gene.1","gene.2","ka","ks","ka.ks","4dtv")
  assign(sp_name,tmp)
  all_kaks_result=rbind(all_kaks_result,get(sp_name))
}


dollop=cbind(summary[,1],dollop)
names(dollop)=c(names(summary)[1],"up","down","all")
dollop$up=paste("+",dollop$up)

gdtype$AABBratio=round(as.numeric(gdtype[[2]]/(gdtype[[2]]+gdtype[[3]]+gdtype[[4]])),2)
t=paste(gdtype$AABBratio*100,"%")
gdtype$describe=t


summary$GD=paste(summary$GD,"of",summary$NUM)
summary$GDratio=as.double(lapply(summary$GDratio, function(x) as.numeric(gsub("\\%", "", x))/100))
summary$GDratio=round(summary$GDratio,2)

tree$edge.length=rep(1,length(tree$edge[,1]))

###################################
dir.create("./R_plotout/")
dir.create("./html_plot_in/")
for(F in kaks_list){
  pdf(paste("./R_plotout/",F,"ks.R.result.pdf",sep = ""),width = 10, height = 10)
  kaks=get(F)
  kaks_p=ggplot(kaks,aes(x=ks,fill=level))+geom_histogram(binwidth = 0.05,alpha = 0.3,colour="black",size=0.1,position ="stack")+xlim(0,4)
  print(kaks_p)
  dev.off()
  test=data.frame()
  node=c()
  for(l in unique(kaks$level)){
    hist_data=table(cut(kaks[which(kaks$level==l&kaks$ks<5),6],breaks=seq(0,5,by=0.02)))
    test=rbind(test,hist_data)
    node=c(node,l)
  }
  rownames(test)=node
  colnames(test)=seq(0,5,by=0.02)[-1]
  write.table(test,paste("./html_plot_in/",F,"_gd2ks.txt",sep = ""),sep = "\t",quote=FALSE)
}



###########################




p=ggtree(tree,size=rel(1),branch.length='none') %<+% summary %<+% dollop %<+% gdtype+geom_tippoint(aes(color=label),size=3)

p2=p+xlim_tree(max(p$data$x)+0.5)+geom_tiplab(offset = 0, hjust = 0,size=7)+geom_rootedge(rootedge = 1)
p3=p2+
  geom_label(aes(x=branch, label=all), fill='grey',label.size = 0.4,size = 3,colour="blue")+
  geom_text(aes(x=branch, label=up), colour='green',size = 3,vjust=-1.5,check_overlap = T)+
  geom_text(aes(x=branch, label=down), colour='red',size = 3,vjust=+2.5)+
  scale_colour_manual(values =palette(rainbow(as.numeric(length(tree$tip.label)))),guide = FALSE)+
  geom_label(aes(label=describe,fill=AABBratio),size = 3,vjust=+4,label.padding = unit(0.03, "lines"),label.size=0.2)+
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(7,'RdYlBu')))(35),position="right",name="Color gradient")+
  geom_label(aes(label=GD,fill= GDratio+0.2),size = 3.5,vjust=2,label.padding = unit(0.1, "lines"),label.size=0.5)+
  theme_classic()

p4=facet_plot(p3,data=all_kaks_result[which(all_kaks_result$ks<3&!is.na(all_kaks_result$ks)),],mapping = aes(x=`ks`,group=label,fill=..density..,height = stat(density)),stat = "density",quantile_lines = TRUE,col="#e94560",size=0.5,scale = 1.5,rel_min_height = 0.01,alpha=15,geom=geom_density_ridges_gradient,panel = 'Ks Plot')%>%facet_labeller(c(Tree = "phylogeny")) %>% facet_widths(c(Tree =1,'Ks Plot'=0.1))



pdf("Tree2GD.result.pdf", width = (max(p$data$x)+1)*3, height = length(tree$tip.label)/2+2)
print(p4+ xlim_expand(c(0, 3), 'Ks Plot'))
dev.off()

