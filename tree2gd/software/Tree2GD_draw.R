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
library("grid")
library("scales")

tree=read.tree("./Tree2GD_out/Phtree.nwk")
summary=read.table("./Tree2GD_out/summarytable.txt",header = T)
gdtype=read.table("./Tree2GD_out/GDtype_stat.txt",header = T)
dollop=read.csv("./dollop.out.scv",header = F)
kaks_file_list=list.files("./sp_kaks_out/")
all_kaks_result=data.frame()
kaks_list=c()

for(F in kaks_file_list){
  sp_name=strsplit(F,'[.]')[[1]][1]
  kaks_list=c(kaks_list,sp_name)
  tmp=read.table(paste('./sp_kaks_out/',F,sep = ""))
  colnames(tmp)=c("sp","level","gene.1","gene.2","ka","ks","ka.ks","4dtv")
  assign(sp_name,tmp)
  t=get(sp_name)
  t$Node=t$level
  t$sum=t$ks
  for(l in summary$Newick_label){
    if(l %in% t$level){
    a=table(cut(t[which(t$level==l&t$ks<3),6],breaks=seq(0,3,by=0.01)))
    tmp=which(a==a[which.max(a)],arr.ind=T)[1]*0.01
    num=length(t[which(t$level==l),]$ks)
    t[which(t$level==l),]$Node=paste(t[which(t$level==l),]$level,"(",num,";",tmp,")")
    t[which(t$level==l),]$sum=num
    }
  }
  all_kaks_result=rbind(all_kaks_result,get(sp_name))
  assign(sp_name,t)
}


dollop=cbind(summary[,1],dollop)
names(dollop)=c(names(summary)[1],"up","down","all")
dollop$up=paste("+",dollop$up)

gdtype$AABBratio=round(as.numeric(gdtype[[2]]/(gdtype[[2]]+gdtype[[3]]+gdtype[[4]])),2)
t=paste(gdtype$AABBratio*100,"%")
gdtype$describe=t


summary$GD=paste(summary$GD,"/",summary$NUM)
summary$GDratio=as.double(lapply(summary$GDratio, function(x) as.numeric(gsub("\\%", "", x))/100))
summary$GDratio=round(summary$GDratio,2)

tree$edge.length=rep(1,length(tree$edge[,1]))

###################################
dir.create("./R_plotout/")
dir.create("./html_plot_in/")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#0f4c75","#d7385e","#d2e603","#ff4b5c","#1f3c88","#03c4a1","#f39233","#7579e7")
for(F in kaks_list){
  pdf(paste("./R_plotout/",F,".ks.R.result.pdf",sep = ""),width = 15, height = 6)
  kaks=get(F)
  kaks_p=ggplot(kaks[which(kaks$ks<3),],aes(x=ks,fill=Node))+geom_histogram(binwidth = 0.05,alpha = 0.8,colour="black",size=0.3,position ="stack",show.legend=FALSE)+
    labs(x = "Ks",y = "Number of Gene Pairs",title = paste(" Ks-plot of paralogs from ",F,sep = ""))+
    scale_fill_manual(values=cbPalette)+
    scale_colour_manual(values=cbPalette)+scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
    theme_classic()+theme(legend.key.size = unit(10, "pt"),legend.position=c(.9,.7),legend.key = element_blank(),legend.background = element_blank(),legend.text = element_text(size = 6),legend.title = element_text(size = 8))
  ks_p=ggplot(kaks[which(kaks$ks<3),],aes(x=ks,y=Node,fill=Node,alpha=rescale(sum,c(0,1))))+
    geom_density_ridges_gradient(scale=3,size=0.3,rel_min_height=0.00,show.legend=FALSE)+
    scale_fill_manual(values=cbPalette)+
    scale_y_discrete(position = "right")+ theme_minimal()+theme(title=element_text(size=7))+
    labs(x = "Ks",title ="Density Ridges Gradient of Gene Pairs Number",y="")
  vp <- viewport(width = 0.2, height = 0.4, x = 0.95,y = 0.9,just=c("right","top"))
  print(kaks_p)
  print(ks_p,vp=vp)
  dev.off()
  test=data.frame()
  node=c()
  for(l in unique(kaks$level)){
    hist_data=table(cut(kaks[which(kaks$level==l&kaks$ks<5),6],breaks=seq(0,5,by=0.05)))
    test=rbind(test,hist_data)
    node=c(node,l)
  }
  rownames(test)=node
  colnames(test)=seq(0,5,by=0.05)[-1]
  write.table(test,paste("./html_plot_in/",F,"_gd2ks.txt",sep = ""),sep = "\t",quote=FALSE)
}



###########################

summary$pho=substr(summary$Newick_label,7,10)
p=ggtree(tree,size=rel(0.5),branch.length='none') %<+% summary %<+% dollop %<+% gdtype+geom_tippoint(color="#56556e",size=3)

p2=p+xlim_tree(max(p$data$x)+1.3)+geom_tiplab(offset = 0, hjust = 0,size=7)+geom_rootedge(rootedge = 1)+
  geom_label2(aes(label=pho,subset = !is.na(as.numeric(pho))))

p3=p2+
  geom_label(aes(x=branch, label=all), fill='grey',label.size = 0.2,size = 3,colour="blue")+
  geom_text(aes(x=branch, label=up), colour='green',size = 3,vjust=-1.5,check_overlap = T)+
  geom_text(aes(x=branch, label=down), colour='red',size = 3,vjust=+2.5)+
  scale_colour_manual(values =palette(rainbow(as.numeric(length(tree$tip.label)))),guide = FALSE)+
  geom_label(aes(label=describe),size = 3,vjust=+4,fill="#a5ecd7",label.padding = unit(0.03, "lines"),label.size=0.2,show.legend=FALSE)+
  #scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(7,'RdYlBu')))(35),position="right",name="Color gradient")+
  geom_label(aes(label=GD),size = 3.5,fill= "#ff9a76",vjust=2,label.padding = unit(0.1, "lines"),label.size=0.2,show.legend=FALSE)+
  theme_classic()

p4=p3+geom_facet(data=all_kaks_result[which(all_kaks_result$ks<2&!is.na(all_kaks_result$ks)),],mapping = aes(x=`ks`,group=label,height = stat(density)),stat = "density",col="#e94560",size=0.5,scale = 2,rel_min_height = 0.01,alpha=15,geom=geom_density_ridges_gradient,panel = 'Ks Plot')+
geom_facet(data=all_kaks_result[which(all_kaks_result$`4dtv`<1&!is.na(all_kaks_result$`4dtv`)),],mapping = aes(x=`4dtv`,group=label,height = stat(density)),stat = "density",col="#e94560",size=0.5,scale = 2,rel_min_height = 0.01,alpha=15,geom=geom_density_ridges_gradient,panel = '4dtv Plot')


p5=p4%>%facet_labeller(c(Tree = "phylogeny")) %>% facet_widths(c(Tree =1,'4dtv Plot'=0.1,'Ks Plot'=0.1))



pdf("Tree2GD.result.pdf", width = (max(p$data$x)+1)*1.7, height = length(tree$tip.label)/2+2)
print(p5)
dev.off()
