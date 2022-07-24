
#limma包基因差异分析----------------------------------------------------------------------------------
logFCcutoff=1.5
PValue=0.05

library(limma)
library(tidyverse)
#首先读入数据，整理数据格式，normalize数据
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=normalizeBetweenArrays(as.matrix(rt))
#rt=log2(rt+1)

#读入临床信息矩阵
clinical<-read.table("clinical.txt",sep="\t",header=T,check.names=F) 

##获取分组信息
group<- clinical %>% dplyr::filter(submitter_id %in% colnames(rt))%>% dplyr::select("submitter_id","group") 

rt<-rt[,group$submitter_id]
identical(group$submitter_id,colnames(rt))

##写出相同submitter_id的表达矩阵
out=rbind(gene=colnames(rt),rt)
write.table(out,file="./data/symbol_filter(normalize).txt",sep="\t",quote=F,col.names=F)

##写出相同submitter_id的临床矩阵
clinical_filter<-clinical %>% filter(submitter_id %in% colnames(rt))

table(clinical_filter$group)#肿瘤和正常组织

out=rbind(gene=colnames(clinical_filter),clinical_filter)
write.table(out,file="./data/clinical_filter(N72T508).txt",sep="\t",quote=F,col.names = F)

##基因样本过滤
table(group$group)
keep <- rowSums(rt > 0) >= 72
rt <- rt[keep, ]

##差异表达分析
Type=group$group
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("Normal","Tumor")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(Tumor-Normal,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

##输出所有基因以及符合logFC值及P值标准的差异表达基因

all=topTable(fit2,number=Inf,adjust.method="holm")
out=rbind(gene=colnames(all),all)
write.table(out,file="./data/all.txt",sep="\t",quote=F,col.names=F)

diffSig = all[all$adj.P.Val < PValue & abs(all$logFC)>logFCcutoff,]
out=rbind(gene=colnames(diffSig),diffSig)
write.table(out, file="./data/diffSig.txt",sep="\t",quote=F,col.names = F)

#### 火山图 -----------------------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggsci)
##准备火山图数据 
DEG_all <- all %>% dplyr::select(adj.P.Val,logFC) %>%
  dplyr::mutate(direction = factor(ifelse(adj.P.Val < PValue & abs(logFC) > logFCcutoff,
        ifelse(logFC >logFCcutoff, "Up", "Down"),"NS"),levels=c('Up','Down','NS')))

##读取特异性相关表达基因
gene<- read.table("GeneName.txt",header = T)

##特异性相关表达基因分别与训练集和验证集基因取交集
test=read.table("validation_set.txt",header = T,check.names = F)

length(intersect(rownames(all), gene$GeneName)) #357
length(intersect(test$GeneName, gene$GeneName)) #351
##特异性相关表达基因与训练集差异表达基因取交集，再与验证集基因取交集

volcano <- DEG_all %>%rownames_to_column("symbol") %>% dplyr::filter(symbol %in% test$GeneName) %>% 
           dplyr::filter(symbol %in% gene$GeneName)
           
#### 绘制火山图
ggplot(data = volcano, aes(x = logFC, y = -log10(adj.P.Val), colour = direction)) + 
  geom_point(alpha = 2,size=3) +
   scale_color_manual(values =c("#BC3C2999", "#0072B599", "#808080")) + 
  theme_bw() +
  theme(panel.grid =element_blank(),
        axis.text.y = element_text(size=17),
        axis.text.x = element_text(size=17),
        axis.title = element_text(size=18),
        plot.title = element_text(hjust = 0.5,size=20),
        strip.text = element_text(size=19),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        text = element_text(family ="sans",face ="bold"))+
  ylab(expression(-log[10]("adjust P Value"))) +
  xlab(expression(log[2]("Fold Change"))) +
  xlim(-4.5, 4.5) +
  ylim(0, 180) +
  geom_vline(xintercept = c(-logFCcutoff,logFCcutoff), 
             lty = 2,
             col = "black",
             lwd = 0.6) +
  geom_hline(yintercept = -log10(PValue),
             lty = 2,
             col = "black",
             lwd = 0.6)  
library(export)
graph2pdf(file="./Figure/volcano.pdf", width = 8,height = 5)       
dev.off()

##差异表达基因与特异性相关表达基因取交集------------------------------------------------------
library(VennDiagram)
library(ggsci)

DEGs <- intersect(rownames(diffSig), test$GeneName) 
length(intersect(DEGs,gene$GeneName))

venn.diagram(x = list("DEGs" = DEGs,
                      "AAMRGs" = gene$GeneName),
             filename = "./Figure/Venn.tiff",
             fill = pal_nejm("default", alpha = 0.6)(2),
             col=pal_nejm("default")(2),
             category.names = c("DEGs","AAM-related genes"),
             height = 380 , 
             width = 380 , 
             resolution = 300,
             compression = "lzw",
             lwd = 1,
             scaled = F,
             main = "Venn Diagram",
             main.cex = 0.3,
             main.pos = c(0.5,0.1),
             cat.cex=0.3,
             cat.pos = c(-10, 10),
             cat.dist = c(0.055, 0.055),cex = 0.3,fontfamily = "sans",
             cat.fontfamily = "sans",main.fontfamily = "sans",
             fontface = "bold",cat.fontface = "bold",main.fontface = "bold")

##特异性相关差异表达基因热图--------------------------------------------
diffName<-intersect(DEGs, gene$GeneName)
write.table(diffName,"./data/diffgene.txt",sep = "\t")

heatmap<-rt %>% as.data.frame()%>% .[diffName,] 
##以logFC排序
res_FC<- DEG_all %>%  rownames_to_column("symbol")%>% dplyr::filter(symbol %in% diffName) %>% arrange(-logFC)
##注释信息
annotation<-group %>% arrange(group)
annotation_col= annotation %>% column_to_rownames("submitter_id")
##表达矩阵
heatmap1<-heatmap %>% .[res_FC$symbol,] %>% .[,annotation$submitter_id] %>% as.matrix()

##画热图
library(ComplexHeatmap)
library(circlize)
heatmap2<-t(scale(t(heatmap1)))
annotation$group<-as.factor(annotation$group)
table(annotation$group)

Cluster <- c('#21b6af','#eeba4d')
names(Cluster) <- levels(annotation$group)

Top = HeatmapAnnotation(Risk=annotation$group,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                     title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col= list(Risk=Cluster),
                        show_annotation_name = F,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))


Heatmap(heatmap2,name='Z-score',
        top_annotation = Top,
        cluster_rows = FALSE,
        col=colorRamp2(c(-2,0,2),c("#0072B599","white","#BC3C2999")),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,72),rep(2,508)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                  title_gp = gpar(fontsize = 10, fontface = "bold"))) 
library(export)
graph2tif(file = "./Figure/heatmap2.tiff",width = 10, height = 8)

##单因素COX回归##--------------------------------------------------------------------------

##获取肿瘤患者相关基因表达以及临床信息矩阵
test1<-read.table("./data/symbol_filter(normalize).txt",sep = "\t",header = T,check.names = F)
test2<-read.table("clinical_T503.txt",sep = "\t",header = T)

##删除正常样本，获取肿瘤患者相关基因表达矩阵
symbol_only_tumor<- test1 %>% as.data.frame()%>%.[,c("gene",test2$submitter_id)]

out<-rbind(ID=colnames(symbol_only_tumor),symbol_only_tumor)
write.table(out,"./data/symbol_T503.txt",sep="\t",quote=F,col.names= F,row.names = F)

KIRC_exp<-symbol_only_tumor %>%column_to_rownames("gene")%>%.[diffName,] %>%
           t() %>% as.data.frame() %>% rownames_to_column("submitter_id")

##获取合并之后肿瘤患者的基因表达以及临床信息矩阵
KIRC_lasso<-test2 %>% inner_join(KIRC_exp,by="submitter_id")

out<-rbind(submitter_id=colnames(KIRC_lasso),KIRC_lasso)
write.table(out,file = "./data/KIRC_lasso.txt",sep="\t",quote = F,col.names = F,row.names = F)

##通过单因素COX回归筛选与预后相关的特异性差异表达基因
library(survival)
library(plyr)
library(tidyverse)

KIRC_lasso<-read.table("./data/KIRC_lasso.txt",sep = "\t",header = T) 

KIRC_lasso$fustat <- as.numeric(as.character(KIRC_lasso$fustat)) ## 设置数据格式
BaSurv <- Surv(time = KIRC_lasso$futime, event = KIRC_lasso$fustat) ## futime为生存时间，fustat生存状态

UniCox <- function(x){ ## 构建一个R function 便于后期调用
  FML <- as.formula(paste0('BaSurv~',x)) ## 构建生存分析公式
  GCox <- coxph(FML, data = data) ## Cox分析
  GSum <- summary(GCox) ## 输出结果
  HR <- round(GSum$coefficients[,2],2) ## 输出HR值
  PValue <- round(GSum$coefficients[,5],3) ## 输出P值
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") ## 输出HR的执行区间
  HR.95L=round(GSum$conf.int[,"lower .95"],2) #分别输出HR的执行区间
  HR.95H=round(GSum$conf.int[,'upper .95'],2)
  Unicox <- data.frame("characteristics" = x, ## 返回结果，并构建数据框
                       "Hazard Ratio" = HR,
                       "CI95" = CI,
                       "HR.95L"=HR.95L,
                       "HR.95H"=HR.95H,
                       "P Value" = PValue)
  return(Unicox)
}

##赋予单因素COX回归函数数据集
data<-KIRC_lasso

VarNames <- colnames(KIRC_lasso)[12:ncol(KIRC_lasso)] ## 输出需要分析的变量名字
UniVar <- lapply(VarNames,UniCox) ## 批量做Cox分析
UniVar <- ldply(UniVar,data.frame) ## 将结果整理为一个数据框
(GetFactors <- UniVar$characteristics[which(UniVar$P.Value < 0.05)] %>% as.character()) ## 筛选其中P值<0.2的变量纳入多因素cox分析。
out<-rbind(gene=colnames(UniVar),UniVar)
write.table(out,"./data/Univarcox.txt",sep = "\t",row.names = F,col.names = F)

##正则化回归筛选变量 ##------------------------------------------------------------------

## 加载R包
library("glmnet") 
library(tidyverse)
library(rms)
library(Hmisc)
library(lattice)
library(Formula)
library(ggplot2)
library(foreign)
library(survMisc)
library(timeROC)
library(survminer)
library(ResourceSelection)
library(caret)
library(export)
## 根据R包的要求，将数据需要筛选的部分提取转换为矩阵,lasso不需要生存时间和生存状态，所以多删除两列
expr<- KIRC_lasso[,GetFactors] %>% as.matrix() 

##正则化lasso

set.seed(123)
cvfit = cv.glmnet(expr,
                  Surv(KIRC_lasso$futime,KIRC_lasso$fustat), 
                  nfold=10,#10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提
                  family = "cox",alpha=1)  #alpha一般默认为1

## 画图
par(lwd=2,cex=1.2)
plot(cvfit) 
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
graph2pdf(file="./Figure/lasso1.pdf", width = 8,height = 5)  
dev.off()

fit <- glmnet(expr, Surv(KIRC_lasso$futime,KIRC_lasso$fustat), family = "cox") 

par(lwd=2,cex=1.2)
plot(fit, xvar="lambda",label=TRUE)
graph2pdf(file="./Figure/lasso2.pdf", width = 8,height = 5)  
dev.off()

coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
active.min = which(coef.min != 0 ) ## 找出那些回归系数没有被惩罚为0的
(lassoGene <- colnames(expr)[active.min]) ## 提取基因名称

##个性化画图-------------------------------------------

cvfit$lambda.min
##提取绘图数据
xx <- data.frame(lambda=cvfit[["lambda"]],cvm=cvfit[["cvm"]],cvsd=cvfit[["cvsd"]],
                 cvup=cvfit[["cvup"]],cvlo=cvfit[["cvlo"]],nozezo=cvfit[["nzero"]])
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo,' vars')
##lasso1图
ggplot(xx,aes(ll,cvm,color=NZERO))+
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),width=0.05,size=1)+
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_point(size=2)+
  xlab("Log Lambda")+ylab('Partial Likelihood Deviance')+
  theme_bw(base_rect_size = 1.5)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_jco()(8)))+
  scale_x_continuous(expand = c(0.02,0.02))+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'bottom')+
  annotate('text',x = -5.5,y=12.9,label='Optimal Lambda = 0.0130',color='black')+
  guides(col=guide_legend(ncol = 3))
graph2pdf(file="./Figure/lasso1.2.pdf", width = 15,height = 10)  
dev.off()

##提取绘图数据
x <- coef(fit)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1] # extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm  

##lasso2图
ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Lambda (log scale)") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  #annotate('text',x = -3.3,y=0.26,label='Optimal Lambda = 0.012',color='black')+
  guides(col=guide_legend(ncol = 1))
graph2pdf(file="./Figure/lasso2.2.pdf", width = 15,height = 10)  
dev.off()

##多因素COX回归分析------------------------------------------------------------------------
library(survival)

rt<-data
fml <-  as.formula(paste0('BaSurv~',paste0(lassoGene,collapse = '+'))) ## 通过多因素Cox回归构建一个评分
MultiCox <- coxph(fml, data = data)
MultiCox=step(MultiCox,direction = "both")
MultiCoxSum <- summary(MultiCox)

outTab=data.frame()
outTab=cbind(
  coef=MultiCoxSum$coefficients[,"coef"],
  HR=MultiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=MultiCoxSum$conf.int[,"lower .95"],
  HR.95H=MultiCoxSum$conf.int[,"upper .95"],
  pvalue=MultiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="./data/multiCox.txt",sep="\t",row.names=F,quote=F)

##计算患者的风险评分---------------------------------------------------------------------------
library("survminer")

riskScore=predict(MultiCox,type="risk",newdata=data)
coxGene=rownames(MultiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("submitter_id","patient_id","futime","fustat",coxGene)

##使用survminer包确定最佳截断值，将患者分为高低风险组

data1<-cbind(data[,outCol],riskScore)
cut <- surv_cutpoint(data1,'futime','fustat','riskScore',minprop = 0.1)
cat <- surv_categorize(cut)
risk_level<-cat[,'riskScore']
out<-cbind(data1,risk_level) 
out<-rbind(id=colnames(out),out)
write.table(out,"./data/risk.txt",sep="\t",quote=F,col.names = F,row.names = F)

##K-M生存曲线 --------------------------------------------------------------------------
library(survival)
library("survminer")
library(ggsci)
library(export)

rt=read.table("./data/risk.txt",header=T,sep="\t",row.names = 1)

kmp <- function(data,legend,main){
  mytheme <- theme_survminer(font.legend = c(12,"plain", "black"),
                             font.x = c(12,"plain", "black"),
                             font.y = c(12,"plain", "black"),
                             legend = "top")
  fit <- survfit(Surv(futime,fustat)~risk_level,data)
  pp <- ggsurvplot(fit,data = data,
                   palette= pal_nejm()(2),
                   conf.int=FALSE,size=1.3,
                   pval=T,pval.method = T,
                   legend.labs=c('High-risk','Low-risk'), 
                   legend.title="",
                   legend=legend,
                   xlab="Time(months)",
                   ylab='Overall survival',
                   ggtheme = mytheme,risk.table = TRUE)
  return(pp)
}

kmp(data=rt,legend = c(0.84,0.97),main = NULL)
graph2pdf(file='./Figure/KM_survival.pdf',width=6,height=6)
dev.off()

##RiskPlot---------------------------------------------------------------------
library(pheatmap)
library(ggsci)
rt=read.table("./data/risk.txt",sep="\t",header=T,check.names=F,row.names = 1)
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk_level"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10

##RiskScore
pdf(file="./Figure/RiskScore.pdf",width = 12,height = 5)
par(lwd=2,cex=1.2)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#0072B599",lowLength),
           rep("#BC3C2999",highLength)))
#abline(h=median(rt$riskScore),v=lowLength,lty=2)
abline(h=max(rt[rt$risk_level=="low",]$riskScore),v=lowLength,lty=2)
dev.off()

#SurvStat
color=as.vector(rt$fustat)
color[color==1]="#BC3C2999"
color[color==0]="#0072B599"
pdf(file="./Figure/SurvStat.pdf",width = 12,height = 5)
par(lwd=2,cex=1.2)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (months)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#Heatmap
rt1=rt[c(4:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
annotation$type<- factor(annotation$type,levels =c("low","high"),labels = c("1","2"))
annotation$type<-factor(annotation$type,levels =c("1","2"),labels = c("low","high"))

##画热图
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
rt2=rt[c(4:(ncol(rt)-2))]
rt2<-t(scale(rt2)) %>% as.matrix()

annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)

identical(rownames(annotation),colnames(rt2))
table(annotation$type)

Cluster <- c('#eeba4d','#21b6af')
names(Cluster) <- levels(as.factor(annotation$type))

Top = HeatmapAnnotation(Risk=annotation$type,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                     title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col= list(Risk=Cluster),
                        show_annotation_name = F,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))


Heatmap(rt2,name='Z-score',
        top_annotation = Top,
        cluster_rows = FALSE,
        col=colorRamp2(c(-2,0,2),c("#0072B599","white","#BC3C2999")),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,299),rep(2,204)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), border = T,
                                  title_gp = gpar(fontsize = 10, fontface = "bold"))) 
library(export)
graph2pdf(file = "./Figure/Heatmap_risk2.pdf",width = 12,height = 5)

##riskScore的ROC曲线----------------------------------------------------------------------------------
library(survival)
library(survminer)
library(timeROC)
library(ggplot2)
library(ggsci)
library(tidyverse)

rt=read.table("./data/risk.txt",header=T,sep="\t",check.names=F,row.names = 1)
risk=rt[,c("futime", "fustat", "riskScore")]

##画时间依赖性ROC曲线-------------------------------------------------------
ROC_rt <- timeROC(T=risk$futime,delta=risk$fustat,marker=risk$riskScore,
                  cause = 1,weighting = 'marginal',times = c(12,36,60),ROC = TRUE)

TP <- ROC_rt$TP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'TP')
FP <- ROC_rt$FP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'FP')

dd <- TP
dd$FP <- FP$FP
dd$time <- ifelse(dd$time =='t=12',paste0("AUC at 1 years = ", sprintf("%.3f", ROC_rt$AUC[[1]])),
                  ifelse(dd$time=='t=36',paste0("AUC at 3 years = ", sprintf("%.3f", ROC_rt$AUC[[2]])),
                         paste0("AUC at 5 years = ", sprintf("%.3f", ROC_rt$AUC[[3]]))))
pdf(file="./Figure/ROC.pdf", width=5, height=5)
ggplot(dd,aes(FP,TP,color=time))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_nejm()+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))
dev.off()

##C-index-------------------------------------------------------------------------------------------
library(survcomp)

rt=read.table("./data/risk.txt",header=T,sep="\t",check.names=F,row.names = 1)
cindex <- concordance.index(x=rt$riskScore,
                            surv.time = rt$futime, 
                            surv.event = rt$fustat,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

##基于纳入临床指标和风险分组的COX回归分析-----------------------------------------------------
library(rms)
library(foreign)
library(survival)
library(survminer)

rt<-read.table("cox_input.txt",header=T,sep="\t")

rt$grade<-as.factor(rt$grade)
rt$T_stage<- as.factor(rt$T_stage)
rt$M_stage<-as.factor(rt$M_stage)
rt$gender<-as.factor(rt$gender)
rt$tumor_stage<-as.factor(rt$tumor_stage)

##下面两步目的是换参考组别
rt$risk <- factor(rt$risk,levels = c("low","high"),labels = c(1,2))
rt$risk <- factor(rt$risk,levels = c("1","2"),labels = c("low","high"))

rt$age<- factor(rt$age,levels = c("Yonger","Older"),labels = c(1,2))
rt$age <- factor(rt$age,levels = c("1","2"),labels = c("Yonger","Older"))

#step1，筛选P <0.1者
f <- as.formula(Surv(futime, fustat) ~gender+age+grade+T_stage+tumor_stage+risk)
cox <- coxph(f,data=rt)
summary(cox)

#ggforest(cox)
##写出cox结果
MultiCoxSum <- summary(cox)
outTab=data.frame()
outTab=cbind(
  coef=MultiCoxSum$coefficients[,"coef"],
  HR=MultiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=MultiCoxSum$conf.int[,"lower .95"],
  HR.95H=MultiCoxSum$conf.int[,"upper .95"],
  pvalue=MultiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="./data/clinical_multiCox.txt",sep="\t",row.names=F,quote=F)

#step2，筛选P<0.05者
f1 <- as.formula(Surv(futime, fustat) ~ age+tumor_stage+risk)
cox1 <- coxph(f1,data=rt)
summary(cox1)
#ggforest(cox1)

##Nomogram列线图--------------------------------------------------------------------------------
library(rms)
library(foreign)
library(survival)
library(regplot)
library(mstate)
library(export)

rt<-read.table("cox_input.txt",header=T,sep="\t")

rt$grade<-as.factor(rt$grade)
rt$T_stage<- as.factor(rt$T_stage)
rt$M_stage<-as.factor(rt$M_stage)
rt$gender<-as.factor(rt$gender)
rt$tumor_stage<-as.factor(rt$tumor_stage)

##下面两步目的是换参考组别
rt$risk <- factor(rt$risk,levels = c("low","high"),labels = c(1,2))
rt$risk <- factor(rt$risk,levels = c("1","2"),labels = c("low","high"))

rt$age<- factor(rt$age,levels = c("Yonger","Older"),labels = c(1,2))
rt$age <- factor(rt$age,levels = c("1","2"),labels = c("Yonger","Older"))

ddist <- datadist(rt)
options(datadist='ddist')

##Nomogram
label(rt$tumor_stage)<-"tumor stage"
cox <- cph(Surv(futime, fustat) ~ age+tumor_stage+risk,surv=T,x=T, y=T,data=rt) 
surv <- Survival(cox)
sur_1_year<-function(x)surv(1*12,lp=x)
sur_3_year<-function(x)surv(1*12*3,lp=x)
sur_5_year<-function(x)surv(1*12*5,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,funlabel=c('1-Year survival','3-Year survival','5-Year survival'),maxscale=100,fun.at= c('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'))
pdf("./Figure/nomogram.pdf")
par(lwd=2,cex=1.2)
plot(nom_sur,xfrac=0.4,col.grid = gray(c(0.8,0.95)))
dev.off()

##计算Nomogram模型的risk score，并将患者分为高低危组----------------------------------------------------------
library(survivalROC)
library(survival)
library(ggsci)
library(timeROC)
library(tidyverse)

rt<-read.table("roc_input.txt",header=T,sep="\t",check.names=F,row.names=1)

cox_m <- coxph(Surv(futime, fustat) ~ age+tumor_stage+risk_score, data = rt)
cox_m1<-step(cox_m,direction = "both")
risk_score<-predict(cox_m1,type="risk",newdata=rt)
risk_level<-as.vector(ifelse(risk_score>median(risk_score),"High","Low"))
write.table(cbind(id=rownames(cbind(rt[,1:3],risk_score,risk_level)),cbind(rt[,1:3],risk_score,risk_level)),
            "./data/risk_score.txt",sep="\t",quote=F,row.names=F)

##绘制ROC曲线
rt=read.table("./data/risk_score.txt",header=T,sep="\t",check.names=F,row.names=1)
risk=rt[,c("futime", "fustat", "risk_score")]
#1年ROC

ROC_rt <- timeROC(T=risk$futime,delta=risk$fustat,marker=risk$risk_score,
                  cause = 1,weighting = 'marginal',times = c(12,36,60),ROC = TRUE)

TP <- ROC_rt$TP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'TP')
FP <- ROC_rt$FP%>%as.data.frame()%>%pivot_longer(cols = 1:3,names_to = 'time',values_to = 'FP')

dd <- TP
dd$FP <- FP$FP
dd$time <- ifelse(dd$time =='t=12',paste0("AUC at 1 years = ", sprintf("%.3f", ROC_rt$AUC[[1]])),
                  ifelse(dd$time=='t=36',paste0("AUC at 3 years = ", sprintf("%.3f", ROC_rt$AUC[[2]])),
                         paste0("AUC at 5 years = ", sprintf("%.3f", ROC_rt$AUC[[3]]))))
pdf(file="./Figure/nomogram_ROC.pdf", width=5, height=5)
ggplot(dd,aes(FP,TP,color=time))+
  geom_line(size=1)+
  labs(x='1-Specificity',y='Sensitivity',color=NULL)+
  theme_bw(base_rect_size = 1.5)+
  geom_abline(slope = 1,color='grey70')+
  theme(panel.grid =element_blank(),
        axis.text = element_text(size=11),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = c(0.995,0.012),
        legend.justification = c(1,0))+
  scale_color_nejm()+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))
dev.off()

##C-index----------------------------------------------------------------------------------------
library(survcomp)

rt=read.table("./data/risk_score.txt",header=T,sep="\t",check.names=F,row.names=1)
cindex <- concordance.index(x=rt$risk_score,
                            surv.time = rt$futime, 
                            surv.event = rt$fustat,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

##Calibration-------------------------------------------------------------------------------------
library(rms)
library(foreign)
library(survival)
library(ggsci)

rt<-read.table("cox_input.txt",header=T,sep="\t")

rt$grade<-as.factor(rt$grade)
rt$T_stage<- as.factor(rt$T_stage)
rt$M_stage<-as.factor(rt$M_stage)
rt$gender<-as.factor(rt$gender)
rt$tumor_stage<-as.factor(rt$tumor_stage)
rt$age <- as.factor(rt$age)
rt$risk<-as.factor(rt$risk)

ddist <- datadist(rt)
options(datadist='ddist')

##合并画图
library(ggsci)

cox1 <- cph(Surv(futime, fustat) ~ age+tumor_stage+risk,surv=T,x=T, y=T,time.inc = 1*12,data=rt) 
cal1 <- calibrate(cox1, cmethod="KM", method="boot", u=1*12,m= nrow(rt)/3,B=1000)
cox3 <- cph(Surv(futime, fustat) ~ age+tumor_stage+risk,surv=T,x=T, y=T,time.inc = 3*12,data=rt) 
cal3 <- calibrate(cox3, cmethod="KM", method="boot", u=3*12, m=nrow(rt)/3, B=1000)
cox5 <- cph(Surv(futime, fustat) ~ age+tumor_stage+risk,surv=T,x=T, y=T,time.inc = 5*12,data=rt) 
cal5 <- calibrate(cox5, cmethod="KM", method="boot", u=5*12, m=nrow(rt)/3, B=1000)

pdf('./Figure/calibration.pdf',width = 8,height = 8, onefile=FALSE)
par(lwd=2,cex=1.2)
plot(cal1[,c('mean.predicted',"KM")],lwd = 2,lty = 1,xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted OS",
     ylab="Actual OS(proportion)",
     type = "b",col =pal_nejm()(3)[1],pch=16) 
abline(0,1,lty = 3, lwd = 2,col = "grey") 
lines(cal3[,c('mean.predicted',"KM")],  type = 'b',lwd = 2,pch = 16,col =pal_nejm()(3)[2]) 
lines(cal5[,c('mean.predicted',"KM")],  type = 'b',lwd = 2,pch = 16,col =pal_nejm()(3)[3]) 
legend("topleft",legend = c("1-year","3-year","5-year"),pch = c(16,16,16),lty = c(1,1,1),col = pal_nejm()(3))
dev.off()

##DCA-------------------------------------------------------------------------------------------------
library(Hmisc)
library(grid)
library(lattice)
library(Formula)
library(ggplot2)
library(rms)
library(survival)
library(rmda)
library(ggsci)

rt<- read.table('DCA_input.txt')

age<- decision_curve(Y~age,data = rt, family = binomial(link ='logit'),
                     thresholds= seq(0,1, by = 0.01),
                     confidence.intervals =0.95,study.design = 'cohort') 
tumor_stage<- decision_curve(Y~tumor_stage,data = rt, family = binomial(link ='logit'),
                         thresholds= seq(0,1, by = 0.01),
                         confidence.intervals =0.95,study.design = 'cohort')
risk_score<- decision_curve(Y~risk_score,data = rt, family = binomial(link ='logit'),
                            thresholds= seq(0,1, by = 0.01),
                            confidence.intervals =0.95,study.design = 'cohort')
nomogram<- decision_curve(Y~nomogram,data = rt, family = binomial(link ='logit'),
                          thresholds= seq(0,1, by = 0.01),
                          confidence.intervals =0.95,study.design = 'cohort')

List<- list(age,tumor_stage,risk_score,nomogram)
pdf(file = "./Figure/DCA.pdf",width = 8,height = 8,)
plot_decision_curve(List,curve.names= c('age','tumor_stage','risk_score','nomogram'),
                    
                    cost.benefit.axis =FALSE,col =pal_nejm()(4),
                    
                    confidence.intervals =FALSE,standardize = FALSE)
dev.off()

