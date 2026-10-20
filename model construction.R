
##单因素cox回归--------
outTab=data.frame()
sigGenes=c("MergedColumn","vital_status.demographic")
outTab <- data.frame(Gene=character(), HR=numeric(), HR.95L=numeric(), HR.95H=numeric(), pvalue=numeric(), stringsAsFactors=FALSE)
# 遍历列
for(i in colnames(exp_merge[,3:ncol(exp_merge)])){
  tryCatch({
    cox <- coxph(Surv(MergedColumn, vital_status.demographic) ~ exp_merge[,i], data = exp_merge)
    coxSummary <- summary(cox)
    coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
    
    # 检查coxP是否是缺失值
    if(!is.na(coxP) && coxP < pFilter){
      sigGenes <- c(sigGenes, i)
      newRow <- data.frame(
        id=i,
        HR=coxSummary$conf.int[,"exp(coef)"],  # 注意：通常这里应该是HR的直接估计值，而不是"exp(coef)"的列名，除非您确实有这样的列名
        HR.95L=coxSummary$conf.int[,"lower .95"],
        HR.95H=coxSummary$conf.int[,"upper .95"],
        pvalue=coxP,
        stringsAsFactors=FALSE
      )
      outTab <- rbind(outTab, newRow)
    }
  }, error=function(e) {
    # 这里可以添加错误处理逻辑，例如打印错误信息或跳过当前迭代
    warning(paste("Error fitting Cox model for column:", i, "Error message:", e$message))
  })
}
#####lasso缩小基因范围#####
#install.packages("glmnet")
library(glmnet)
#outTab中有87个基因（单因素cox）
write.table(outTab,file="UniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=exp_merge[,sigGenes]
uniSigExp_data=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp_data,file="UniSigExp.txt",sep="\t",row.names=F,quote=F)
exptimefile="UniSigExp.txt"
unifile="UniCox.txt"

yifeng<-read.table(exptimefile,header = T,sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F)
biodata=read.table(unifile,header = T,sep = "\t",check.names = F)
yifeng$MergedColumn=yifeng$MergedColumn/365
yifeng$MergedColumn=as.numeric(yifeng$MergedColumn)
yifeng$vital_status.demographic=as.numeric(yifeng$vital_status.demographic)
#yifeng=yifeng[which(yifeng$MergedColumn !=0)]
v1=as.matrix(yifeng[,c(3:ncol(yifeng))])
v2=as.matrix(Surv(yifeng$MergedColumn,yifeng$vital_status.demographic))#患者的生存时间和状态

myfit<-glmnet(v1,v2,family = "cox",alpha = 1)
pdf("lambda08.pdf")
plot(myfit, xvar = "lambda", label = TRUE)
dev.off()

myfit2<-cv.glmnet(v1, v2,  family="cox",alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit08.pdf",width=6,height=5.5)
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.lse)),lty="dashed")
dev.off()

coe <- coef(myfit, s=myfit2$lambda.min)
ant_index<-which(coe !=0)
act_coe<-coe[ant_index]
lassoGene=row.names(coe)[ant_index]

prognostic_genes_df <- data.frame(
  Gene = lassoGene,
  Coefficient = act_coe
)

rownames(outTab)=outTab[[1]]#设置行名
#####绘制单因素Cox森林图#####
rt=outTab[lassoGene,]
gene <- rownames(rt)
hr <- sprintf("%.3f",as.numeric(rt$"HR"))
hrLow  <- sprintf("%.3f",as.numeric(rt$"HR.95L"))
hrHigh <- sprintf("%.3f",as.numeric(rt$"HR.95H"))
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(as.numeric(rt$pvalue)<0.001, "<0.001", sprintf("%.3f", as.numeric(rt$pvalue)))
#输出图形
pdf(file="unforest.pdf", width = 6,height = 6)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))


xlim = c((min(as.numeric(hrLow))-0.1),(max(as.numeric(hrHigh))+0.1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#####多因素Cox森林图####
rt=data.frame(uniSigExp[,1:2],uniSigExp[,lassoGene])
#rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)#对表达量取log
cox <- coxph(Surv(MergedColumn,vital_status.demographic) ~ ., data = rt)
cox=step(cox,direction = "both")
riskScore=predict(cox,type="risk",newdata=rt)#得到风险值
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))#分组，中位数分，高低风险组
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="risk.txt",sep="\t",quote=F,row.names=F)


#####绘制多因素Cox森林图#
coxSummary = summary(cox)
outTab=cbind(HR=coxSummary$conf.int[,"exp(coef)"],
             HR.95L=coxSummary$conf.int[,"lower .95"],
             HR.95H=coxSummary$conf.int[,"upper .95"],
             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
risk=data.frame(rt[,1:2],rt[,rownames(outTab)],riskScore,risk)
rt=as.data.frame(outTab)
gene <- rownames(rt)
hr <- sprintf("%.3f",as.numeric(rt$"HR"))
hrLow  <- sprintf("%.3f",as.numeric(rt$"HR.95L"))
hrHigh <- sprintf("%.3f",as.numeric(rt$"HR.95H"))
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(as.numeric(rt$pvalue)<0.001, "<0.001", sprintf("%.3f", as.numeric(rt$pvalue)))
#输出图形
pdf(file="mulforest.pdf", width = 6,height = 6)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c((min(as.numeric(hrLow))-0.1),(max(as.numeric(hrHigh))+0.1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#####患者分布可视化（风险曲线，生存图，风险热图）#####
rt=risk   
rt=rt[order(rt$riskScore),]                                            
#绘制风险曲线
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore.pdf",width = 10,height = 3.5)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("green",lowLength),rep("red",highLength)) )
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
dev.off()

#绘制生存图
color=as.vector(rt$vital_status.demographic)
color[color==1]="red"
color[color==0]="green"
pdf(file="survstat.pdf",width = 10,height = 3.5)
plot(rt$MergedColumn, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#绘制风险热图
rt1=rt[c(3:(ncol(rt)-2))]
rt1=log2(rt1+1)
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="risk_heatmap.pdf",width = 10,height = 10)
p=pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           fontsize_row=11,
           show_colnames = F,
           fontsize_col=3,
           scale = "row",
           color = colorRampPalette(c("green", "black", "red"))(100) )
dev.off()




