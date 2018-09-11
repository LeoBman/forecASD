library(randomForest)
setwd('/wdata/jmichaelson/SPARK/WG_network/forecASD/')
##############################################
### the BrainSpan RF
##############################################
load("brainspan_expr_matrix_ensembl.Rdata")
load("training_labels.Rdata")

Xbs = bs[rownames(bs)%in%c(pos,neg),]
ybs = as.factor(rownames(Xbs)%in%pos)

set.seed(5136)
rfbs = randomForest(y=ybs,x=Xbs,ntree=1000,importance=T,strata=ybs,sampsize=c(76,76),do.trace=10,proximity=T)

prdBS = predict(rfbs,na.roughfix(bs[!rownames(bs)%in%rownames(Xbs),]),type="prob")
prdBS = rbind(rfbs$vote,prdBS)

save(rfbs,prdBS,file="brainspan_RF.Rdata")

##############################################
### the STRING RF
##############################################
load("STRING_graph.Rdata")
load("training_labels.Rdata")
X = sp[rownames(sp)%in%c(pos,neg),]
y = rownames(X)%in%pos
y = as.factor(rownames(X)%in%pos)

set.seed(2176)
rf = randomForest(y=y,x=X,importance=T,strata=y,sampsize=c(77,77),ntree=500,do.trace=10,proximity=T)

strata = rf$y:rf$pred
sampsize = rep(min(table(strata)),4)

set.seed(679)
while(sum(varUsed(rf)<1)>0){
 rf = randomForest(y=y,x=X[,rownames(rf$importance)[varUsed(rf)>0]],importance=T,strata=strata,sampsize=sampsize,ntree=500,do.trace=100,proximity=T)
 print(paste(sum(varUsed(rf)==0),nrow(rf$importance)))
}


prd = predict(rf,sp[!rownames(sp)%in%rownames(X),],type="prob")
prd = rbind(prd,rf$vote)
prd = prd[order(prd[,2],decreasing=T),]

save(rf,prd,file="optimized_feat_sel_STRING_RF.Rdata")


### integrating the STRING RF score, the BrainSpan RF score, and the TADA score
e2e = read.table("/sdata/STRING/entrez_gene_id.vs.string.v10.28042015.tsv",sep="\t",stringsAsFactors=F)
e2e[,2] = gsub("9606.","",e2e[,2],fixed=T)
e2e = structure(e2e[,2],names=as.character(e2e[,1]))
missing = read.table("../brainspan_missing_ids.txt",sep="\t",header=T,stringsAsFactors=F)
missing = missing[missing[[4]]%in%rownames(sp),]
missing = missing[!duplicated(missing[[3]]),]
missing = structure(missing[[4]],names=as.character(missing[[3]]))
e2e = c(e2e,missing)

nn = rownames(prd)[rownames(prd)%in%rownames(prdBS)]
load("/sdata/NCBI/entrezgene/entrezgene2symbol.Rdata")
symbol = eg_map[names(e2e)[match(nn,e2e)]]

tada = read.table("/wdata/jmichaelson/SPARK/WG_network/tada_BFs.txt",sep="\t",header=T,stringsAsFactors=F)
tada = structure(tada[[2]],names=tada[[1]])

Xpr = data.frame(string=prd[nn,2],brainspan=prdBS[nn,2],tada=tada[symbol])
Xpr = na.roughfix(Xpr)

Xtr = Xpr[rownames(Xpr)%in%c(pos,neg),]
ytr = as.factor(rownames(Xtr)%in%pos)

rfm = randomForest(y=ytr,x=Xtr,importance=T,strata=ytr,sampsize=c(76,76),do.trace=50,ntree=1000)

prdm = predict(rfm,Xpr[!rownames(Xpr)%in%rownames(rfm$votes),],type="prob")
prdm = rbind(prdm,rfm$votes)

save(rfm,prdm,file="forecASD_rf_pred.Rdata")


set.seed(5393)
rfm_notada = randomForest(y=ytr,x=Xtr[,1:2],importance=T,strata=ytr,sampsize=c(76,76),do.trace=50,ntree=1000)

prdm_notada = predict(rfm_notada,Xpr[!rownames(Xpr)%in%rownames(rfm$votes),1:2],type="prob")
prdm_notada = rbind(prdm_notada,rfm_notada$votes)

save(rfm_notada,prdm_notada,file="forecASD_rf_pred_noTADA.Rdata")

#############################################################################
### output table for distribution: ensembl, entrez, symbol, pLI, Iowa score, 
### Krishnan score, TADA bayes factor, SFARI gene score, #SPARK pilot DNMs
#############################################################################
load("/sdata/NCBI/entrezgene/entrezgene2symbol.Rdata")
load("/sdata/STRING/combined_score.Rdata")
pli = read.table("/sdata/ExAC/pLI/pLI.txt",sep="\t",header=T,stringsAsFactors=F)
pli = structure(pli$pLI,names=pli$gene)

sfari = read.table('SFARI-Gene_genes_export01-11-2017.csv',sep=",",header=T,stringsAsFactors=F)
sscore = structure(sfari$gene.score,names=sfari$gene.symbol)
tada = read.table("/wdata/jmichaelson/SPARK/WG_network/tada_BFs.txt",sep="\t",header=T,stringsAsFactors=F)
tada_mut = structure(tada[[3]],names=tada[[1]])
tada = structure(tada[[2]],names=tada[[1]])


krishnan = read.table("/wdata/jmichaelson/SPARK/WG_network/princeton_asd.txt",sep="\t",stringsAsFactors=F,header=T)
krishnan = structure(krishnan[[3]],names=krishnan[[2]])

s = structure(prdm[,2],names=rownames(prdm))
s2 = structure(prdm_notada[,2],names=rownames(prdm_notada))
df_out = data.frame(ensembl_string=names(s),entrez=names(e2e)[match(names(s),e2e)],
	symbol=eg_map[names(e2e)[match(names(s),e2e)]],forecASD=s,STRING_score=prd[names(s),2],
	BrainSpan_score=prdBS[names(s),2],forecASD_no_TADA=s2,
	stringsAsFactors=F)

df_out$mutation_rate = tada_mut[df_out$symbol]
df_out$SFARI_listed = df_out$symbol%in%sfari$gene.symbol
df_out$SFARI_score = sscore[df_out$symbol] 
df_out$pLI = pli[df_out$symbol]
df_out$TADA_BF = tada[df_out$symbol]
df_out$krishnan_post = krishnan[df_out$symbol]

save(df_out,file="/wdata/jmichaelson/SPARK/WG_network/forecASD/forecASD_table.Rdata")

write.table(df_out[order(df_out$forecASD,decreasing=T),],file="ASD_gene_prioritization_21Dec2017.txt",sep="\t",row.names=F)

### look at how the RF carves out the decision space in the 3d of the input features
p3d(df_out$STRING_score,df_out$BrainSpan_score,log10(df_out$TADA),col=ifelse(df_out$forecASD>0.45,'red','grey'))

flood = cbind(runif(6e4,0,1),runif(6e4,0,1),sample(df_out$TADA_BF,6e4,replace=T))
pflood = predict(rfm,na.roughfix(flood),type="prob")
#p3d(flood[pflood[,2]<0.45,],size=10,col='orangered')

ct = cut(flood[,3],breaks=quantile(df_out$TADA,c(0,0.2,0.4,0.6,0.8,1),na.rm=T))
par(mfrow=c(1,5))
for(i in 1:5){
 idx = as.integer(ct)==i
 plot(flood[idx,1:2],col=ifelse(pflood[idx,2]>0.45,'red','grey'),pch=16,cex=2.5,
	xlab="STRING score",ylab="BrainSpan score")
}


### training an ensemble classifier with forecASD(noTADA), Krishnan, 
### DAMAGES, and 2014/2015 TADA scores
load("forecASD_table_w_extra_stuff.Rdata")
tada = read.table("sanders_tada_2015.txt",sep="\t",header=T,stringsAsFactors=F)
tada = tada[,c(2,16,18,20,21)]
rownames(tada) = tada[,1]

load("training_labels.Rdata")
Xtr = foo[foo$ensembl%in%c(pos,neg),]
Xtr_tada = tada[Xtr$symbol,]

y = as.factor(Xtr$ensembl_string%in%pos)
Xtr = cbind(Xtr[,c(5,6,12,13,15,16)],Xtr_tada[,-1])

set.seed(43775)
Xtr = na.roughfix(Xtr)
rf = randomForest(y=y,x=Xtr,importance=T,do.trace=10,strata=y,sampsize=c(76,76))

Xtst = foo[!foo$ensembl%in%c(pos,neg),]
Xtst_tada = tada[Xtst$symbol,]

Xtst = cbind(Xtst[,c(5,6,12,13,15,16)],Xtst_tada[,-1])
prd = predict(rf,na.roughfix(Xtst),type="prob")

foo$meta = prd[foo$ensembl,2]
foo = cbind(foo,tada[foo$symbol,-1])

### the "meta" score has lower AIC and better deviance explained (1-(res/null))
### than any other TADA variable (when modeling DNMs occurring in non-SFARI genes)
summary(glm(Y[,1]~mutation_rate+meta,family="poisson",data=foo[!foo$SFARI_listed,]))
summary(glm(Y[,1]~mutation_rate+log10(tadaFdrAscSscExomeSscAgpSmallDel),family="poisson",data=foo[!foo$SFARI_listed,]))


######################################################################
### Feb 21, 2018 - add in DAWN metrics as features
######################################################################
load("forecASD_w_extra_plus_DAWN.Rdata")
library(randomForest)

load("forecASD_table_w_extra_stuff.Rdata")
tada = read.table("sanders_tada_2015.txt",sep="\t",header=T,stringsAsFactors=F)
tada = tada[,c(2,16,18,20,21)]
rownames(tada) = tada[,1]

load("training_labels.Rdata")
Xtr = foo[foo$ensembl%in%c(pos,neg),]
Xtr_tada = tada[Xtr$symbol,]

y = as.factor(Xtr$ensembl_string%in%pos)
Xtr = cbind(Xtr[,c(5,6,12,13,15,16)],Xtr_tada[,-1])
Xtr = cbind(Xtr,df[rownames(Xtr),23:25])

set.seed(43775)
Xtr = na.roughfix(Xtr)
rf = randomForest(y=y,x=Xtr,importance=T,do.trace=10,strata=y,sampsize=c(76,76))

Xtst = foo[!foo$ensembl%in%c(pos,neg),]
Xtst_tada = tada[Xtst$symbol,]

Xtst = cbind(Xtst[,c(5,6,12,13,15,16)],Xtst_tada[,-1])
Xtst = cbind(Xtst,df[rownames(Xtst),23:25])

prd = predict(rf,na.roughfix(Xtst),type="prob")
prd = rbind(rf$votes,prd)

df$meta2 = prd[df$ensembl,2]
save(df,file="forecASD_w_extra_plus_DAWN_v2.Rdata")


