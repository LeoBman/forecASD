# # # # # # # # # # # # 
#### Load Packages ####
# # # # # # # # # # # #
require(randomForest)


# # # # # # # # # # # #
#### BrainSpan RF #####
# # # # # # # # # # # #
load("01_BrainSpan_matrix.Rdata")
load("01_training_labels.Rdata")

pobs.x = bs[ rownames(bs) %in% c(pos,neg) , ]
bs.y = as.factor(rownames(bs.x) %in% pos)

## train random forest
set.seed(5136)
bs.rf = randomForest(y=bs.y,
                     x=bs.x,
                     ntree=1000,
                     importance=T,
                     strata=bs.y,
                     sampsize=c(76,76),
                     do.trace=10,
                     proximity=T)

## get predictions on all remaining genes
bs.z =  na.roughfix(bs[ !rownames(bs) %in% rownames(bs.x) , ])

bs.prd = predict(bs.rf, bs.z, type="prob")
bs.prd = rbind(bs.rf$vote, bs.prd)

save(bs.rf, bs.prd, file="02_brainspan_RF.Rdata")

# # # # # # # # # # # # 
#### the STRING RF ####
# # # # # # # # # # # # 
load("01_STRING_graph.Rdata")

string.x = string.path[rownames(string.path) %in% c(pos,neg), ]
string.y = as.factor( rownames(string.x) %in% pos )

set.seed(2176)
string.rf = randomForest(
  y = string.y,
  x = string.x,
  importance = T,
  strata = string.y,
  sampsize = c(77,77),
  ntree=500, 
  do.trace=10,
  proximity=T)

## stratify by predicted class lables
strata = string.rf$y:string.rf$pred
## sample size is smallest class in strata
sampsize = rep(min(table(strata)),4)

## fit a model while omitting variables which are not used in any tree
set.seed(679)
while( sum( varUsed(string.rf) < 1 ) > 0 ){
 string.rf = randomForest(
   y = string.y,
   x = string.x[, rownames(string.rf$importance)[varUsed(string.rf) > 0] ],
   importance = T,
   strata = strata,
   sampsize = sampsize,
   ntree = 500,
   do.trace = 100,
   proximity = T)
 
 print(paste(sum(varUsed(string.rf)==0), "out of", nrow(string.rf$importance), "variables unused"))
}


string.prd = predict(string.rf, string.path[!rownames(string.path) %in% rownames(string.x), ], type="prob")
string.prd = rbind(string.prd, string.rf$vote)
string.prd = string.prd[ order(string.prd[,2], decreasing=T), ]

save(string.rf,string.prd,file="02_optimized_STRING_rf.Rdata")

# # # # # # # # # # # # # #
#### integrate scores #####
# # # # # # # # # # # # # #
## the STRING RF score, the BrainSpan RF score, and the TADA score

load("01_id_conversion.Rdata")

nn = rownames(prd)[rownames(prd)%in%rownames(bs.prd)]
load("/sdata/NCBI/entrezgene/entrezgene2symbol.Rdata")
symbol = eg_map[names(e2e)[match(nn,e2e)]]

tada = read.table("/wdata/jmichaelson/SPARK/WG_network/tada_BFs.txt",sep="\t",header=T,stringsAsFactors=F)
tada = structure(tada[[2]],names=tada[[1]])

Xpr = data.frame(string=prd[nn,2],brainspan=bs.prd[nn,2],tada=tada[symbol])
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
	BrainSpan_score=bs.prd[names(s),2],forecASD_no_TADA=s2,
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


