### get the STRING shortest paths matrix
library(Matrix)
library(RBGL)
dat = read.table("/sdata/STRING/9606.protein.links.v10.txt.gz",stringsAsFactors=F,header=T)
dat[,1] = gsub("9606.","",dat[,1],fixed=T)
dat[,2] = gsub("9606.","",dat[,2],fixed=T)

d2 = dat[dat[,3]>400,]
dim(d2)
g = ftM2graphNEL(as.matrix(d2[,1:2]))
sp = johnson.all.pairs.sp(g)

dim(sp)
summary(colSums(sp))
sp[1:5,1:5]

save(g,sp,file="STRING_graph.Rdata")
#savehistory("STRING_shortest_paths_and_graph.Rhistory") ## in /sdat/STRING

### BrainSpan data
load("/sdata/BRAINSPAN/genes_matrix_csv/lowess_smoothed_expression_scaled.Rdata")
e2e = read.table("/sdata/STRING/entrez_gene_id.vs.string.v10.28042015.tsv",sep="\t",stringsAsFactors=F)
e2e[,2] = gsub("9606.","",e2e[,2],fixed=T)
e2e = structure(e2e[,2],names=as.character(e2e[,1]))
bs = sapply(expr_smth_sc,as.numeric)
bs = t(bs)
bs = bs[!is.na(rownames(bs)),]
missing = read.table("brainspan_missing_ids.txt",sep="\t",header=T,stringsAsFactors=F)
missing = missing[missing[[4]]%in%rownames(sp),]
missing = missing[!duplicated(missing[[3]]),]
missing = structure(missing[[4]],names=as.character(missing[[3]]))
e2e = c(e2e,missing)
rownames(bs) = e2e[rownames(bs)]
colnames(bs) = 1:800
bs = na.roughfix(bs)
save(bs,file="brainspan_expr_matrix_ensembl.Rdata")

### the training labels

sfari = read.table('SFARI-Gene_genes_export01-11-2017.csv',sep=",",header=T,stringsAsFactors=F)
sid = read.table("/sdata/STRING/sfari_gene_ids.txt",sep="\t",stringsAsFactors=F,header=T)

top = sfari$gene.symbol[sfari$gene.score%in%c("1","2")]
top = sid[sid[[2]]%in%top,][[4]]
top = unique(top[top%in%rownames(sp)])

pos = rownames(sp)[rownames(sp)%in%top]
#pos = c(pos,as.vector(out$ensembl[log10(out$tada_bf)>2]))
pos = unique(pos)
pos = pos[!is.na(pos)]
#neg = sample(unique(as.vector(out$ensembl[log10(out$tada_bf) < -0.5])),1000)
#neg = neg[!neg%in%as.vector(sfari$ensembl)]
rn = rownames(sp)[rownames(sp)%in%rownames(bs)]
set.seed(3716359)
neg = sample(rn[!rn%in%sid[[4]]],1000)

save(pos,neg,file="training_labels.Rdata")



