#### load packages ####
require(Matrix)
require(RBGL)

#### STRING network shortest paths matrix ####
download.file(url = "http://version10.string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz", 
              destfile = "./ext_data/9606.protein.links.v10.txt.gz")

string.dat = read.table("./ext_data/9606.protein.links.v10.txt.gz", 
                 stringsAsFactors=F,
                 header=T)

## clean up protein IDs
string.dat[,1] = gsub("9606." , "", string.dat[,1], fixed=T)
string.dat[,2] = gsub("9606.", "", string.dat[,2], fixed=T)

## keep interactions with scores over 400
string.keep = string.dat[string.dat[,3]>400,]
string.graph = ftM2graphNEL(as.matrix(string.keep[,1:2]))
string.path = johnson.all.pairs.sp(string.graph)

save(string.graph, string.path, file="01_STRING_graph.Rdata")

### BrainSpan data ####
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

sfari = read.table("./ext_data/SFARI-Gene_genes_export01-11-2017.csv",sep=",",header=T,stringsAsFactors=F)
sid = read.table("./ext_data/sfari_gene_ids.txt",sep="\t",stringsAsFactors=F,header=T)

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