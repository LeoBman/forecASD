#### load packages ####
require(Matrix)
require(RBGL)
require(snow)
require(randomForest)

#### STRING NETWORK (shortest paths matrix) ####
## "9606.protein.links.v10.txt.gz" is already in repo
# download.file(url = "http://version10.string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz", 
#               destfile = "./ext_data/9606.protein.links.v10.txt.gz")

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

### BRAINSPAN ####
## "genes_matrix_csv.zip" already in repo
## download.file(url = "http://www.brainspan.org/api/v2/well_known_file_download/267666529",
##               destfile = "./ext_data/genes_matrix_csv.zip") 

unzip("./ext_data/genes_matrix_csv.zip", exdir = "ext_data/brainspan/")

## Import brainspan data
m = read.table("./ext_data/brainspan/expression_matrix.csv",
               sep=",",header=F,stringsAsFactors=F)
ann = read.table("./ext_data/brainspan/rows_metadata.csv",
                 sep=",",header=T,stringsAsFactors=F)
fac = read.table("./ext_data/brainspan/columns_metadata.csv",
                 sep=",",header=T,stringsAsFactors=F)
M = as.matrix(m[,-1])

## try to map missing entrez gene IDs via gene symbols
eg.map = read.table("./ext_data/entrezgene2symbol.csv",
                    sep=",", header = T, stringsAsFactors = F)
eg.map = structure(eg.map$symbol, names = eg.map$entrez)
ann[is.na(ann$entrez_id),"entrez_id"] = names(eg.map)[match(ann[is.na(ann$entrez_id),]$gene_symbol,eg.map)]

rownames(M) = ann$entrez

## convert all ages into weeks
age = fac$age
age.tmp = strsplit(age,split=" ")
age.wk = sapply(age.tmp,function(x) {
  if(x[2]=="pcw") out=as.numeric(x[1])*1
  if(x[2]=="mos") out=as.numeric(x[1])*4.33 + 38
  if(x[2]=="yrs") out=as.numeric(x[1])*52 + 38
  return(out)
}
)

## keep only structures found in most samples
str = fac$structure_acronym
these.str = names(table(str)[table(str)>20])

## smooth and interpolate data ##
#### Lower number of parallel cores if fewer than 4 are available
cl = makeSOCKcluster(4)

## function to smooth and interpolate data ##
lw.smooth = function(y,x_age){
  lw = lowess(log(x_age),y,f=1/3)
  apx = approx(lw ,
               xout = seq(2, log(2118), length.out = 50),
               rule = 2)
  return(apx$y)
}

expr = list()

## expr is a list where each element corresponds to the expression in 
## one brain region. genes are columns, timepoints are rows
for(i in 1:length(these.str)){
  here = these.str[i]
  X = M[,str==here]
  x_age = age.wk[str==here]
  expr[[i]] = parApply(cl,X,1,lw.smooth,x_age)
  names(expr)[i] = these.str[i]
  print(i)
}

## transpose expression matricies so genes are rows
expr = lapply(expr,t)

## rbind rows of the same gene together, so that there is one
## matix per gene where rows are regions and columns are timepoints
expr = parLapply(cl, 1:nrow(expr[[1]]), function(i, y){
  do.call('rbind', lapply(y, function(x) {
    x[i,]
  }))
}, expr)

## relabel matricies by entrezID
names(expr) = rownames(M)
rn = rownames(expr[[1]])

## function to scale 
sc = function(x){
  matrix(
    scale(as.numeric(x)),
    nrow = 16,
    ncol = 50,
    dimnames = list( rownames(x) , colnames(x) )
  )
  }

expr.sc = lapply(expr,sc)


## entrezID to ensembl protein ID mapping
e2e = read.table("./ext_data/entrez_gene_id.vs.string.v10.28042015.tsv",sep="\t",stringsAsFactors=F)
## fix protein names
e2e[,2] = gsub("9606.","",e2e[,2],fixed=T)
e2e = structure(e2e[,2],names=as.character(e2e[,1]))

## vectorize the gene-wise expression matricies so
## rows are genes
bs = sapply(expr.sc,as.numeric)
bs = t(bs)
bs = bs[!is.na(rownames(bs)),]

## some IDs are missing from brainspan data, try to map them
missing = read.table("./ext_data/brainspan_missing_ids.txt",sep="\t",header=T,stringsAsFactors=F)
missing = missing[missing$Protein.stable.ID%in%rownames(string.path),]
missing = missing[!duplicated(missing$NCBI.gene.ID),]
missing = structure(missing$Protein.stable.ID,names=as.character(missing$NCBI.gene.ID))
e2e = c(e2e,missing)
rownames(bs) = e2e[rownames(bs)]
colnames(bs) = 1:800

## replace NAs with median
bs = na.roughfix(bs)
save(bs, file="01_BrainSpan_matrix.Rdata")

###### Training Labels ####

## read in SFARI genes and gene IDs
sfari = read.table("./ext_data/SFARI-Gene_genes_export01-11-2017.csv",sep=",",header=T,stringsAsFactors=F)
sfari.id = read.table("./ext_data/sfari_gene_ids.txt",sep="\t",stringsAsFactors=F,header=T)

## top are SFARI 1 and 2 genes
top = sfari$gene.symbol[sfari$gene.score %in% c("1","2")]
top = sfari.id$Protein.stable.ID[sfari.id$Gene.name %in% top]
top = unique( top[ top %in% rownames(string.path) ] )
## pos are top genes also in STRING
pos = rownames(string.path)[ rownames(string.path) %in% top ]
pos = unique(pos)
pos = pos[ !is.na(pos) ]
rn = rownames(string.path)[ rownames(string.path) %in% rownames(bs) ]

## negative genes are random background genes
set.seed(3716359)
neg = sample(rn[ !(rn %in% sfari.id$Protein.stable.ID) ] , 1000)

save(pos,neg,file="01_training_labels.Rdata")
