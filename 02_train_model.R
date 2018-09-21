# # # # # # # # # # # # 
#### Load Packages ####
# # # # # # # # # # # #
require(randomForest)

# # # # # # # # # # # #
#### BrainSpan RF #####
# # # # # # # # # # # #
load("01_BrainSpan_matrix.Rdata")
load("01_training_labels.Rdata")

bs.x = bs[ rownames(bs) %in% c(pos,neg) , ]
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
 
 print(paste(sum(varUsed(string.rf) == 0), "out of", nrow(string.rf$importance), "variables unused"))
}


string.prd = predict(string.rf, string.path[!rownames(string.path) %in% rownames(string.x), ], type="prob")
string.prd = rbind(string.prd, string.rf$vote)
string.prd = string.prd[ order(string.prd[,2], decreasing=T), ]

save(string.rf,string.prd,file="02_STRING_rf.Rdata")

# # # # # # # # # # # # # #
#### integrate scores #####
# # # # # # # # # # # # # #
## the STRING RF score, the BrainSpan RF score, and the TADA score

load("01_id_conversion.Rdata")

nn = rownames(string.prd)[rownames(string.prd)%in%rownames(bs.prd)]
symbol = eg.map[names(e2e)[match(nn,e2e)]]



x.pred = data.frame(
  string = string.prd[nn, "TRUE"],
  brainspan = bs.prd[nn, "TRUE"]
  )
x.pred = na.roughfix(x.pred)

x.train = x.pred[rownames(x.pred) %in% c(pos,neg), ]
y.train = as.factor(rownames(x.train) %in% pos)


set.seed(5393)
network.rf = 
  randomForest(
    y = y.train,
    x = x.train[,c("string", "brainspan")],
    importance = T,
    strata = y.train,
    sampsize = c(76,76),
    do.trace = 50,
    ntree = 1000)

network.prd = predict(network.rf,
                          x.pred[!rownames(x.pred) %in% rownames(network.rf$votes), c("string", "brainspan")],
                          type=  "prob")
network.prd = rbind(network.prd, network.rf$votes)

save(network.rf, network.prd, file="02_network_rf.Rdata")
