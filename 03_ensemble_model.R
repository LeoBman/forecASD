# # # # # # # # # # # # #
#### Ensembel Model #####
# # # # # # # # # # # # #

library(randomForest)

load("01_training_labels.Rdata")
load("02_STRING_rf.Rdata")
load("02_brainspan_RF.Rdata")
load("02_network_rf.Rdata")

# # # # # # # # # # # # 
#### Load in data #####
# # # # # # # # # # # # 

meta = read.csv("./ext_data/composite_table.csv", stringsAsFactors = F, row.names = 1)

string.prd = string.prd[rownames(meta), ]
bs.prd = bs.prd[rownames(meta), ]

## combine other predictors with network scores
meta = cbind(
  data.frame(
    STRING_score = string.prd[rownames(meta) , "TRUE"],
    BrainSpan_score = bs.prd[rownames(meta), "TRUE"]
  ),
  meta[rownames(meta), ]
)

# # # # # # # # # # # # 
#### train forest #####
# # # # # # # # # # # # 

meta.train = na.roughfix(
  meta[meta$ensembl_string %in% c(pos,neg), -(3:9)]   ## remove gene identifiers, etc.
) 
y = as.factor(rownames(meta.train) %in% pos)

set.seed(43775)
rf = randomForest(
  y = y,
  x = meta.train, 
  importance = T,
  do.trace = 10,
  strata = y,
  sampsize = c(76,76)
  )

meta.test = na.roughfix(
  meta[!(meta$ensembl_string %in% c(pos,neg)), -(3:9)]   ## remove gene identifiers, etc.
) 

meta.prd <- predict(rf, 
        meta.test, 
        type = "prob")

meta.score <- rbind(rf$votes, meta.prd)

final.data <- cbind(
  data.frame(
    forecASD = meta.score[rownames(meta),"TRUE"],
    'STRING+BrainSpan_RF' = network.prd[rownames(meta),"TRUE"]
  ), 
  meta
)

write.csv(final.data, 
          file = "forecASD_table.csv",
          quote = F, 
          row.names = F)

