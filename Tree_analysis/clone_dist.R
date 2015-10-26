#cells <- c("7PB", "14hi", "14lo", "90hi", "90lo")
library(ape)
library(dplyr)

#setup database 
my_db <- src_sqlite("~/Documents/RepSeq2/trees/clone_dist.sqlite", create = T)
dist <- c("NA", "NA", "NA", "NA", "NA")
names(dist) <- c("cell1", "cell2", "dist", "clone", "sample")
dist <- as.data.frame(t(dist))
copy_to(my_db, dist, temporary = FALSE)


#writes pairwise distances for sequences in a clone to database
get_dist <- function(f, db, sample){
  tree <- read.tree(f)
  x <- cophenetic.phylo(tree)
  x <- as.data.frame(x)

  #make all unique combinations of seqIDs
  combo <- t(combn(names(x),2))
  #pull distances for each combo from distance matrix
  dist <- x[combo]
  #make pairwise dataframe
  df <- data.frame(combo, dist)
  names(df) <- c("cell1", "cell2", "dist")
  #trim seqIDs to just cell type
  df$cell1 <- substring(df$cell1, 0, 3)
  df$cell2 <- substring(df$cell2, 0, 3)
  #add clonename
  clone <- unlist(strsplit(f, "[.]"))[2]
  df$clone <- clone
  #add sample name
  df$sample <- sample
  #write to database
  db_insert_into(con = my_db$con, table = "dist", values = df)
}

#get all bestTree filenames
setwd("~/Documents/RepSeq2/trees/bestTree_011/")
filenames <- dir(getwd(), "*.phy")
#apply get_dist function to all bestTree files 
lapply(filenames, get_dist, db=my_db, sample="011")

#get 007 files
setwd("~/Documents/RepSeq2/trees/bestTree_007/")
filenames <- dir(getwd(), "*.phy")
lapply(filenames, get_dist, db=my_db, sample="007")

#get 012 files
setwd("~/Documents/RepSeq2/trees/bestTree_012/")
filenames <- dir(getwd(), "*.phy")
lapply(filenames, get_dist, db=my_db, sample="012")


