#make sql table with only nececssary data

#submit IMGT 1, 2, 4 tables
#refer to RepSeq Data Dict for all data available

require(reshape2)

mk_big_df <- function(fn.IMGT_1, fn.IMGT_2, fn.IMGT_4, cell.type){
  IMGT_1 <- read.delim(fn.IMGT_1, header = TRUE)
  IMGT_2 <- read.delim(fn.IMGT_2, header = TRUE)
  IMGT_4 <- read.delim(fn.IMGT_4, header = TRUE)
  
  #select required columns
  IMGT_1 <- IMGT_1[,c(3,4,10, 14, 18, 21, 29)]
  IMGT_2 <- IMGT_2[,15]
  IMGT_4 <- IMGT_4[,15]
  
  #fix colnames
  df <- cbind(IMGT_1, IMGT_2, IMGT_4)
  colnames <- c("Functionality", "VGENE", "JGENE", "DGENE", "CDR3_len", "Junction", "Sequence", "CDR3_NT", "CDR3_AA")
  names(df) <- colnames
  
  #add cell type
  df$cell_type <- cell.type
  
  return(df)
}

#011
df.d7 <- mk_big_df("~/Documents/RepSeq/011_d7-PB/1_Summary_011_d7_070712.txt", 
          "~/Documents/RepSeq/011_d7-PB/2_IMGT-gapped-nt-sequences_011_d7_070712.txt", 
          "~/Documents/RepSeq/011_d7-PB/4_IMGT-gapped-AA-sequences_011_d7_070712.txt", 
          "d7.PB")

df.d14lo <- mk_big_df("~/Documents/RepSeq/011-D14-Low/1_Summary_011-D14A_190413.txt",
                      "~/Documents/RepSeq/011-D14-Low/2_IMGT-gapped-nt-sequences_011-D14A_190413.txt",
                      "~/Documents/RepSeq/011-D14-Low/4_IMGT-gapped-AA-sequences_011-D14A_190413.txt",
                      "d14.lo")

df.d14hi <- mk_big_df("~/Documents/RepSeq/011-D14-Hi/1_Summary_011-D14B_190413.txt",
                      "~/Documents/RepSeq/011-D14-Hi/2_IMGT-gapped-nt-sequences_011-D14B_190413.txt",
                      "~/Documents/RepSeq/011-D14-Hi/4_IMGT-gapped-AA-sequences_011-D14B_190413.txt",
                      "d14.hi")

df.d90lo <- mk_big_df("~/Documents/RepSeq/011-D90-Low/1_Summary_011-D90A_190413.txt",
                      "~/Documents/RepSeq/011-D90-Low/2_IMGT-gapped-nt-sequences_011-D90A_190413.txt",
                      "~/Documents/RepSeq/011-D90-Low/4_IMGT-gapped-AA-sequences_011-D90A_190413.txt",
                      "d90.lo")

df.d90hi <- mk_big_df("~/Documents/RepSeq/011-D90-Hi/1_Summary_011-D90B_190413.txt",
                      "~/Documents/RepSeq/011-D90-Hi/2_IMGT-gapped-nt-sequences_011-D90B_190413.txt",
                      "~/Documents/RepSeq/011-D90-Hi/4_IMGT-gapped-AA-sequences_011-D90B_190413.txt",
                      "d90.hi")

df.011 <- rbind(df.d7, df.d14hi, df.d14lo, df.d90hi, df.d90lo)
df.011$subject <- "011"

#012

df.d7.012 <- mk_big_df("~/Documents/RepSeq/012-d7PB/1_Summary_012-d7PB-1_190413.txt",
                       "~/Documents/RepSeq/012-d7PB/2_IMGT-gapped-nt-sequences_012-d7PB-1_190413.txt",
                       "~/Documents/RepSeq/012-d7PB/4_IMGT-gapped-AA-sequences_012-d7PB-1_190413.txt",
                       "d7.PB")

df.d14lo.012 <- mk_big_df("~/Documents/RepSeq/012-d-D14-Low/1_Summary_012-d-D14A_190413.txt",
                          "~/Documents/RepSeq/012-d-D14-Low/2_IMGT-gapped-nt-sequences_012-d-D14A_190413.txt",
                          "~/Documents/RepSeq/012-d-D14-Low/4_IMGT-gapped-AA-sequences_012-d-D14A_190413.txt",
                          "d14.lo")

df.d14hi.012 <- mk_big_df("~/Documents/RepSeq/012-d-D14-Hi/1_Summary_012-d-D14B_190413.txt",
                          "~/Documents/RepSeq/012-d-D14-Hi/2_IMGT-gapped-nt-sequences_012-d-D14B_190413.txt",
                          "~/Documents/RepSeq/012-d-D14-Hi/4_IMGT-gapped-AA-sequences_012-d-D14B_190413.txt",
                          "d14.hi")

df.d90lo.012 <- mk_big_df("~/Documents/RepSeq/012-d-D90-Low/1_Summary_012-d-D90A_190413.txt",
                          "~/Documents/RepSeq/012-d-D90-Low/2_IMGT-gapped-nt-sequences_012-d-D90A_190413.txt",
                          "~/Documents/RepSeq/012-d-D90-Low/4_IMGT-gapped-AA-sequences_012-d-D90A_190413.txt",
                          "d90.lo")

df.d90hi.012 <- mk_big_df("~/Documents/RepSeq/012-d-D90-High/1_Summary_012-d-D90B_190413.txt",
                          "~/Documents/RepSeq/012-d-D90-High/2_IMGT-gapped-nt-sequences_012-d-D90B_190413.txt",
                          "~/Documents/RepSeq/012-d-D90-High/4_IMGT-gapped-AA-sequences_012-d-D90B_190413.txt",
                          "d90.hi")
df.012 <- rbind(df.d7.012, df.d14hi.012, df.d14lo.012, df.d90hi.012, df.d90lo.012)

df.012$subject <- "012"

#007

df.d7.007 <- mk_big_df("~/Documents/RepSeq/SFA-007-D7-PB/1_Summary_SFA-007-D7-PB_261012.txt",
                       "~/Documents/RepSeq/SFA-007-D7-PB/2_IMGT-gapped-nt-sequences_SFA-007-D7-PB_261012.txt",
                       "~/Documents/RepSeq/SFA-007-D7-PB/4_IMGT-gapped-AA-sequences_SFA-007-D7-PB_261012.txt",
                       "d7.PB")

df.d14lo.007 <- mk_big_df("~/Documents/RepSeq/SFA-14-Low/1_Summary_SFA-14-Low_261012.txt",
                          "~/Documents/RepSeq/SFA-14-Low/2_IMGT-gapped-nt-sequences_SFA-14-Low_261012.txt",
                          "~/Documents/RepSeq/SFA-14-Low/4_IMGT-gapped-AA-sequences_SFA-14-Low_261012.txt",
                          "d14.lo")

df.d14hi.007 <- mk_big_df("~/Documents/RepSeq/SFA-14-High/1_Summary_SFA-14-igh_261012.txt",
                          "~/Documents/RepSeq/SFA-14-High/2_IMGT-gapped-nt-sequences_SFA-14-igh_261012.txt",
                          "~/Documents/RepSeq/SFA-14-High/4_IMGT-gapped-AA-sequences_SFA-14-igh_261012.txt",
                          "d14.hi")
df.d90hi.007 <- mk_big_df("~/Documents/RepSeq/SFA-90-High/1_Summary_SFA-90-High_261012.txt",
                          "~/Documents/RepSeq/SFA-90-High/2_IMGT-gapped-nt-sequences_SFA-90-High_261012.txt",
                          "~/Documents/RepSeq/SFA-90-High/4_IMGT-gapped-AA-sequences_SFA-90-High_261012.txt",
                          "d90.hi")

df.007 <- rbind(df.d7.007, df.d14lo.007, df.d14hi.007, df.d90hi.007)

df.007$subject <- "007"

big.df <- rbind(df.007, df.011, df.012)

#Format gene names
format.csv <- function(df){
  v.split <- colsplit(df$VGENE, "F", c("VGENE", "Vjunk"))
  j.split <- colsplit(df$JGENE, "F", c("JGENE", "Jjunk"))
  d.split <- colsplit(df$DGENE, "F", c("DGENE", "Djunk"))
  df$VGENE <- v.split$VGENE
  df$JGENE <- j.split$JGENE
  df$DGENE <- d.split$DGENE
  return(df)
}

big.df <- format.csv(big.df)

#remove sequences with not CDR3 identified
df.mod <- subset(big.df, big.df$CDR3_AA != "")

#remove sequences with no JGENE
df.mod <- subset(df.mod, df.mod$JGENE != "")

#remove sequences with CDR3 < 5
#df.mod$CDR3_len <- as.numeric(df.mod$CDR3_len)
#df.mod <- subset(df.mod, df.mod$CDR3_len < 5)

#make subgroup ID column
df.mod$subgroup <- with(df.mod, paste0(VGENE, JGENE, CDR3_len))

#remove sequences with no CDR3 length identified
df.out <- subset(df.mod, df.mod$CDR3_len != "X")

#csv.011 <- subset(df.mod, df.mod$subject == "011")
#write.csv(csv.011, file = "~/Documents/RepSeq/subgroup011.csv")

#csv.012 <- subset(df.mod, df.mod$subject == "012")
#write.csv(csv.012, file = "~/Documents/RepSeq/subgroup012.csv")

#csv.007 <- subset(df.mod, df.mod$subject == "007")
#write.csv(csv.007, file = "~/Documents/RepSeq/subgroup007.csv")

csv.011 <- subset(df.out, df.out$subject == "011")
write.csv(csv.011, file = "~/Documents/RepSeq2/IMGT_011.csv")

csv.012 <- subset(df.out, df.out$subject == "012")
write.csv(csv.012, file = "~/Documents/RepSeq2/IMGT_012.csv")

csv.007 <- subset(df.out, df.out$subject == "007")
write.csv(csv.007, file = "~/Documents/RepSeq2/IMGT_007.csv")


#Make sqlite table
db <- dbConnect(SQLite(), dbname="~/Documents/RepSeq2/IMGT_parsed.sqlite")
dbWriteTable(conn = db, name = "IMGT_007", value = csv.007,
             row.names = FALSE, overwrite = TRUE)
dbWriteTable(conn = db, name = "IMGT_011", value = csv.011,
             row.names = FALSE, overwrite = TRUE)
dbWriteTable(conn = db, name = "IMGT_012", value = csv.012,
             row.names = FALSE, overwrite = TRUE)

dbDisconnect(db)


