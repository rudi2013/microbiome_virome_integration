

source("~/alberts.R")


d <- cha(read.table("a_beforePEG.tsv",sep="\t", header=F))
rownames(d) <- d[,1]
d <- d[,-1]
d$total <- apply(d>0,1,sum)
d <- d[d$total>=7,]
dim(d)

collection <- rownames(d)


d2 <- cha(read.table("../fast_b/b_afterPEG.tsv",sep="\t", header=F))
rownames(d2) <- d2[,1]
d2 <- d2[,-1]
d2$total <- apply(d2>0,1,sum)
d2 <- d2[d2$total>=7,]
dim(d2)

collection <- c(collection, rownames(d2))



d3 <- cha(read.table("../fast_c/c_6months.tsv",sep="\t", header=F))
rownames(d3) <- d3[,1]
d3 <- d3[,-1]
d3$total <- apply(d3>0,1,sum)
d3 <- d3[d3$total>=7,]
dim(d3)

collection <- c(collection, rownames(d3))



d4 <- cha(read.table("../fast_d/d_healthy.tsv",sep="\t", header=F))
rownames(d4) <- d4[,1]
d4 <- d4[,-1]
d4$total <- apply(d4>0,1,sum)
d4 <- d4[d4$total>=13,]
dim(d4)

collection <- c(collection, rownames(d4))

collection <- sort(unique(collection))
length(collection)


# read in again to take the AVSs in the collection in each cohort

d <- cha(read.table("a_beforePEG.tsv",sep="\t", header=F))
rownames(d) <- d[,1]
d <- d[,-1]

d2 <- cha(read.table("../fast_b/b_afterPEG.tsv",sep="\t", header=F))
rownames(d2) <- d2[,1]
d2 <- d2[,-1]

d3 <- cha(read.table("../fast_c/c_6months.tsv",sep="\t", header=F))
rownames(d3) <- d3[,1]
d3 <- d3[,-1]

d4 <- cha(read.table("../fast_d/d_healthy.tsv",sep="\t", header=F))
rownames(d4) <- d4[,1]
d4 <- d4[,-1]


d <- d[collection, ]
d2 <- d2[collection, ]
d3 <- d3[collection, ]
d4 <- d4[collection, ]

dim(d)
dim(d2)
dim(d3)
dim(d4)

write.table(d,"a_beforePEG_present_50.txt",sep="\t",quote=F)

write.table(d2,"../fast_b/b_afterPEG_present_50.txt",sep="\t",quote=F)

write.table(d3,"../fast_c/c_6months_present_50.txt",sep="\t",quote=F)

write.table(d4,"../fast_d/d_healthy_present_50.txt",sep="\t",quote=F)

## add #OTU ID by hand in first cell
















