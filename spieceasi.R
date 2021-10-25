
source("~/alberts.R")

library(library(SpiecEasi))



g1 <- cha(read.table("fast_a/a_beforePEG_present_50.txt",sep="\t", header=F))
rownames(g1) <- g1[,1]
g1 <- g1[,-1]


g2 <- cha(read.table("fast_b/b_afterPEG_present_50.txt",sep="\t", header=F))
rownames(g2) <- g2[,1]
g2 <- g2[,-1]
 
g3 <- cha(read.table("fast_c/c_6months_present_50.txt",sep="\t", header=F))
rownames(g3) <- g3[,1]
g3 <- g3[,-1]

g4 <- cha(read.table("fast_d/d_healthy_present_50.txt",sep="\t", header=F))
rownames(g4) <- g4[,1]
g4 <- g4[,-1]


## samples go into the rows
## OTU or ASV go into the columns

g1 <- t(g1)
g2 <- t(g2)
g3 <- t(g3)
g4 <- t(g4)



### normalization

depths <- rowSums(g1)
g1.n  <- t(apply(g1, 1, norm_to_total))
g1.cs <- round(g1.n * min(depths))

depths <- rowSums(g2)
g2.n  <- t(apply(g2, 1, norm_to_total))
g2.cs <- round(g2.n * min(depths))

depths <- rowSums(g3)
g3.n  <- t(apply(g3, 1, norm_to_total))
g3.cs <- round(g3.n * min(depths))

depths <- rowSums(g4)
g4.n  <- t(apply(g4, 1, norm_to_total))
g4.cs <- round(g4.n * min(depths))



### apply SpiecEasi glasso method

s1 <- spiec.easi(as.matrix(g1.cs), method='glasso', lambda.min.ratio=1e-2,
                          nlambda=30, pulsar.params=list(rep.num=50))
s2 <- spiec.easi(as.matrix(g2.cs), method='glasso', lambda.min.ratio=1e-2,
                          nlambda=30, pulsar.params=list(rep.num=50))
s3 <- spiec.easi(as.matrix(g3.cs), method='glasso', lambda.min.ratio=1e-2,
                          nlambda=30, pulsar.params=list(rep.num=50))
s4 <- spiec.easi(as.matrix(g4.cs), method='glasso', lambda.min.ratio=1e-2,
                          nlambda=30, pulsar.params=list(rep.num=50))

gs1 <- adj2igraph(getRefit(s1))
gs2 <- adj2igraph(getRefit(s2))
gs3 <- adj2igraph(getRefit(s3))
gs4 <- adj2igraph(getRefit(s4))



### apply sparCC
library(Matrix)

## Define arbitrary threshold for SparCC correlation matrix for the graph

g1.sparcc <- sparcc(g1.cs)
g1.graph <- abs(g1.sparcc$Cor) >= 0.75
diag(g1.graph) <- 0
g1.graph <- Matrix(g1.graph, sparse=TRUE)

g2.sparcc <- sparcc(g2.cs)
g2.graph <- abs(g2.sparcc$Cor) >= 0.75
diag(g2.graph) <- 0
g2.graph <- Matrix(g2.graph, sparse=TRUE)

g3.sparcc <- sparcc(g3.cs)
g3.graph <- abs(g3.sparcc$Cor) >= 0.75
diag(g3.graph) <- 0
g3.graph <- Matrix(g3.graph, sparse=TRUE)

g4.sparcc <- sparcc(g4.cs)
g4.graph <- abs(g4.sparcc$Cor) >= 0.75
diag(g4.graph) <- 0
g4.graph <- Matrix(g4.graph, sparse=TRUE)



g5.sparcc <- sparcc(g1.cs)
g5.graph <- abs(g5.sparcc$Cor) >= 0.8
diag(g5.graph) <- 0
g5.graph <- Matrix(g5.graph, sparse=TRUE)

g6.sparcc <- sparcc(g2.cs)
g6.graph <- abs(g6.sparcc$Cor) >= 0.8
diag(g6.graph) <- 0
g6.graph <- Matrix(g6.graph, sparse=TRUE)

g7.sparcc <- sparcc(g3.cs)
g7.graph <- abs(g7.sparcc$Cor) >= 0.8
diag(g7.graph) <- 0
g7.graph <- Matrix(g7.graph, sparse=TRUE)

g8.sparcc <- sparcc(g4.cs)
g8.graph <- abs(g8.sparcc$Cor) >= 0.8
diag(g8.graph) <- 0
g8.graph <- Matrix(g8.graph, sparse=TRUE)




g9.sparcc <- sparcc(g1.cs)
g9.graph <- abs(g9.sparcc$Cor) >= 0.85
diag(g9.graph) <- 0
g9.graph <- Matrix(g9.graph, sparse=TRUE)

g10.sparcc <- sparcc(g2.cs)
g10.graph <- abs(g10.sparcc$Cor) >= 0.85
diag(g10.graph) <- 0
g10.graph <- Matrix(g10.graph, sparse=TRUE)

g11.sparcc <- sparcc(g3.cs)
g11.graph <- abs(g11.sparcc$Cor) >= 0.85
diag(g11.graph) <- 0
g11.graph <- Matrix(g11.graph, sparse=TRUE)

g12.sparcc <- sparcc(g4.cs)
g12.graph <- abs(g12.sparcc$Cor) >= 0.85
diag(g12.graph) <- 0
g12.graph <- Matrix(g12.graph, sparse=TRUE)




ig1 <- adj2igraph(g1.graph)
ig2 <- adj2igraph(g2.graph)
ig3 <- adj2igraph(g3.graph)
ig4 <- adj2igraph(g4.graph)

ig5 <- adj2igraph(g5.graph)
ig6 <- adj2igraph(g6.graph)
ig7 <- adj2igraph(g7.graph)
ig8 <- adj2igraph(g8.graph)

ig9  <- adj2igraph(g9.graph)
ig10 <- adj2igraph(g10.graph)
ig11 <- adj2igraph(g11.graph)
ig12 <- adj2igraph(g12.graph)



library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(g1, 1))+6
am.coord <- layout.fruchterman.reingold(gs1)

par(mfrow=c(2,4))
plot(gs1, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="SpiecEasi glasso beforePEG")
plot(gs2, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="afterPEG")
plot(gs3, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="6months")
plot(gs4, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="healthy")

plot(ig1, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparCC  0.75  beforePEG")
plot(ig2, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="afterPEG")
plot(ig3, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="6months")
plot(ig4, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="healthy")








par(mfrow=c(2,4))

plot(ig5, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparCC  0.80  beforePEG")
plot(ig6, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="afterPEG")
plot(ig7, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="6months")
plot(ig8, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="healthy")

plot(ig9, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparCC  0.85   beforePEG")
plot(ig10,layout=am.coord, vertex.size=vsize, vertex.label=NA, main="afterPEG")
plot(ig11,layout=am.coord, vertex.size=vsize, vertex.label=NA, main="6months")
plot(ig12,layout=am.coord, vertex.size=vsize, vertex.label=NA, main="healthy")



### get the annotation of the edges



asv <- cha(read.table("fast_a/edited_ASVs_taxonomy.pick.txt",sep="\t", header=T))
rownames(asv) <- asv[,1]
asv <- asv[,-1]
asv$all <- NA
for (j in 1:nrow(asv)) {
    asv$all[j] <- paste(asv[j,1], asv[j,2], asv[j,3], asv[j,4], asv[j,5], asv[j,6], asv[j,7], sep="_")
}





names1 <- colnames(g1.cs)
names2 <- colnames(g2.cs)
names3 <- colnames(g3.cs)
names4 <- colnames(g4.cs)

n1 <- cha(data.frame(net="SpeacEasi before PEG", as_edgelist(gs1)))
n1$asv1 <- names1[n1[,2]]
n1$asv2 <- names1[n1[,3]]
n1$species1 <- asv[n1$asv1,8]
n1$species2 <- asv[n1$asv2,8]

n2 <- cha(data.frame(net="SpeacEasi after PEG", as_edgelist(gs2)))
n2$asv1 <- names2[n2[,2]]
n2$asv2 <- names2[n2[,3]]
n2$species1 <- asv[n2$asv1,8]
n2$species2 <- asv[n2$asv2,8]

n3 <- cha(data.frame(net="SpeacEasi 6 months", as_edgelist(gs3)))
n3$asv1 <- names3[n3[,2]]
n3$asv2 <- names3[n3[,3]]
n3$species1 <- asv[n3$asv1,8]
n3$species2 <- asv[n3$asv2,8]

n4 <- cha(data.frame(net="SpeacEasi healthy", as_edgelist(gs4)))
n4$asv1 <- names4[n4[,2]]
n4$asv2 <- names4[n4[,3]]
n4$species1 <- asv[n4$asv1,8]
n4$species2 <- asv[n4$asv2,8]

result <- cha(rbind(n1,n2,n3,n4))



n1 <- cha(data.frame(net="sparCC 0.75 before PEG", as_edgelist(ig1)))
n1$asv1 <- names1[n1[,2]]
n1$asv2 <- names1[n1[,3]]
n1$species1 <- asv[n1$asv1,8]
n1$species2 <- asv[n1$asv2,8]

n2 <- cha(data.frame(net="sparCC 0.75 after PEG", as_edgelist(ig2)))
n2$asv1 <- names2[n2[,2]]
n2$asv2 <- names2[n2[,3]]
n2$species1 <- asv[n2$asv1,8]
n2$species2 <- asv[n2$asv2,8]

n3 <- cha(data.frame(net="sparCC 0.75 6 months", as_edgelist(ig3)))
n3$asv1 <- names3[n3[,2]]
n3$asv2 <- names3[n3[,3]]
n3$species1 <- asv[n3$asv1,8]
n3$species2 <- asv[n3$asv2,8]

n4 <- cha(data.frame(net="sparCC 0.75 healthy", as_edgelist(ig4)))
n4$asv1 <- names4[n4[,2]]
n4$asv2 <- names4[n4[,3]]
n4$species1 <- asv[n4$asv1,8]
n4$species2 <- asv[n4$asv2,8]


result <- cha(rbind(result,n1,n2,n3,n4))

write.table(result,"network_components.txt",sep="\t",quote=F,row.names=F)


