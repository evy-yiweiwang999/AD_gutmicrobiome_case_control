geu <- read.csv("/Users/yw794/Yale_Posd/AD_gutmicrobiome/csvfiles/genus.csv", header = T, row.names = 1)
geu_otu <- geu[, -c(1,2)]
geu_meta <- geu[,c(1,2)]
geu_otu <- geu_otu[rownames(geu_meta),]
zerocol <- which(colSums(geu_otu)==0)
if(length(zerocol) >0 ){
  geu_otu <- geu_otu[,-zerocol];
}
dim(geu_otu)

X1 <- geu_otu
group1 <-   geu_meta$Group
Z1 <- ifelse(group1=="AD",1,0)

geu_meta2 <- geu_meta[geu_meta$Group=="AD",]
geu_otu2 <- geu_otu[rownames(geu_meta2),]

zerocol <- which(colSums(geu_otu2)==0)
if(length(zerocol) >0 ){
  geu_otu2 <- geu_otu2[,-zerocol];
}
dim(geu_otu2)

X2 <- geu_otu2
Z2 <- c(rep(0,round(length(geu_meta2$Group)/2)),rep(1,length(geu_meta2$Group)-round(length(geu_meta2$Group)/2)))
group2 <- ifelse(Z2==0,"AD","Normal")

group <- list(group1,group2);
X <- list(as.matrix(X1),as.matrix(X2))
Z <- list(Z1,Z2)
n <- length(X)
## Fitting models  
ZINB <- list(data1=NULL,data2=NULL); ZINB_cov <- list(data1=NULL,data2=NULL)
for(i in 1:n){
  ZINB[[i]] <- tryCatch({ZIPPCApn(X[[i]],rank = T)},
                        error=function(e){ NULL})
  ZINB_cov[[i]] <- tryCatch({ZIPPCApn(X[[i]],Z[[i]],rank = T)},
                            error=function(e){ NULL})
}
adonis_result_dis = adonis2(dist(ZINB[[i]]$lvs$factor_scores2)~group[[i]],method = "euclidean")
p1 <- ggplot(data.frame(ZINB[[i]]$lvs$factor_scores2),aes(x=ZINB[[i]]$lvs$factor_scores2[,1], y=ZINB[[i]]$lvs$factor_scores2[,2],colour=as.factor(group[[i]]))) +
  geom_point( size = 2) + scale_colour_brewer(palette = "Set1")+
  xlab("mbDenoise-zinb F1")+ylab("mbDenoise-zinb F2")+ labs(colour="Group") +
  theme_bw()+
  stat_ellipse(aes(fill=factor(group[[i]],level=c("AD","Normal"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)
X_new <- log2(ZINB[[i]]$muz+1)
ZINB_pca[[i]] <- prcomp(X_new)
dist_bray <- vegan::vegdist(X_new, method="bray")
ZINB_pcoa[[i]] <- ape::pcoa(dist_bray)
set.seed(4)
ZINB_tsne[[i]] <- Rtsne::Rtsne(X_new,initial_dims = ncol(X_new),perplexity=round((nrow(X_new)-2)/3))

adonis_result_dis = adonis2(dist(ZINB_pca[[i]]$x[,1:2])~group[[i]],method = "euclidean")
p2 <- ggplot(data.frame(ZINB_pca[[i]]$x[,1:2]),aes(ZINB_pca[[i]]$x[,1], y=ZINB_pca[[i]]$x[,2],colour=as.factor(group[[i]]))) +
  geom_point( size = 2) + scale_colour_brewer(palette = "Set1")+
  xlab("mbDenoise-zinb_pca F1")+ylab("mbDenoise-zinb_pca F2")+ labs(colour="Group") +
  theme_bw()+
  stat_ellipse(aes(fill=factor(group[[i]],level=c("AD","Normal"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)
adonis_result_dis = adonis2(dist(ZINB_pcoa[[i]]$vectors[,1:2])~group[[i]],method = "euclidean")
p3 <- ggplot(data.frame(ZINB_pcoa[[i]]$vectors[,1:2]),aes(ZINB_pcoa[[i]]$vectors[,1], y=ZINB_pcoa[[i]]$vectors[,2],colour=as.factor(group[[i]]))) +
  geom_point( size = 2) + scale_colour_brewer(palette = "Set1")+
  xlab("mbDenoise-zinb_pcoa F1")+ylab("mbDenoise-zinb_pcoa F2")+ labs(colour="Group") +
  theme_bw()+
  stat_ellipse(aes(fill=factor(group[[i]],level=c("AD","Normal"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)
adonis_result_dis = adonis2(dist(ZINB_tsne[[i]]$Y)~group[[i]],method = "euclidean")
p4 <- ggplot(data.frame(ZINB_tsne[[i]]$Y),aes(ZINB_tsne[[i]]$Y[,1], y=ZINB_tsne[[i]]$Y[,2],colour=as.factor(group[[i]]))) +
  geom_point(size = 2) + scale_colour_brewer(palette = "Set1")+
  xlab("mbDenoise-zinb_tsne F1")+ylab("mbDenoise-zinb_tsne F2")+ labs(colour="Group") +
  theme_bw()+
  stat_ellipse(aes(fill=factor(group[[i]],level=c("AD","Normal"))),type = "norm", geom = "polygon",alpha= 0.1,show.legend = F,linetype=2)
