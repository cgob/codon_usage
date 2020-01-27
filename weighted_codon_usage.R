require('xlsx')
require('stringi')
require('plyr')
require('Biostrings')
require('ggplot2')
require('ggrepel')
require('RColorBrewer')
require('gridExtra')

dat=read.xlsx("/home/cgobet/codon_usage/Data/2wk_Ribotag_GarsC201R_NegCtrlvsWT.xlsx",1)
rownames(dat)=dat$Ensembl.Gene.ID
dat.2=read.xlsx("/home/cgobet/codon_usage/Data/8wk_Ribotag_Gars249_NegCtrlvsWT.xlsx",1)
rownames(dat.2)=dat.2$Ensembl.Gene.ID
dat.3=read.xlsx("/home/cgobet/codon_usage/Data/GSE45684_B6maleLiverSamples_rawCounts_processedmatrix_forRobBurgess.xlsx",1)
rownames(dat.3)=dat.3$Ensembl.Gene.ID
modi=0

codon = readDNAStringSet("/home/cgobet/codon_usage/Data/Mus_musculus.GRCm38.cds.all.fa")
codon=as.data.frame(codon)
seq_name = sapply(strsplit(sapply(strsplit(rownames(codon),"gene:"),"[[",2),".",fixed=T),"[[",1)
if(modi !=0){
df <- data.frame(seq_name, substr(codon$x,4,modi+3),stringsAsFactors = F)
df=df[nchar(df[,2])==modi,]
}else{ 
df <- data.frame(seq_name, codon$x,stringsAsFactors = F)
}

names(df)[2]='codon.x'
b=sapply(df$codon.x,function(x) unlist(stri_extract_all_regex(x, '.{1,3}')))
if(modi!=0){
  b.table=apply(b,2,function(x) data.frame(t(as.matrix(table(x)))))
}else{ 
  b.table=lapply(b,function(x) data.frame(t(as.matrix(table(x)))))
}

names(b.table)=""

b.df=do.call(rbind.fill,b.table)
b.df[is.na(b.df)]=0
b.df=b.df[,nchar(colnames(b.df))==3]
b.df$gene.name=df$seq_name
ss.name=split(1:nrow(b.df),b.df$gene.name)
b.df.gen= sapply(ss.name,function(x) colMeans(b.df[x,-ncol(b.df)]))
b.df.gen=t(b.df.gen)

dat$mean.neg=log2(1+rowMeans(dat[,grep('Neg',names(dat))]))
dat$mean.Rpl=log2(1+rowMeans(dat[,grep('Rpl',names(dat))]))

dat.2$mean.neg=log2(1+rowMeans(dat.2[,grep('Neg',names(dat.2))]))
dat.2$mean.Rpl=log2(1+rowMeans(dat.2[,grep('Rpl',names(dat.2))]))
 
dat.3$mean.Rpl=log2(1+rowMeans(dat.3[,-c(1,2)]))
rownames(dat.3)=dat.3$Ensembl.Gene.ID
cod.ss=as.data.frame(b.df.gen)

cod.ss=cod.ss[,-grep('TAG|TGA',names(cod.ss))]
GG.df=cod.ss

compute_cu = function(dati){
 ii=intersect(rownames(GG.df),rownames(dati))
 GG.df.sel=GG.df[ii,] 
 dati=dati[ii,]
 tt=GG.df.sel[,1:61]
 tt=as.data.frame(tt)
 tt$Mean=dati$mean.Rpl
 tt=sweep(tt[,1:61],1,2^tt$Mean,FUN="*")

 cod.u=colMeans(tt[,1:61])
 co.du=cod.u/sum(cod.u)

}

c.1=compute_cu(dat)
c.2=compute_cu(dat.2)
c.3=compute_cu(dat.3)

aa=GENETIC_CODE[names(c.1)]

c1.aa=tapply(c.1,aa,sum)
c2.aa=tapply(c.2,aa,sum)
c3.aa=tapply(c.3,aa,sum)



V=data.frame(c.1,c.2,c.3,aa=GENETIC_CODE[names(c.1)],codon=names(c.1))

g1=ggplot(V,aes(x=c.1/c.3,y=c.2/c.3)) + geom_point(size=3,aes(col=aa)) +
  geom_abline(intercept=0,slope=1) +
  geom_text_repel(aes(label=codon,col=aa)) +
  theme_bw()+theme(axis.title=element_text(size=12),axis.text=element_text(size=10,angle=90),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),aspect.ratio=1) +
  xlab("Weighted Codon Usage 2wk_Ribotag_GarsC201R") + ylab("Weighted Codon Usage 8wk_Ribotag_Gars249") + scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(20))
g2=ggplot(V,aes(x=c.1/c.3,y=c.2/c.3)) + geom_point(size=3,aes(col=aa)) +
  geom_abline(intercept=0,slope=1) +
  geom_text_repel(aes(label=codon,col=aa)) +
  theme_bw()+theme(axis.title=element_text(size=12),axis.text=element_text(size=10,angle=90),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),aspect.ratio=1) +
  xlab("Weighted Codon Usage 2wk_Ribotag_GarsC201R/liver") + ylab("Weighted Codon Usage 8wk_Ribotag_Gars249/liver") + scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(20))


V=data.frame(c1.aa,c2.aa,c3.aa,aa=names(c1.aa))
g3=ggplot(V,aes(x=c1.aa,y=c2.aa)) + geom_point(size=3,aes(col=aa)) +
  geom_abline(intercept=0,slope=1) +
  geom_text_repel(aes(label=aa,col=aa)) +
  theme_bw()+theme(axis.title=element_text(size=12),axis.text=element_text(size=10,angle=90),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),aspect.ratio=1) +
  xlab("Weighted AA Usage 2wk_Ribotag_GarsC201R") + ylab("Weighted AA Usage 8wk_Ribotag_Gars249") + scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(20))
g4=ggplot(V,aes(x=c1.aa/c3.aa,y=c2.aa/c3.aa)) + geom_point(size=3,aes(col=aa)) +
  geom_abline(intercept=0,slope=1) +
  geom_text_repel(aes(label=aa,col=aa)) +
  theme_bw()+theme(axis.title=element_text(size=12),axis.text=element_text(size=10,angle=90),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),aspect.ratio=1) +
  xlab("weighted AA Usage (2wk_Ribotag_GarsC201R/Liver)") + ylab("weighted AA Usage (8wk_Ribotag_Gars249/liver)") + scale_color_manual(values=colorRampPalette(brewer.pal(12, "Paired"))(20))
 
pdf("codon_usage.pdf")
grid.arrange(g1,g2,ncol=2)
grid.arrange(g3,g4,ncol=2)
dev.off()
