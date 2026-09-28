#!/usr/bin/env Rscript
# Paper-style Figure 6d: draw all significant terms from the regenerated table.
# Geometry, palette, font and count-size rule follow the provided final figure.
suppressPackageStartupMessages({library(grid);library(scales);library(ragg);library(svglite)})
argv <- commandArgs(TRUE)
if(length(argv)<2)stop("Usage: Rscript draw_Figure6d.R GO_table.tsv output_stem")
go <- read.delim(argv[1],check.names=FALSE)
out <- argv[2]; dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
dd <- go[go$html_displayed,,drop=FALSE]
dd$class <- ifelse(dd$SourceType=='Common','Shared',dd$SourceType)
dd$name <- dd$Description;dd$ontology <- dd$ONTOLOGY;dd$display_Count<-dd$Count
stopifnot(nrow(dd)==131,all(is.finite(dd$display_Count)))
dd$diameter_mm <- scales::rescale(sqrt(dd$display_Count),to=c(1,3.35))
s<-list(D=list(significant_Protein=sum(go$significant_Protein),significant_RNA=sum(go$significant_RNA),shared=sum(go$significant_in_both)))
SX<-.84;SY<-.74
U<-function(x,units,...)grid::unit(if(units=='mm')x*SX else x,units,...)
Y<-function(x,units,...)grid::unit(if(units=='mm')x*SY else x,units,...)
fontgp<-function(col='#202020',size=8)gpar(fontfamily='Arial',fontface='bold',fontsize=size,col=col,lineheight=1.05)
tx<-function(label,x,y,just='centre',col='#202020',size=8,bg=FALSE){g<-textGrob(label,x=U(x,'mm'),y=Y(y,'mm'),just=just,gp=fontgp(col,size));if(bg)grid.rect(U(x,'mm'),Y(y,'mm'),width=grobWidth(g)+U(.8,'mm'),height=grobHeight(g)+U(.5,'mm'),just=just,gp=gpar(fill='white',col=NA));grid.draw(g)}
ln<-function(x,y,x1,y1,col='#c8c8c8',lwd=.6,arrow=FALSE){grid.segments(U(x,'mm'),Y(y,'mm'),U(x1,'mm'),Y(y1,'mm'),gp=gpar(col=col,lwd=lwd),arrow=if(arrow)grid::arrow(length=U(1.2,'mm'),type='closed')else NULL)}
dot<-function(x,y,r=1.4,fill='#afd0de',col='#647781',lwd=.7)grid.circle(U(x,'mm'),Y(y,'mm'),r=U(r,'mm'),gp=gpar(fill=fill,col=col,lwd=lwd))
wrap<-function(x,n=25)paste(strwrap(x,width=n),collapse='\n')
letter<-function(z,h)tx(z,3,h-4,just=c('left','top'),size=12)
wt<-function(d,n)write.table(d,paste0(out,'_',n,'.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
# D: multi-layer arcs; size reproduces the documented original Count rule.
dp<-data.frame();category<-data.frame()
for(o in c('Protein','RNA'))for(on in c('BP','CC','MF')){
 z<-dd[dd$class==o&dd$ontology==on,];if(!nrow(z))next;z<-z[order(z$ID),];center<-if(o=='Protein')44 else 195
 angle<-switch(on,BP=90,MF=270,CC=if(o=='Protein')180 else 0);spread<-switch(on,BP=90,MF=80,CC=100)
 radii<-if(on=='CC')c(17,26,35) else c(20,30,40);capacity<-c(6,10,16);remaining<-nrow(z);counts<-pmin(capacity,pmax(0,remaining-c(0,cumsum(capacity)[1:2])))
 if(sum(counts)<remaining){counts[3]<-remaining-sum(counts[1:2])};j<-0
 for(k in which(counts>0)){n<-counts[k];ang<-(angle+if(n==1)0 else seq(-spread/2,spread/2,length.out=n))*pi/180;ii<-seq.int(j+1,j+n);dp<-rbind(dp,data.frame(ID=z$ID[ii],name=z$name[ii],class=o,ontology=on,x=center+radii[k]*cos(ang),y=57+radii[k]*sin(ang),diameter_mm=z$diameter_mm[ii],display_Count=z$display_Count[ii]));j<-j+n}
 category<-rbind(category,data.frame(class=o,ontology=on,x=if(on=='CC'){if(o=='Protein')5 else 235}else center,y=switch(on,BP=110,MF=7,CC=57),n=nrow(z)))
}
for(on in c('BP','CC')){z<-dd[dd$class=='Shared'&dd$ontology==on,];z<-z[order(z$ID),];ij<-seq_len(nrow(z));dp<-rbind(dp,data.frame(ID=z$ID,name=z$name,class='Shared',ontology=on,x=if(on=='BP')c(98,105)[ij]else c(102,97,103)[ij],y=if(on=='BP')c(94,77)[ij]else c(43,28,13)[ij],diameter_mm=z$diameter_mm,display_Count=z$display_Count));category<-rbind(category,data.frame(class='Shared',ontology=on,x=139,y=if(on=='BP')110 else 7,n=nrow(z)))}
stopifnot(nrow(dp)==nrow(dd),!anyDuplicated(dp$ID),all(is.finite(dp$x)),all(is.finite(dp$y)));wt(dp,'FIG6D_NODE_POSITIONS');wt(category,'FIG6D_CATEGORY_POSITIONS')
lnD<-function(x,y,x1,y1,...)ln(12+.9*x,y,12+.9*x1,y1,...)
dotD<-function(x,y,...)dot(12+.9*x,y,...)
txD<-function(label,x,y,...)tx(label,12+.9*x,y,...)
plotD<-function(){letter('d',135)
 for(i in seq_len(nrow(dp))){z<-dp[i,];ca<-category[category$class==z$class&category$ontology==z$ontology,];hx<-switch(z$class,Protein=44,RNA=195,Shared=120);lnD(hx,57,z$x,z$y,col='#edbbc1',lwd=.65);lnD(ca$x,ca$y,z$x,z$y,col='#a3a3a3',lwd=.65)}
 for(i in seq_len(nrow(dp))){z<-dp[i,];dotD(z$x,z$y,z$diameter_mm/2);if(z$class=='Shared')txD(wrap(z$name,29),z$x+4,z$y,just='left',bg=TRUE)}
 for(i in seq_len(nrow(category))){z<-category[i,];dotD(z$x,z$y,1.8,'#F3AD28');col<-if(z$class=='Shared')'#BE2737'else'#305D96';label<-sprintf('%s\n(n = %d)',z$ontology,z$n)
  if(z$ontology=='BP')txD(label,z$x,z$y+8,col=col)
  else if(z$ontology=='MF')txD(label,z$x,z$y-8,col=col)
  else if(z$class=='Shared')txD(label,149,7,just='left',col=col)
  else txD(label,if(z$class=='Protein')-2 else 242,49,just='centre',col=col)}
 for(o in c('Protein','RNA','Shared')){x<-switch(o,Protein=44,RNA=195,Shared=120);dotD(x,57,2.8,switch(o,Protein='#399359',RNA='#7656AA',Shared='#B82D3C'));txD(o,x,63,col=if(o=='Shared')'#BE2737'else'#305D96')}
 txD(sprintf('Protein %d / RNA %d / shared %d / unique union %d',s$D$significant_Protein,s$D$significant_RNA,s$D$shared,nrow(dd)),120,-5)
 txD('GO hit count (shared: assay sum)',75,-13,just='right')
 for(i in 1:3){v<-c(10,100,300)[i];diam<-scales::rescale(sqrt(v),to=c(1,3.35),from=range(sqrt(dd$display_Count)));x<-95+(i-1)*27;dotD(x,-13,diam/2);txD(as.character(v),x+4,-13,just='left')}
}

render <- function(){grid.newpage();grid.rect(gp=gpar(fill='white',col=NA));pushViewport(viewport(x=unit(4.2,'mm'),y=unit(13.5,'mm'),width=unit(201.6,'mm'),height=unit(100,'mm'),just=c('left','bottom'),clip='off'));plotD();popViewport()}
cairo_pdf(paste0(out,'.pdf'),width=210/25.4,height=118/25.4,family='Arial');render();dev.off()
svglite::svglite(paste0(out,'.svg'),width=210/25.4,height=118/25.4);render();dev.off()
ragg::agg_png(paste0(out,'.png'),width=210,height=118,units='mm',res=600,background='white');render();dev.off()
wt(dd,'source_data')
cat('131 statistically eligible terms redrawn; all source values originate in regenerated GO table.\n')
