#!/usr/bin/env Rscript
# Draw Figure 1d from the case-insensitive identifier recount.
# Optional third argument assembles the panel into the supplied Figure 1 PNG.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop('Usage: Rscript draw_Figure1d.R counts.tsv output_directory [Figure1_base.png]')
d <- read.delim(args[1],check.names=FALSE)
out <- args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(identical(as.character(d$molecular_category),c('mRNA','lncRNA','miRNA','Proteome','Phosphoproteome','Metabolites')))
library(grid)
W<-554;H<-540;ppi<-1269/(183/25.4);cols<-c('#89B3BE','#A3B5CD','#D6A2A5','#C9B7AD','#E5AF81','#97B78D')
# Coordinates preserve the supplied panel's size, category order and visual style.
draw_panel<-function(){
 grid.newpage();pushViewport(viewport(xscale=c(0,W),yscale=c(H,0)))
 grid.rect(gp=gpar(fill='white',col=NA))
 textp<-function(label,x,y,size=8.2,angle=0,just='centre')grid.text(label,x=unit(x,'native'),y=unit(y,'native'),rot=angle,just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface='bold'))
 linep<-function(x,y)grid.lines(unit(x,'native'),unit(y,'native'),gp=gpar(col='black',lwd=.65))
 base<-370;top<-55;left<-123;right<-543;maximum<-22000
 y<-function(v)base-(base-top)*v/maximum
 linep(c(left,left),c(y(20000),base));linep(c(left,right),c(base,base))
 for(v in seq(0,20000,5000)){linep(c(left-4,left),rep(y(v),2));textp(format(v,big.mark=',',scientific=FALSE,trim=TRUE),left-8,y(v),just='right')}
 textp('Number of distinct identifiers',24,222,angle=90)
 centers<-seq(left+40,right-40,length.out=6)
 for(i in seq_len(nrow(d))){
  v<-d$distinct_identifiers[i];h<-base-y(v)
  grid.rect(x=unit(centers[i],'native'),y=unit(base-h/2,'native'),width=unit(49,'native'),height=unit(h,'native'),gp=gpar(fill=cols[i],col=NA))
  textp(format(v,big.mark=',',scientific=FALSE,trim=TRUE),centers[i],y(v)-13)
  linep(rep(centers[i],2),c(base,base+4))
  textp(d$molecular_category[i],centers[i]+10,base+16,angle=45,just='right')
 }
 popViewport()
}
for(ext in c('png','tiff','pdf','svg')){
 f<-file.path(out,paste0('Figure1d.',ext))
 if(ext=='png')ragg::agg_png(f,width=W/ppi,height=H/ppi,units='in',res=600)
 if(ext=='tiff')ragg::agg_tiff(f,width=W/ppi,height=H/ppi,units='in',res=600,compression='lzw')
 if(ext=='pdf')cairo_pdf(f,width=W/ppi,height=H/ppi,family='Arial')
 if(ext=='svg')svg(f,width=W/ppi,height=H/ppi,family='Arial')
 draw_panel();dev.off()
}
# Render a native-resolution panel so the other assembled panels stay pixel-identical.
native<-file.path(out,'panel_native.png')
ragg::agg_png(native,width=W,height=H,units='px',res=ppi);draw_panel();dev.off()
if(length(args)>=3){
 im<-png::readPNG(args[3]);patch<-png::readPNG(native)
 stopifnot(dim(im)[1]==1795,dim(im)[2]==1269)
 # Zero-based panel rectangle: x=715, y=570, width=554, height=540.
 im[571:1110,716:1269,1:3]<-patch[,,1:3]
 png::writePNG(im,file.path(out,'Figure1_updated.png'),dpi=ppi)
 # High-resolution TIFF and PDF preserve the inherited artwork and plot placement.
 ragg::agg_tiff(file.path(out,'Figure1_updated.tiff'),width=183/25.4,height=1795/ppi,units='in',res=600,compression='lzw')
 grid.newpage();grid.raster(im,interpolate=FALSE);dev.off()
 cairo_pdf(file.path(out,'Figure1_updated.pdf'),width=183/25.4,height=1795/ppi,family='Arial')
 grid.newpage();grid.raster(im,interpolate=FALSE);dev.off()
}
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
