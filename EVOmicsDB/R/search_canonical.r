#!/usr/bin/env Rscript
# Read-only presentation adapter. No hypothesis tests, model fits, or save().
suppressPackageStartupMessages(library(jsonlite))
requests <- fromJSON(paste(readLines(file("stdin"), warn = FALSE), collapse = "\n"), simplifyVector = FALSE)
finite <- function(x) unname(as.numeric(x[is.finite(x)]))
scalar <- function(x) if (length(x) == 1L && is.finite(x)) unname(as.numeric(x)) else NA_real_
summary_values <- function(x, prefix) {
  vals <- if (length(x)) c(min(x), quantile(x,.25), median(x), quantile(x,.75), max(x), mean(x), sqrt(mean((x-mean(x))^2))) else rep(NA_real_,7)
  setNames(as.list(unname(vals)), paste0(prefix,c("min","p25","median","p75","max","mean","std")))
}
rows <- lapply(requests, function(req) {
  e <- new.env(parent=emptyenv()); load(req$rda,envir=e)
  stopifnot(all(c("deg","exp_model","group_list","preprocessing") %in% ls(e)))
  m <- as.matrix(e$exp_model); g <- e$group_list
  if (!is.null(names(g))) g <- g[match(colnames(m),names(g))]
  stopifnot(length(g)==ncol(m), !anyNA(g), setequal(as.character(g),c("case","control")))
  ctrl <- g=="control"; cases <- g=="case"; gene <- req$gene_name
  v <- if (gene %in% rownames(m)) as.numeric(m[gene,]) else rep(NA_real_,ncol(m))
  nc <- sum(is.finite(v[ctrl])); nt <- sum(is.finite(v[cases]))
  cm <- mean(v[ctrl],na.rm=TRUE); cs <- sd(v[ctrl],na.rm=TRUE)
  z <- if (is.finite(cs) && cs > sqrt(.Machine$double.eps)) (v-cm)/cs else rep(NA_real_,length(v))
  cvals <- finite(z[ctrl]); tvals <- finite(z[cases])
  deg <- e$deg; d <- if (gene %in% rownames(deg)) deg[gene,,drop=FALSE] else NULL
  pv <- scalar(d$P.Value); adj <- scalar(d$adj.P.Val); fc <- scalar(d$logFC)
  reason <- if (is.null(d)) "Feature not present in canonical model results" else if(nc<2 || nt<2) "Insufficient feature-valid samples (fewer than two in a group)" else if(!is.finite(pv)) "Canonical raw P is not available" else NULL
  if (!is.null(reason)) { pv<-NA_real_; adj<-NA_real_; fc<-NA_real_ }
  pre <- e$preprocessing
  c(list(gene_name=gene,dataset_id=req$dataset_id),summary_values(cvals,"ctrl_"),summary_values(tvals,"tumor_"),
    list(ctrl_values=as.list(cvals),tumor_values=as.list(tvals),log2fc=fc,pvalue=pv,adjusted_pvalue=adj,
         statistical_status=if(is.null(reason)) "available" else "not_available", statistical_reason=reason,
         statistical_source=list(canonical_rda=basename(req$rda),model=pre$model_method,profile=pre$analysis_profile,
           design=pre$design_formula,contrast=pre$contrast,effect_size_type=pre$effect_size_type,
           p_field="deg.P.Value",adjusted_p_field="deg.adj.P.Val",adjustment="BH within dataset modelled features",
           log2fc_field="deg.logFC",sample_source="exp_model + group_list",display_transform="(exp_model - control mean) / control SD",
           display_status=if(length(cvals)&&length(tvals)) "available" else "not_available: control SD or expression unavailable"),
         analysis_n_control=sum(ctrl),analysis_n_case=sum(cases),feature_n_control=nc,feature_n_case=nt,
         display_n_control=length(cvals),display_n_case=length(tvals)))
})
cat(toJSON(rows,auto_unbox=TRUE,na="null",null="null",digits=16))
