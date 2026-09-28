# Stable pair keys must not depend on the host LC_COLLATE setting.
kegg_joint_pair_key <- function(type1,type2) {
  types<-c(type1,type2)
  priority<-c('Metabolite','mRNA','miRNA','lncRNA','Protein')
  if(any(!types %in% priority) || anyDuplicated(types)) return('unsupported')
  paste(types[order(match(types,priority))],collapse='_')
}
# Assay-controlled KEGG ORA for metabolite-containing pairs. No live requests.
# Both annotation layers are frozen together and use KEGG mapNNNNN IDs.
kegg_joint_load_snapshot <- function(path) {
  if (is.null(path) || !dir.exists(path)) stop('A local KEGG joint snapshot directory is required')
  manifest <- jsonlite::fromJSON(file.path(path, 'manifest.json'))
  for (name in names(manifest$files_sha256)) {
    file <- file.path(path, name)
    if (!file.exists(file)) stop('Snapshot file missing: ', name,
      '. KEGG reference tables are external licensed inputs, excluded from the public code package. ',
      'Restore the version-matched local snapshot described in docs/REPRODUCING.md#external-annotations; SHA256 checks must pass.')
    actual <- digest::digest(file=file, algo='sha256')
    if (!identical(actual, manifest$files_sha256[[name]])) stop('Snapshot SHA256 mismatch: ', name)
  }
  read <- function(name) utils::read.delim(file.path(path,name), check.names=FALSE, stringsAsFactors=FALSE)
  list(manifest=manifest, id=manifest$snapshot_id,
       fingerprint=digest::digest(file=file.path(path,'manifest.json'),algo='sha256'),
       pathways=read('pathways.tsv'), gene=read('gene_pathway.tsv'), compound=read('compound_pathway.tsv'),
       symbols=read('gene_symbols.tsv'), synonyms=read('compound_synonyms.tsv'))
}
kegg_joint_name_key <- function(x) tolower(gsub('[[:space:]]+', ' ', trimws(as.character(x))))
kegg_joint_map_compounds <- function(ids, snapshot, curated_path=NULL) {
  ids <- unique(as.character(ids)); synonyms <- snapshot$synonyms
  by_name <- split(as.character(synonyms$compound_id), kegg_joint_name_key(synonyms$synonym))
  curated <- NULL
  if (!is.null(curated_path) && nzchar(curated_path)) {
    if (!file.exists(curated_path)) stop('Curated metabolite mapping file not found: ',curated_path)
    curated <- utils::read.table(curated_path,sep=if(grepl('\\.csv$',curated_path,ignore.case=TRUE)) ',' else '\t', header=TRUE,quote='"',comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
    if (!all(c('Metabolite','KEGG_ID') %in% names(curated))) stop('Curated metabolite mapping requires Metabolite and KEGG_ID')
  }
  known <- unique(as.character(synonyms$compound_id))
  if (!length(ids)) return(data.frame(feature_id=character(),mapped_id=character(),mapping_status=character(),candidates=character()))
  do.call(rbind,lapply(ids,function(id) {
    matches <- character(); method <- 'exact_unique'
    if (!is.null(curated)) {
      hits <- curated[curated$Metabolite==id & !is.na(curated$KEGG_ID),,drop=FALSE]
      if (nrow(hits)) {matches <- unique(sub('^cpd:','',as.character(hits$KEGG_ID)));method <- 'curated_unique'}
    }
    if (!length(matches) && grepl('^(cpd:)?C[0-9]{5}$',id)) {matches<-sub('^cpd:','',id);method<-'explicit_KEGG_ID'}
    if (!length(matches)) matches <- unique(by_name[[kegg_joint_name_key(id)]])
    status <- if(!length(matches)) 'unmapped' else if(length(matches)>1) 'ambiguous' else if(!matches %in% known) 'invalid_KEGG_ID' else method
    mapped <- if(status %in% c('exact_unique','curated_unique','explicit_KEGG_ID')) matches else NA_character_
    data.frame(feature_id=id,mapped_id=mapped,mapping_status=status,candidates=paste(matches,collapse=';'),stringsAsFactors=FALSE)
  }))
}
kegg_joint_map_genes <- function(ids, snapshot) {
  ids <- unique(trimws(as.character(ids))); ann<-snapshot$symbols
  by_name <- split(as.character(ann$gene_id),as.character(ann$gene_symbol))
  if (!length(ids)) return(data.frame(feature_id=character(),mapped_id=character(),mapping_status=character(),candidates=character()))
  do.call(rbind,lapply(ids,function(id) {
    hit<-unique(by_name[[toupper(id)]])
    status<-if(!length(hit)) 'unmapped' else if(length(hit)>1) 'ambiguous' else 'exact_symbol'
    data.frame(feature_id=id,mapped_id=if(length(hit)==1) hit else NA_character_,mapping_status=status,candidates=paste(hit,collapse=';'),stringsAsFactors=FALSE)
  }))
}
kegg_joint_ora <- function(foreground,universe,term2feature,pathways,p_adjust_method='BH',min_size=10L,max_size=500L) {
  if (!p_adjust_method %in% c('BH','fdr','bonferroni')) stop('Unsupported p_adjust_method')
  mapped <- unique(as.character(term2feature[[2]]))
  universe <- intersect(unique(as.character(universe)),mapped)
  foreground <- intersect(unique(as.character(foreground)),universe)
  sets<-split(as.character(term2feature[[2]]),term2feature[[1]])
  sets<-lapply(sets,function(x) intersect(unique(x),universe))
  sets<-sets[lengths(sets)>=min_size & lengths(sets)<=max_size]
  out<-data.frame(ID=character(),Description=character(),Count=integer(),GeneRatio=character(),BgRatio=character(),pvalue=numeric(),p.adjust=numeric(),GeneID=character(),universe_n=integer(),foreground_n=integer(),pathway_universe_n=integer(),stringsAsFactors=FALSE)
  if (!length(sets) || !length(universe)) return(out)
  N<-length(universe);n<-length(foreground)
  out<-do.call(rbind,lapply(names(sets),function(id) {
    members<-sets[[id]];hit<-intersect(members,foreground);M<-length(members);k<-length(hit)
    data.frame(ID=id,Description=pathways$description[match(id,pathways$pathway_id)],Count=k,
      GeneRatio=paste0(k,'/',n),BgRatio=paste0(M,'/',N),pvalue=if(n==0) 1 else stats::phyper(k-1,M,N-M,n,lower.tail=FALSE),
      GeneID=paste(sort(hit),collapse='/'),universe_n=N,foreground_n=n,pathway_universe_n=M,stringsAsFactors=FALSE)
  }))
  # A pre-specified complete family, including zero-hit pathways (P=1).
  out$p.adjust<-stats::p.adjust(out$pvalue,method=p_adjust_method);out$tested_family_n<-nrow(out);out$p_adjust_method<-p_adjust_method
  out
}
kegg_joint_combine <- function(gene,metabolite,gene_label,pathways,threshold=.05,p_adjust_method='BH') {
  names(gene)[names(gene)!='ID']<-paste0(names(gene)[names(gene)!='ID'],'_Gene')
  names(metabolite)[names(metabolite)!='ID']<-paste0(names(metabolite)[names(metabolite)!='ID'],'_Metab')
  x<-merge(gene,metabolite,by='ID',all=TRUE,sort=TRUE)
  if(!nrow(x)) return(x)
  x$Description<-pathways$description[match(x$ID,pathways$pathway_id)]
  x$tested_Gene<-is.finite(x$pvalue_Gene);x$tested_Metab<-is.finite(x$pvalue_Metab)
  x$tested_in_both<-x$tested_Gene & x$tested_Metab
  x$significant_Gene<-is.finite(x$p.adjust_Gene) & x$p.adjust_Gene<=threshold
  x$significant_Metab<-is.finite(x$p.adjust_Metab) & x$p.adjust_Metab<=threshold
  x$shared_significant<-x$significant_Gene & x$significant_Metab
  x$SourceType<-ifelse(x$shared_significant,'Common',ifelse(x$significant_Gene,gene_label,ifelse(x$significant_Metab,'Metabolite','Not significant')))
  x$joint_p<-NA_real_;x$joint_FDR<-NA_real_;i<-which(x$tested_in_both)
  if(length(i)) {
    x$joint_p[i]<-stats::pchisq(-2*(log(pmax(x$pvalue_Gene[i],.Machine$double.xmin))+log(pmax(x$pvalue_Metab[i],.Machine$double.xmin))),df=4,lower.tail=FALSE)
    x$joint_FDR[i]<-p.adjust(x$joint_p[i],method=p_adjust_method)
  }
  x$combined_p_exploratory<-x$joint_p;x$combined_FDR_exploratory<-x$joint_FDR
  x$joint_tested_family_n<-length(i);x$joint_p_adjust_method<-p_adjust_method
  x$joint_statistic_note<-'Exploratory Fisher combination assumes independent null P values; cross-omics dependence may violate calibration. Not used to label shared significance.'
  x$p.adjust<-ifelse(x$shared_significant,pmax(x$p.adjust_Gene,x$p.adjust_Metab),ifelse(x$significant_Gene,x$p.adjust_Gene,ifelse(x$significant_Metab,x$p.adjust_Metab,NA_real_)))
  x$p.adjust_role<-'Display only: assay adjusted P; Common uses maximum of both assay adjusted P values'
  x$Count<-ifelse(x$shared_significant,x$Count_Gene+x$Count_Metab,ifelse(x$significant_Gene,x$Count_Gene,ifelse(x$significant_Metab,x$Count_Metab,0)))
  x$direction_test_contract<-'One non-directional union foreground ORA per assay; Up/Down only annotations'
  x
}
kegg_joint_tested_ids <- function(df, universe) {
  ids <- as.character(df$ID[is.finite(df$P.Value)])
  if (is.null(universe) || !length(universe)) stop('An explicit tested feature universe is required')
  universe <- if(length(universe)==1L && file.exists(universe)) readLines(universe,warn=FALSE) else as.character(universe)
  unique(intersect(ids,trimws(universe)))
}
# Return the same gene union as an expanded lncRNA x miRNA x gene join,
# without materializing every path. The two annotation edges remain auditable.
kegg_joint_lnc_reachability <- function(links, associations, selected_lnc) {
  links<-unique(links[,c('lncRNA','target'),drop=FALSE])
  links$miRNA_key<-normalize_mirna_exact(links$target)
  links$selected_lncRNA<-links$lncRNA %in% selected_lnc
  keys<-normalize_mirna_exact(associations$mature_mirna_id)
  keep<-keys %in% links$miRNA_key
  associations<-associations[keep,,drop=FALSE];keys<-keys[keep]
  selected_keys<-unique(links$miRNA_key[links$selected_lncRNA])
  out<-unique(data.frame(source_feature=as.character(associations$mature_mirna_id),
    target_gene=as.character(associations$target_symbol),
    projection_step=rep('linked_miRNA_to_gene',nrow(associations)),is_measured_source=rep(FALSE,nrow(associations)),
    selected=keys %in% selected_keys,stringsAsFactors=FALSE))
  list(targets=out,links=links)
}
kegg_joint_regulator_targets <- function(ids,type,master_lnc_targets=NULL,selected_ids=character()) {
  if(type=='miRNA') {
    if(is.null(args$mirna_resolution) || !nzchar(args$mirna_resolution)) stop('Explicit miRNA measurement resolution required')
    ann<-get_mirna_targets(data.frame(ID=ids),max_targets_per_mirna=NULL)
    return(unique(data.frame(source_feature=as.character(ann$mature_mirna_id),target_gene=as.character(ann$target_symbol))))
  }
  ann<-master_lnc_targets[master_lnc_targets$lncRNA %in% ids,,drop=FALSE]
  direct<-ann[!grepl('miRNA',ann$target_type,ignore.case=TRUE) & grepl('protein|pcg|gene|tf|RNA',ann$target_type,ignore.case=TRUE),,drop=FALSE]
  out<-unique(data.frame(source_feature=as.character(direct$lncRNA),target_gene=as.character(direct$target),
      projection_step=rep('direct_lncRNA_to_gene',nrow(direct)),is_measured_source=rep(TRUE,nrow(direct)),
      selected=direct$lncRNA %in% selected_ids,stringsAsFactors=FALSE))
  intermediate<-unique(ann[grepl('miRNA',ann$target_type,ignore.case=TRUE),c('lncRNA','target'),drop=FALSE])
  cat('  lncRNA direct annotation rows:',nrow(out),' intermediary links:',nrow(intermediate),'\n')
  links<-intermediate
  if(nrow(intermediate)) {
    saved_resolution<-args$mirna_resolution;args$mirna_resolution<<-'mirbase_mature'
    on.exit({args$mirna_resolution<<-saved_resolution},add=TRUE)
    projected<-get_mirna_targets(data.frame(ID=unique(intermediate$target)),max_targets_per_mirna=NULL)
    reachable<-kegg_joint_lnc_reachability(intermediate,projected,selected_ids)
    out<-unique(rbind(out,reachable$targets));links<-reachable$links
  }
  attr(out,'lncrna_mirna_links')<-links
  out
}
build_kegg_metabolite_pair <- function(partner_df,metabolite_df,partner_type,partner_universe,metabolite_universe,master_lnc_targets=NULL) {
  if(!identical(args$organism,'hsa')) stop('The packaged joint KEGG snapshot is human (hsa) only')
  if(is.null(args$kegg_joint_snapshot_dir)) stop('Missing frozen KEGG joint snapshot')
  snapshot<-kegg_joint_load_snapshot(args$kegg_joint_snapshot_dir)
  tested_partner<-kegg_joint_tested_ids(partner_df,partner_universe)
  tested_metabolite<-kegg_joint_tested_ids(metabolite_df,metabolite_universe)
  pcol<-if(isTRUE(args$use_padj_for_features)) 'adj.P.Val' else 'P.Value'
  selected_partner<-intersect(tested_partner,as.character(partner_df$ID[is.finite(partner_df[[pcol]]) & partner_df[[pcol]]<=feature_pvalue_threshold]))
  selected_metabolite<-intersect(tested_metabolite,as.character(metabolite_df$ID[is.finite(metabolite_df[[pcol]]) & metabolite_df[[pcol]]<=feature_pvalue_threshold]))
  compound_audit<-kegg_joint_map_compounds(tested_metabolite,snapshot,args$metabolite_mapping_file)
  compound_audit$selected<-compound_audit$feature_id %in% selected_metabolite
  compound_audit$in_pathway_annotation<-!is.na(compound_audit$mapped_id) & compound_audit$mapped_id %in% snapshot$compound$compound_id
  compound_audit$logFC<-metabolite_df$logFC[match(compound_audit$feature_id,metabolite_df$ID)]
  cpd_universe<-unique(na.omit(compound_audit$mapped_id));cpd_fg<-unique(na.omit(compound_audit$mapped_id[compound_audit$selected]))
  target_audit<-NULL
  if(partner_type %in% c('mRNA','Protein')) {
    gene_audit<-kegg_joint_map_genes(tested_partner,snapshot)
    gene_audit$selected<-gene_audit$feature_id %in% selected_partner
    label<-if(partner_type=='mRNA') 'RNA' else 'Protein'
  } else {
    if(is.null(args$mirna_target_db) || !file.exists(args$mirna_target_db)) stop('A frozen local miRNA target database is required for regulatory partner projection')
    target_audit<-kegg_joint_regulator_targets(tested_partner,partner_type,master_lnc_targets,selected_partner)
    selected_keys<-if(partner_type=='miRNA') resolve_mirna_feature_key(selected_partner,args$mirna_resolution) else selected_partner
    target_keys<-if(partner_type=='miRNA') resolve_mirna_feature_key(target_audit$source_feature,args$mirna_resolution) else target_audit$source_feature
    if(!"selected" %in% names(target_audit)) target_audit$selected<-target_keys %in% selected_keys
    gene_audit<-kegg_joint_map_genes(unique(target_audit$target_gene),snapshot)
    gene_audit$selected<-gene_audit$feature_id %in% target_audit$target_gene[target_audit$selected]
    label<-if(partner_type=='miRNA') 'miRNA_target' else 'lncRelatedGene'
  }
  gene_audit$in_pathway_annotation<-!is.na(gene_audit$mapped_id) & gene_audit$mapped_id %in% snapshot$gene$gene_id
  gene_universe<-unique(na.omit(gene_audit$mapped_id));gene_fg<-unique(na.omit(gene_audit$mapped_id[gene_audit$selected]))
  gene_result<-kegg_joint_ora(gene_fg,gene_universe,snapshot$gene,snapshot$pathways,args$p_adjust_method)
  compound_result<-kegg_joint_ora(cpd_fg,cpd_universe,snapshot$compound,snapshot$pathways,args$p_adjust_method)
  # Direction is annotation only. Each unique compound is counted once per ORA.
  up<-unique(na.omit(compound_audit$mapped_id[compound_audit$selected & is.finite(compound_audit$logFC) & compound_audit$logFC>0]))
  down<-unique(na.omit(compound_audit$mapped_id[compound_audit$selected & is.finite(compound_audit$logFC) & compound_audit$logFC<0]))
  compound_result$Count_up<-vapply(strsplit(compound_result$GeneID,'/',fixed=TRUE),function(x) as.integer(length(intersect(x,up))),integer(1))
  compound_result$Count_down<-vapply(strsplit(compound_result$GeneID,'/',fixed=TRUE),function(x) as.integer(length(intersect(x,down))),integer(1))
  compound_result$Count_conflicting_direction<-vapply(strsplit(compound_result$GeneID,'/',fixed=TRUE),function(x) as.integer(length(intersect(x,intersect(up,down)))),integer(1))
  prefix<-sub('\\.[^.]+$','',args$output_csv)
  write_audit<-function(x,suffix) utils::write.table(x,paste0(prefix,suffix),sep='\t',row.names=FALSE,quote=FALSE,na='')
  write_audit(compound_audit,'_compound_mapping.tsv');write_audit(gene_audit,'_gene_mapping.tsv')
  if(!is.null(target_audit)) {
    write_audit(target_audit,'_regulator_targets.tsv')
    if(!is.null(attr(target_audit,'lncrna_mirna_links'))) write_audit(attr(target_audit,'lncrna_mirna_links'),'_lncrna_mirna_links.tsv')
  }
  write_audit(gene_result,'_gene_ORA.tsv');write_audit(compound_result,'_compound_ORA.tsv')
  hashes<-function(path) if(!is.null(path) && file.exists(path)) digest::digest(file=path,algo='sha256') else NULL
  counts<-list(tested_partner_features=length(tested_partner),selected_partner_features=length(selected_partner),tested_metabolite_features=length(tested_metabolite),selected_metabolite_features=length(selected_metabolite),mapped_gene_universe=length(intersect(gene_universe,snapshot$gene$gene_id)),mapped_gene_foreground=length(intersect(gene_fg,snapshot$gene$gene_id)),mapped_compound_universe=length(intersect(cpd_universe,snapshot$compound$compound_id)),mapped_compound_foreground=length(intersect(cpd_fg,snapshot$compound$compound_id)),gene_tested_pathways=nrow(gene_result),compound_tested_pathways=nrow(compound_result))
  audit<-list(snapshot_id=snapshot$id,snapshot_sha256=snapshot$fingerprint,counts=counts,p_adjust_method=args$p_adjust_method,mirna_resolution=args$mirna_resolution,mirna_evidence=args$mirna_evidence,analytical_target_cap=NULL,legacy_display_max_targets_per_mirna=args$max_targets_per_mirna,mirna_target_db_sha256=hashes(args$mirna_target_db),lncrna_target_db_sha256=hashes(args$lncrna_target_db),metabolite_mapping_sha256=hashes(args$metabolite_mapping_file),min_pathway_size=10,max_pathway_size=500,universe_contract='Finite-P canonical tested features, restricted to explicit universe, then mapped with identical foreground/background rules',hypothesis_family='All background-eligible pathways including zero-hit P=1; adjustment separately per assay (BH default)',direction_contract='One ORA on union of selected compounds; direction annotates only',joint_note='Fisher/BH is exploratory because null-P independence is not guaranteed; network shared classification requires both assay FDRs')
  jsonlite::write_json(audit,paste0(prefix,'_annotation_audit.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
  combined<-kegg_joint_combine(gene_result,compound_result,label,snapshot$pathways,pathway_fdr_threshold,args$p_adjust_method)
  if(!nrow(combined)) {
    utils::write.table(combined,args$output_csv,sep='\t',row.names=FALSE,quote=FALSE)
    kegg_joint_empty_artifacts(args$output_html,args$output_png,args$output_pdf,
      title=if(nzchar(args$main_title)) args$main_title else 'KEGG Pathway Network',
      message='No pathway has 10–500 mapped tested members in either assay.')
    kegg_empty('no_eligible_pathways','No pathway has 10–500 mapped tested members in either assay; mapping and ORA audits were exported.',counts=counts,artifacts=list(csv=args$output_csv,html=args$output_html,png=args$output_png,pdf=args$output_pdf))
  }
  for (key in names(counts)) combined[[paste0('audit_',key)]] <- counts[[key]]
  combined$snapshot_id<-snapshot$id;combined$snapshot_sha256<-snapshot$fingerprint
  combined$assay_gene_label<-label
  cat('  KEGG joint snapshot:',snapshot$id,'\n');print(counts)
  combined
}
kegg_joint_empty_artifacts <- function(output_html,output_png=NULL,output_pdf=NULL,
                                      title='KEGG Pathway Network',message='No pathways passed the selected assay FDR cutoff.') {
  # An explicitly labelled empty result is a valid output, never a synthetic network.
  html<-paste0('<!DOCTYPE html><html><head><meta charset="utf-8"><title>',htmltools::htmlEscape(title),
    '</title></head><body style="background:white;font-family:Arial,sans-serif"><h3>',htmltools::htmlEscape(title),
    '</h3><p>',htmltools::htmlEscape(message),'</p><p>Full ORA tables and mapping audits accompany this result.</p></body></html>')
  writeLines(html,output_html,useBytes=TRUE)
  if(!is.null(output_png) || !is.null(output_pdf)) {
    p<-ggplot2::ggplot()+ggplot2::annotate('text',x=0,y=0,label=message,size=4)+
      ggplot2::labs(title=title,subtitle='Full ORA tables and mapping audits accompany this result.')+
      ggplot2::theme_void(base_size=11)+ggplot2::theme(plot.title=ggplot2::element_text(face='bold'),plot.margin=ggplot2::margin(12,16,12,16),plot.background=ggplot2::element_rect(fill='white',colour=NA))
    if(!is.null(output_png)) ggplot2::ggsave(output_png,p,width=10,height=7,dpi=300,bg='white')
    if(!is.null(output_pdf)) ggplot2::ggsave(output_pdf,p,width=10,height=7,device=grDevices::cairo_pdf,bg='white')
  }
}
