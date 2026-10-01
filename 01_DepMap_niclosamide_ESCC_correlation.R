# Auditable DepMap niclosamide / ESCC expression correlations.
# Base R only. No automatic downloads: use a fingerprinted input snapshot.
# Set DEPMAP_INPUT_DIR and DEPMAP_OUTPUT_DIR; see AUDIT_2026-10-01.md.

if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', '.UTF-8'))

require_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, ' missing columns: ', paste(missing, collapse=', '), call.=FALSE)
}
assert_unique <- function(x, keys, label) {
  if (anyNA(x[keys]) || anyDuplicated(x[keys]))
    stop(label, ' has missing or duplicate keys: ', paste(keys, collapse=', '), call.=FALSE)
}
run_cor_tests <- function(dependency, expression) {
  if (length(dependency) != length(expression)) stop('Correlation vectors must have equal lengths.')
  if (!is.numeric(dependency) || !is.numeric(expression)) stop('Correlation inputs must be numeric.')
  ok <- is.finite(dependency) & is.finite(expression)
  x <- dependency[ok]; y <- expression[ok]
  result <- c(n=length(x), Spearman_Correlation=NA_real_, Spearman_P_Value=NA_real_,
              Pearson_Correlation=NA_real_, Pearson_P_Value=NA_real_)
  if (length(x)<3L || sd(x)==0 || sd(y)==0) return(result)
  sp <- suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE))
  pe <- cor.test(x,y,method='pearson')
  result[c('Spearman_Correlation','Spearman_P_Value','Pearson_Correlation','Pearson_P_Value')] <-
    c(unname(sp$estimate),sp$p.value,unname(pe$estimate),pe$p.value)
  result
}
read_pinned_inputs <- function(input_dir) {
  manifest_path <- file.path(input_dir,'input_manifest.csv')
  if (!file.exists(manifest_path)) stop('Missing input_manifest.csv; unverified or live inputs are not accepted.')
  manifest <- read.csv(manifest_path,stringsAsFactors=FALSE,check.names=FALSE)
  require_columns(manifest,c('file','md5','source'),'Input manifest')
  assert_unique(manifest,'file','Input manifest')
  objects <- list()
  for (name in c('metadata','drug_sensitivity','TPM_data')) {
    candidates <- paste0(name,c('.csv','.rds'))
    chosen <- candidates[file.exists(file.path(input_dir,candidates))]
    if (length(chosen)!=1L) stop('Provide exactly one CSV or RDS for ',name)
    row <- manifest[manifest$file==chosen,,drop=FALSE]
    if(nrow(row)!=1L || is.na(row$source) || !nzchar(row$source)) stop('Missing provenance for ',chosen)
    path <- file.path(input_dir,chosen)
    if (tolower(unname(tools::md5sum(path))) != tolower(row$md5)) stop('Input checksum mismatch: ',chosen)
    objects[[name]] <- if(endsWith(chosen,'.csv')) read.csv(path,stringsAsFactors=FALSE,check.names=FALSE) else as.data.frame(readRDS(path))
  }
  objects$manifest <- manifest
  objects
}
prepare_cohort <- function(metadata, drug, expression, mode=Sys.getenv('DEPMAP_COHORT_MODE','curated_escc')) {
  if (!mode %in% c('curated_escc','historical_metadata_only')) stop('Unknown cohort mode.')
  require_columns(metadata,c('depmap_id','primary_disease','lineage_subtype','cell_line_name'),'Metadata')
  require_columns(drug,c('depmap_id','compound','dependency'),'Drug response')
  # Entrez IDs are not unique in the archived release (17 annotation collisions).
  # The original DepMap 'gene' profile identifier is the analysis key.
  if ('gene' %in% names(expression)) {
    if ('gene_id' %in% names(expression) && !identical(as.character(expression$gene_id),as.character(expression$gene))) stop('gene_id must preserve the original gene profile identifier; Entrez IDs can collide.')
    expression$gene_id <- as.character(expression$gene)
  }
  require_columns(expression,c('depmap_id','gene_name','gene_id','rna_expression'),'Expression')
  assert_unique(metadata,'depmap_id','Metadata')
  if (!is.numeric(drug$dependency) || !is.numeric(expression$rna_expression)) stop('Non-numeric response/expression input.')
  meta <- metadata[!is.na(metadata$primary_disease) & !is.na(metadata$lineage_subtype) &
    metadata$primary_disease=='Esophageal Cancer' & metadata$lineage_subtype=='esophagus_squamous',,drop=FALSE]
  if (!nrow(meta)) stop('No eligible ESCC metadata rows.')
  require_columns(meta, 'Cellosaurus_NCIt_id', 'Histology cross-check')
  histology_ok <- !is.na(meta$Cellosaurus_NCIt_id) & meta$Cellosaurus_NCIt_id == 'C4024'
  histology_audit <- meta[c('depmap_id','cell_line_name','lineage_subtype','Cellosaurus_NCIt_id')]
  histology_audit$histology_concordant <- histology_ok
  histology_audit$included_by_histology <- if (mode=='curated_escc') histology_ok else TRUE
  histology_audit$reason <- ifelse(histology_ok, '', 'Conflicting or missing Cellosaurus histology; excluded in curated mode')
  if (mode=='curated_escc') meta <- meta[histology_ok,,drop=FALSE]
  if (mode=='historical_metadata_only') warning('Historical replay includes histology conflicts; do not report as a curated ESCC cohort.')
  selected <- drug[!is.na(drug$compound) & drug$compound=='BRD-K35960502-001-20-0::2.5::HTS' & drug$depmap_id %in% meta$depmap_id,,drop=FALSE]
  if (!nrow(selected)) stop('Exact niclosamide compound is absent from the eligible cohort.')
  assert_unique(selected,'depmap_id','Selected niclosamide response; do not aggregate doses or replicates implicitly')
  selected <- selected[order(selected$depmap_id),c('depmap_id','compound','dependency'),drop=FALSE]
  selected$cell_line_name <- meta$cell_line_name[match(selected$depmap_id,meta$depmap_id)]
  if (anyNA(selected$cell_line_name) || any(!nzchar(selected$cell_line_name))) stop('Missing cell-line labels.')
  expression <- expression[expression$depmap_id %in% selected$depmap_id,,drop=FALSE]
  assert_unique(expression,c('depmap_id','gene_id'),'Expression')
  if (anyNA(expression$gene_name) || any(!nzchar(expression$gene_name))) stop('Missing gene symbols.')
  gene_map <- unique(expression[c('gene_id','gene_name')])
  assert_unique(gene_map,'gene_id','Gene ID to symbol mapping')
  assert_unique(gene_map,'gene_name','Gene symbols; resolve ambiguous IDs before analysis')
  finite_ids <- unique(expression$depmap_id[is.finite(expression$rna_expression)])
  selected$included_in_expression_analysis <- is.finite(selected$dependency) & selected$depmap_id %in% finite_ids
  selected$exclusion_reason <- ifelse(!is.finite(selected$dependency),'nonfinite drug response',
    ifelse(!selected$depmap_id %in% finite_ids,'no finite expression profile',''))
  cohort <- selected[selected$included_in_expression_analysis,,drop=FALSE]
  if(nrow(cohort)<3L) stop('Fewer than three matched cells.')
  joined <- merge(expression,cohort[c('depmap_id','dependency','cell_line_name')],by='depmap_id',sort=TRUE)
  assert_unique(joined,c('depmap_id','gene_id'),'Joined analysis')
  list(all_drug=selected,cohort=cohort,joined=joined,histology_audit=histology_audit,mode=mode)
}
main <- function() {
  input_dir <- Sys.getenv('DEPMAP_INPUT_DIR')
  output_dir <- Sys.getenv('DEPMAP_OUTPUT_DIR',file.path('Data','DepMap'))
  if (!nzchar(input_dir)) stop('Set DEPMAP_INPUT_DIR to the fingerprinted CSV/RDS input directory.')
  if (dir.exists(output_dir) && length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
    stop('Use a new empty DEPMAP_OUTPUT_DIR; previous outputs will not be overwritten.')
  inputs <- read_pinned_inputs(input_dir)
  cohort <- prepare_cohort(inputs$metadata,inputs$drug_sensitivity,inputs$TPM_data)
  joined <- cohort$joined
  by_gene <- split(seq_len(nrow(joined)),joined$gene_id)
  results <- do.call(rbind,lapply(by_gene,function(ii) {
    x <- joined[ii,,drop=FALSE]
    data.frame(gene_name=x$gene_name[1],gene_id=x$gene_id[1],as.list(run_cor_tests(x$dependency,x$rna_expression)),stringsAsFactors=FALSE)
  }))
  rownames(results) <- NULL
  results <- results[order(results$gene_name),,drop=FALSE]
  results$Spearman_FDR <- p.adjust(results$Spearman_P_Value,method='BH')
  results$Pearson_FDR <- p.adjust(results$Pearson_P_Value,method='BH')
  results$selection_rule <- 'both_nominal_p_lt_0.05'
  primary <- with(results,is.finite(Pearson_P_Value)&is.finite(Spearman_P_Value)&Pearson_P_Value<0.05&Spearman_P_Value<0.05)
  significant <- results[primary,,drop=FALSE]
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  figure_dir <- file.path(output_dir,'Figures');dir.create(figure_dir)
  write_out <- function(x,name) write.csv(x,file.path(output_dir,name),row.names=FALSE,na='')
  write_out(inputs$manifest,'input_manifest.csv')
  write_out(cohort$histology_audit,'histology_audit.csv')
  writeLines(paste0('cohort_mode=',cohort$mode),file.path(output_dir,'cohort_mode.txt'))
  write_out(cohort$all_drug,'drug_cohort_and_exclusions.csv')
  write_out(cohort$cohort,'analysis_cell_lines.csv')
  write_out(results,'Niclosamide_ESCC_allgene_correlations.csv')
  write_out(significant,'Niclosamide_ESCC_allgene_significant_results.csv')
  pearson_only <- results[is.finite(results$Pearson_P_Value)&results$Pearson_P_Value<0.05,,drop=FALSE]
  write_out(pearson_only,'Pearson_only_nominal_sensitivity.csv')
  summary <- data.frame(metric=c('drug_available_cells','matched_expression_cells','expression_profiles','testable_gene_profiles','untestable_gene_profiles',
      'dual_nominal_selected','dual_negative_pearson','dual_positive_pearson','pearson_only_nominal',
      'Pearson_BH_lt_0.05','Spearman_BH_lt_0.05'),
    value=c(nrow(cohort$all_drug),nrow(cohort$cohort),nrow(results),sum(is.finite(results$Pearson_P_Value)),sum(!is.finite(results$Pearson_P_Value)),nrow(significant),
      sum(significant$Pearson_Correlation<0),sum(significant$Pearson_Correlation>0),nrow(pearson_only),
      sum(results$Pearson_FDR<0.05,na.rm=TRUE),sum(results$Spearman_FDR<0.05,na.rm=TRUE)))
  write_out(summary,'analysis_summary.csv')
  # Every figure uses the same matched cell IDs, and scatter statistics use the same finite pairs.
  response_label <- 'PRISM log2 viability fold change (niclosamide 2.5 uM)'
  d <- cohort$cohort[order(cohort$cohort$dependency),]
  png(file.path(figure_dir,'Niclosamide_ESCC_response_barplot.png'),width=2100,height=1800,res=300)
  par(mar=c(5,8,3,1));barplot(d$dependency,names.arg=d$cell_line_name,horiz=TRUE,las=1,col='lightblue',
    xlab=response_label,main=paste(if(cohort$mode=='curated_escc') 'Curated ESCC cohort: n =' else 'Historical cohort (includes histology conflict): n =',nrow(d)),cex.main=0.95,cex.names=0.7)
  abline(v=0,lty=2);dev.off()
  for (method in c('Pearson','Spearman')) {
    r <- results[[paste0(method,'_Correlation')]];p <- results[[paste0(method,'_P_Value')]]
    ok <- is.finite(r)&is.finite(p)
    png(file.path(figure_dir,paste0(method,'_all_gene_correlations.png')),width=2100,height=1800,res=300)
    plot(r[ok],-log10(pmax(p[ok],.Machine$double.xmin)),pch=16,cex=0.45,
      col=ifelse(primary[ok],'#B23B3B','#A0A0A0'),xlab=paste(method,'correlation'),ylab=paste0('-log10(',method,' p value)'),
      main='All tested genes; red = both nominal p < 0.05')
    abline(h=-log10(0.05),lty=2);abline(v=0,lty=3);dev.off()
  }
  for (gene in c('EHMT2','MEN1')) {
    d <- joined[joined$gene_name==gene & is.finite(joined$dependency)&is.finite(joined$rna_expression),,drop=FALSE]
    rr <- results[results$gene_name==gene,,drop=FALSE]
    if(nrow(rr)!=1L || nrow(d)<3L || !is.finite(rr$Pearson_Correlation)) stop('Candidate gene missing or untestable: ',gene)
    check <- run_cor_tests(d$dependency,d$rna_expression)
    stopifnot(isTRUE(all.equal(unname(check['Pearson_Correlation']),rr$Pearson_Correlation,tolerance=1e-12)),nrow(d)==rr$n)
    write_out(d,paste0(gene,'_plotted_pairs.csv'))
    png(file.path(figure_dir,paste0(gene,'_response_expression_scatter.png')),width=2100,height=1800,res=300)
    par(mar=c(6,5,4,2));plot(d$dependency,d$rna_expression,pch=19,
      ylim=range(d$rna_expression)+c(-0.07,0.10)*diff(range(d$rna_expression)),
      xlab=response_label,ylab='Baseline RNA expression log2(TPM + 1)',main=gene,
      sub=sprintf('n = %d; Pearson r = %.3f, p = %.4g, BH q = %.5g',rr$n,rr$Pearson_Correlation,rr$Pearson_P_Value,rr$Pearson_FDR))
    abline(lm(rna_expression~dependency,data=d),col='red')
    # Exact cell labels and plotted coordinates are exported in the paired CSV.
    dev.off()
  }
  write_out(results[results$gene_name %in% c('EHMT2','MEN1'),],'candidate_gene_statistics.csv')
  writeLines(c('The historical selection requires BOTH Pearson and Spearman nominal p < 0.05.',
    'Pearson-only selection is exported separately and must not be called the historical selected set.',
    'BH q values are reported across all tested genes; nominal selection is exploratory.',
    'A negative correlation means higher baseline expression is associated with a lower response score.',
    'This observational association does not establish causality or a validated response biomarker.'),file.path(output_dir,'interpretation.txt'))
  writeLines(capture.output(sessionInfo()),file.path(output_dir,'sessionInfo_DepMap.txt'))
  print(summary)
}
if (Sys.getenv('DEPMAP_TEST_MODE')!='1') main()
