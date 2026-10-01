# Independent base-R validation of the DepMap cohort and correlation contracts.
# Run: Rscript tests/depmap_validation.R from the repository root.
Sys.setenv(DEPMAP_TEST_MODE = '1', DEPMAP_COHORT_MODE = 'curated_escc')
source('01_DepMap_niclosamide_ESCC_correlation.R')
assert_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = function(e) e)
  stopifnot(inherits(err, 'error'), grepl(pattern, conditionMessage(err)))
}
assert_close <- function(a,b) stopifnot(isTRUE(all.equal(a,b,tolerance=1e-12,check.attributes=FALSE)))
ids <- paste0('ACH-',sprintf('%06d',1:4))
meta <- data.frame(depmap_id=ids,primary_disease='Esophageal Cancer',
  lineage_subtype='esophagus_squamous',cell_line_name=paste0('CELL',1:4),Cellosaurus_NCIt_id='C4024')
drug <- data.frame(depmap_id=ids,compound='BRD-K35960502-001-20-0::2.5::HTS',dependency=c(-3,-2,-1,0))
expr <- data.frame(depmap_id=rep(ids,2),gene_name=rep(c('EHMT2','MEN1'),each=4),
  gene_id=rep(c('EHMT2_ID','MEN1_ID'),each=4),rna_expression=c(4,3,2,1,1,2,3,4))
cohort <- prepare_cohort(meta,drug,expr)
stopifnot(nrow(cohort$cohort)==4,nrow(cohort$joined)==8)
rneg <- run_cor_tests(drug$dependency,4:1)
rpos <- run_cor_tests(drug$dependency,1:4)
stopifnot(rneg['Pearson_Correlation']<0,rneg['Spearman_Correlation']<0,
          rpos['Pearson_Correlation']>0,rpos['Spearman_Correlation']>0)
assert_close(rneg['Pearson_Correlation'],-1)
# Finite pairs are removed jointly; pair count cannot include Inf/NA.
rfinite <- run_cor_tests(c(-3,Inf,-2,-1,NA),c(4,99,3,2,1))
stopifnot(rfinite['n']==3)
assert_close(rfinite['Pearson_Correlation'],-1)
stopifnot(is.na(run_cor_tests(1:4,rep(2,4))['Pearson_Correlation']),
          is.na(run_cor_tests(1:2,2:1)['Spearman_P_Value']))
assert_error(run_cor_tests(c('1','2','3'),1:3),'numeric')
assert_error(run_cor_tests(1:4,1:3),'equal lengths')
# Duplicated keys must never create pseudoreplicated correlations.
assert_error(prepare_cohort(rbind(meta,meta[1,]),drug,expr),'duplicate')
assert_error(prepare_cohort(meta,rbind(drug,drug[1,]),expr),'duplicate')
assert_error(prepare_cohort(meta,drug,rbind(expr,expr[1,])),'duplicate')
wrong <- drug; wrong$compound <- 'BRD-K35960502-001-20-0::10::HTS'
assert_error(prepare_cohort(meta,wrong,expr),'Exact niclosamide')
wrong$compound <- 'BRD-K35960502-001-20-0::2.5::OTHER'
assert_error(prepare_cohort(meta,wrong,expr),'Exact niclosamide')
# Extra doses may coexist in raw input, but must not enter the selected cohort.
extra <- drug; extra$compound <- 'BRD-K35960502-001-20-0::10::HTS';extra$dependency<-999
with_extra <- prepare_cohort(meta,rbind(drug,extra),expr)
assert_close(with_extra$cohort$dependency,cohort$cohort$dependency)
# Order of either input table cannot alter correlation estimates or cell labels.
set.seed(184)
shuffled <- prepare_cohort(meta[sample(4),],drug[sample(4),],expr[sample(8),])
stopifnot(identical(cohort$cohort$depmap_id,shuffled$cohort$depmap_id),
          identical(cohort$cohort$cell_line_name,shuffled$cohort$cell_line_name))
for(g in unique(expr$gene_id)) {
 a<-cohort$joined[cohort$joined$gene_id==g,];b<-shuffled$joined[shuffled$joined$gene_id==g,]
 assert_close(run_cor_tests(a$dependency,a$rna_expression),run_cor_tests(b$dependency,b$rna_expression))
}
# A nonfinite response is excluded before statistical and plotted cohorts form.
bad_response<-drug;bad_response$dependency[1]<-NA_real_
excluded<-prepare_cohort(meta,bad_response,expr)
stopifnot(nrow(excluded$cohort)==3,!ids[1]%in%excluded$joined$depmap_id,
          excluded$all_drug$exclusion_reason[1]=='nonfinite drug response')
cat('PASS: sign, finite pairs, duplicate keys, exact dose/screen, extra dose isolation, row-order invariance and exclusions.\n')
# Official raw TPM schema must preserve original profile keys despite Entrez collisions.
raw_expr <- expr
raw_expr$gene <- paste0(raw_expr$gene_name,' (',rep(12345,nrow(raw_expr)),')')
raw_expr$entrez_id <- rep(12345,nrow(raw_expr))
raw_expr$gene_id <- NULL
raw_cohort <- prepare_cohort(meta,drug,raw_expr)
stopifnot(length(unique(raw_cohort$joined$gene_id))==2,
          all(table(raw_cohort$joined$gene_id)==4))
wrong_key <- raw_expr;wrong_key$gene_id <- as.character(wrong_key$entrez_id)
assert_error(prepare_cohort(meta,drug,wrong_key),'original gene profile')
# Ambiguous symbol-to-profile mapping cannot silently merge distinct profiles.
ambiguous <- expr; ambiguous$gene_name <- 'EHMT2'
assert_error(prepare_cohort(meta,drug,ambiguous),'Gene symbols')
cat('PASS: official raw gene-profile keys, shared Entrez IDs, wrong-key rejection and symbol ambiguity.\n')

# Lineage tags alone cannot certify ESCC: OE19 is adenocarcinoma C4025.
conflict_meta <- meta
conflict_meta$cell_line_name[4] <- 'OE-19'
conflict_meta$Cellosaurus_NCIt_id[4] <- 'C4025'
curated <- prepare_cohort(conflict_meta,drug,expr)
stopifnot(curated$mode=='curated_escc',nrow(curated$cohort)==3,
          !ids[4]%in%curated$joined$depmap_id,
          !curated$histology_audit$included_by_histology[4],
          !curated$histology_audit$histology_concordant[4])
warned <- FALSE
historical <- withCallingHandlers(prepare_cohort(conflict_meta,drug,expr,mode='historical_metadata_only'),
  warning=function(w) { warned <<- grepl('histology conflicts',conditionMessage(w));invokeRestart('muffleWarning') })
stopifnot(warned,nrow(historical$cohort)==4,historical$mode=='historical_metadata_only',
          historical$histology_audit$included_by_histology[4],
          !historical$histology_audit$histology_concordant[4])
unknown_histology <- meta;unknown_histology$Cellosaurus_NCIt_id[4]<-NA_character_
stopifnot(nrow(prepare_cohort(unknown_histology,drug,expr)$cohort)==3)
no_histology <- meta;no_histology$Cellosaurus_NCIt_id<-NULL
assert_error(prepare_cohort(no_histology,drug,expr),'Histology cross-check')
assert_error(prepare_cohort(meta,drug,expr,mode='unknown'),'Unknown cohort mode')
cat('PASS: curated histology excludes OE19/unknown disease; explicit historical replay warns and audits conflicts.\n')
