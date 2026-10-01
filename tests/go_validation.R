Sys.setenv(DEPMAP_TEST_MODE="1")
source("02_GO_enrichR_negative_correlated_genes.R")
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(error), grepl(pattern, error))
}
x <- data.frame(gene_name = c("A", "B"), Pearson_Correlation = c(-.6, .6), Spearman_Correlation = c(-.6, .6), Pearson_P_Value = c(.01, .01), Spearman_P_Value = c(.01, .01))
stopifnot(identical(validate_query(x)$gene_name, "A"))
expect_error(validate_query(x[, -3]), "missing required")
x$Spearman_P_Value[1] <- .1
expect_error(validate_query(x), "BOTH")
x$Spearman_P_Value[1] <- .01
x$Pearson_Correlation[1] <- 2
expect_error(validate_query(x), "Out-of-range")
o <- "GO_Biological_Process_2023"
d <- data.frame(Term="Histone Methylation", Overlap="1/2", P.value=.01, Adjusted.P.value=.4, Odds.Ratio=2, Combined.Score=3, Genes="A", check.names=FALSE)
stopifnot(normalize_go(d, o, "A")$Adjusted.P.value == .4)
names(d) <- paste0(o, ".", names(d))
stopifnot(normalize_go(d, o, "A")$Term == "Histone Methylation")
expect_error(normalize_go(d, o, "B"), "provenance mismatch")
case_result <- d; case_result[[paste0(o, ".Genes")]] <- "C1ORF226"
stopifnot(nrow(normalize_go(case_result, o, "C1orf226")) == 1)
expect_error(normalize_go(case_result, o, "CXorf49B"), "provenance mismatch")
expect_error(normalize_go(d[, -4], o, "A"), "Malformed")
base <- tempfile("go_test_")
dir.create(base)
cache <- file.path(base, "cache")
root <- file.path(base, "output")
dir.create(cache); dir.create(root)
q <- data.frame(gene_name="A",Pearson_Correlation=-.6,Spearman_Correlation=-.6,Pearson_P_Value=.01,Spearman_P_Value=.01)
write.csv(q,file.path(root,"Niclosamide_ESCC_allgene_significant_results.csv"),row.names=FALSE)
write.csv(q,file.path(cache,"negative_pearson_genes.csv"),row.names=FALSE)
for (o in c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")) {
  d <- data.frame(Term="Histone Methylation", Overlap="1/2", P.value=.01, Adjusted.P.value=.4, Odds.Ratio=2, Combined.Score=3, Genes="A", check.names=FALSE)
  write.csv(d,file.path(cache,paste0("enrichR_result_",o,"_Negative.csv")),row.names=FALSE)
}
Sys.setenv(DEPMAP_OUTPUT_DIR=root,DEPMAP_GO_CACHE_DIR=cache)
run_go()
stopifnot(nrow(read.csv(file.path(root,"GO","Methylation_GO_results.csv"))) == 3,
          nrow(read.csv(file.path(root,"GO","Methylation_GO_FDR_significant.csv"))) == 0,
          nrow(read.csv(file.path(root,"GO","EHMT2_GO_terms.csv"))) == 0)
q$gene_name <- "B"
write.csv(q,file.path(cache,"negative_pearson_genes.csv"),row.names=FALSE)
expect_error(run_go(),"not empty")
unlink(file.path(root,"GO"),recursive=TRUE)
expect_error(run_go(),"does not match")
q$gene_name <- "A"
write.csv(q,file.path(cache,"negative_pearson_genes.csv"),row.names=FALSE)
writeLines(c("mode=live_official_enrichr_api", "query_readback=verified", "query_gene_count=1", "query_md5=wrong"),file.path(cache,"GO_query_metadata.txt"))
expect_error(run_go(), "metadata fingerprint mismatch")
stopifnot(!dir.exists(file.path(root,"GO")))
unlink(base,recursive=TRUE)
cat("GO validation and offline integration tests passed.\n")
