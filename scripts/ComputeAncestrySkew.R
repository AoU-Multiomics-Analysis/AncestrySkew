library(tidyverse)
library(data.table)
library(optparse)


option_list <- list(
    optparse::make_option(c("--AnnotationData"), type="character", default=NULL,
                        help="Annotation data for variants containing ancestry specific AC/AN and cohort level AC/AN", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Output file name", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
AnnotationPath <- opt$AnnotationData 
OutputPrefix <- opt$OutputPrefix
OutputName <- paste0(opt$OutputPrefix, ".AncestrySkew.tsv.gz")

####### LOAD DATA ########
AnnotationDf <- fread(AnnotationPath)

# load annotation data and filter to variants that 
# are identifiable as causal for eQTL signals
SkewInput <- AnnotationDf %>% 
    filter(pip >.9 ) %>% 
    select(variant,contains('gvs')) %>% 
    mutate(across(contains('gvs'),~as.numeric(.))) %>% 
    mutate(gvs_ac_skew_input = gvs_all_ac - gvs_max_ac,gvs_an_skew_input = gvs_all_an -  gvs_max_an)  %>% 
    select(variant,gvs_max_ac,gvs_max_an,contains('skew_input')) %>% 
    na.omit()

# compute odds ratio for observing a variant in the population with the highest allele frequency 
# compared to the rest of the cohort
res <- apply(SkewInput[,c("gvs_max_ac","gvs_max_an","gvs_ac_skew_input","gvs_an_skew_input")], 1, function(x) {
  test <- fisher.test(matrix(x, nrow = 2, byrow = TRUE))
  c(odds_ratio = unname(test$estimate),
    p_value = test$p.value)
})

# merge results 
merged <- cbind(SkewInput, t(res))
merged %>% write_tsv(OutputName)
