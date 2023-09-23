################################################################################
# This scripts convert the complete results table to two other formats:
#      1- Samples Vs. ASVs table
#      2- Samples Vs. IDs table
# It can be used after performing manual curation on the complete results table.
################################################################################
# Support: Heron OH - heronoh@gmail.com
################################################################################
#read original complete table ----
wdir <- "/home/..."

smp_abd_ID <- readxl::read_xlsx(path = paste0(wdir,"ecomol-iSeq16_03112022_CIAT-complete_analysis_results-2022-11-21.xlsx"))


################################################################################
#function to either sum or unique by column type ----
sum_uniq <- function(vec=vec){

  if (is.character(vec)==TRUE) {
    suniq <- BiocGenerics::unique(vec)
  }
  if (is.numeric(vec)==TRUE) {
    suniq <- sum(vec)
  }
  return(suniq)
}
################################################################################
# 1- generate Samples Vs. ASVs table from complete table ----

smp_abd_ID_eco <-
  smp_abd_ID %>%
  mutate(Identification = if_else(Identification %in% c(NA,"NA"),"NA",Identification)) %>%
  select(-c("Relative abundance to all samples",
            "Sample total abundance",
            "ASV absolute abundance","metadata_1","metadata_2","metadata_3","metadata_4","metadata_5","obs",
            "Curated ID",
            "Final ID",
            "PCR control",
            "Prop. to PCR control",
            "Prop. to Ext control",
            "Prop. to Filt control",
            "Possible contamination",
            "Primer expected length",
            "Type")) %>%
  pivot_wider(
    id_cols = c("Identification","Primer","Primer name","Read origin",
                "Possible bacteria",
                "Superkingdom (BLASTn)",
                "Kingdom (BLASTn)",
                "Phylum (BLASTn)",
                "Subphylum (BLASTn)",
                "Class (BLASTn)",
                "Subclass (BLASTn)",
                "Order (BLASTn)",
                "Suborder (BLASTn)",
                "Family (BLASTn)",
                "Subfamily (BLASTn)",
                "Genus (BLASTn)",
                "max_tax","BLAST ID","Blast pseudo-score",
                "1_subject header","1_subject","1_indentity","1_qcovhsp",
                "1_length","1_mismatches","1_gaps","1_query start","1_query end",
                "1_subject start","1_subject end","1_e-value","1_bitscore","2_subject header",
                "2_subject","2_indentity","2_qcovhsp","2_length","2_mismatches",
                "2_gaps","2_query start","2_query end","2_subject start","2_subject end",
                "2_e-value","2_bitscore","3_subject header","3_subject","3_indentity",
                "3_qcovhsp","3_length","3_mismatches","3_gaps","3_query start",
                "3_query end","3_subject start","3_subject end","3_e-value","3_bitscore",
                "Remove","ASV Size (pb)","ASV header","ASV (Sequence)","OTU"),
    values_from ="Relative abundance on sample",
    values_fn = sum_uniq,
    names_from = File_name,
    names_sort = TRUE,
    names_prefix = "SAMPLE ") %>%
  relocate(c("Primer","Primer name",
             "Read origin",
             "Superkingdom (BLASTn)",
             "Kingdom (BLASTn)",
             "Phylum (BLASTn)",
             "Subphylum (BLASTn)",
             "Class (BLASTn)",
             "Subclass (BLASTn)",
             "Order (BLASTn)",
             "Suborder (BLASTn)",
             "Family (BLASTn)",
             "Subfamily (BLASTn)",
             "Genus (BLASTn)",
             "max_tax","BLAST ID",
             "Possible bacteria", "Identification","Blast pseudo-score",
             starts_with("SAMPLE "))) %>%
  mutate_if(is.numeric , replace_na, replace = 0)


writexl::write_xlsx(x = smp_abd_ID_eco,
                    path = paste0(wdir,"ecomol-iSeq16_03112022_CIAT-Samples_x_ASVs",Sys.Date(),".xlsx"),
                    col_names = TRUE,format_headers = TRUE)


################################################################################
# 2 - generate Samples Vs. ASVs table from complete table ----
################################################################################

smp_abd_ID_eco_ID <- smp_abd_ID %>%
  mutate(Identification = if_else(Identification %in% c(NA,"NA"),"NA",Identification)) %>%
  select(-c("Relative abundance to all samples",
            "Sample total abundance",
            "ASV absolute abundance",
            "metadata_1","metadata_2","metadata_3","metadata_4","metadata_5","obs",
            "Curated ID",
            "Final ID",
            "PCR control",
            "Prop. to PCR control",
            "Prop. to Ext control",
            "Prop. to Filt control",
            "Possible contamination",
            "Primer expected length",
            "Type",
            "blast ID Origin",
            # "Read origin",
            "1_subject header","1_subject","1_indentity","1_qcovhsp",
            "1_length","1_mismatches","1_gaps","1_query start","1_query end",
            "1_subject start","1_subject end","1_e-value","1_bitscore","2_subject header",
            "2_subject","2_indentity","2_qcovhsp","2_length","2_mismatches",
            "2_gaps","2_query start","2_query end","2_subject start","2_subject end",
            "2_e-value","2_bitscore","3_subject header","3_subject","3_indentity",
            "3_qcovhsp","3_length","3_mismatches","3_gaps","3_query start",
            "3_query end","3_subject start","3_subject end","3_e-value","3_bitscore",
            "Remove","ASV Size (pb)","ASV header","ASV (Sequence)","OTU"
  )) %>%
  pivot_wider(
    id_cols = c("Primer","Primer name","Read origin",
                "Identification",
                "Possible bacteria",
                "Kingdom (BLASTn)",
                "Phylum (BLASTn)",
                "Subphylum (BLASTn)",
                "Class (BLASTn)",
                "Subclass (BLASTn)",
                "Order (BLASTn)",
                "Suborder (BLASTn)",
                "Family (BLASTn)",
                "Subfamily (BLASTn)",
                "Genus (BLASTn)",
                "BLAST ID","max_tax"),
    values_from = c("Relative abundance on sample"),
    values_fn = sum_uniq,
    names_from = "File_name",
    names_prefix = "SAMPLE ",
    names_sort = TRUE) %>%
  relocate(c("Primer","Primer name",
             "Read origin",
             "Kingdom (BLASTn)",
             "Phylum (BLASTn)",
             "Subphylum (BLASTn)",
             "Class (BLASTn)",
             "Subclass (BLASTn)",
             "Order (BLASTn)",
             "Suborder (BLASTn)",
             "Family (BLASTn)",
             "Subfamily (BLASTn)",
             "Genus (BLASTn)",
             "max_tax","BLAST ID",
             "Possible bacteria", "Identification",
             # "Blast pseudo-score",
             starts_with("SAMPLE ")
  )) %>%
  mutate_if(is.numeric , replace_na, replace = 0)



writexl::write_xlsx(x = smp_abd_ID_eco_ID,
                    path = paste0(wdir,"ecomol-iSeq16_03112022_CIAT-Samples_x_IDs",Sys.Date(),".xlsx"),
                    col_names = TRUE,format_headers = TRUE)
