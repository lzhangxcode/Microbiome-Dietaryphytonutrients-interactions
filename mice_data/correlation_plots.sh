#1. DNA heatmap plot 
#1.1 significantly differentially abundant EC between DSS+Vehicle vs DSS+Strawberry at D18 
Rscript sig_EC_specific_day_comparison.R metagenomics_join_enzyme_clr.tsv D18 G2 G3 mice_metatable.txt EC_strawberry_both_rename_column.tsv \
./sigec_day18/

#1.2 significantly differentially abundant EC between DSS+Vehicle vs DSS+Strawberry at D6 
Rscript sig_EC_specific_day_comparison.R metagenomics_join_enzyme_clr.tsv D6 G2 G3 mice_metatable.txt EC_strawberry_both_rename_column.tsv \
./sigec_day6/

#1.3 correlation with DAI score 
Rscript correlation_target_ec_clr.R metagenomics_join_enzyme_clr.tsv \
EC_strawberry_both.tsv mice_metatable_G2G3.txt summary_SPF_result.xlsx ./G2G3_correlation/

#1.4 correlation with HE score 
Rscript correlation_target_ec_HE_clr.R metagenomics_join_enzyme_clr.tsv \
EC_strawberry_both.tsv.tsv mice_metatable_G2G3.txt HE_score.xlsx Total,score ./HE_correlation/ G2,G3

#1.5 heatmap 
strawberry_summary_sig.EC_day6andday18_DNA.R
strawberry_correlation_DNA.R

#2. RNA heatmap plot 
#2.1 significantly differentially abundant EC between DSS+Vehicle vs DSS+Strawberry at D18 
Rscript sig_EC_specific_day_comparison.R metatranscriptomics_join_enzyme_clr.tsv D18 G2 G3 \
mice_metatable.txt EC_strawberry_both_metat_rename_column.tsv ./sigec_day18/

#2.2 significantly differentially abundant EC between DSS+Vehicle vs DSS+Strawberry at D6 
Rscript sig_EC_specific_day_comparison.R metatranscriptomics_join_enzyme_clr.tsv D6 G2 G3 \
mice_metatable.txt EC_strawberry_both_metat_rename_column.tsv ./sigec_day6/

#2.3 correlation with DAI score 
Rscript correlation_target_ec_clr.R metatranscriptomics_join_enzyme_clr.tsv EC_strawberry_both_metat.tsv \
mice_metatable_G2G3.txt summary_SPF_result.xlsx ./G2G3_correlation/

#2.4 correlation with HE score 
Rscript correlation_target_ec_HE_clr.R metatranscriptomics_join_enzyme_clr.tsv EC_strawberry_both_metat.tsv \
mice_metatable_G2G3.txt HE_score.xlsx Total,score ./HE_correlation/ G2,G3
 
#2.5 draw heatmap 

strawberry_correlation_RNA.R



