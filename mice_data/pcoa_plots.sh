# DNA D6 : DSS + Vehicle vs DSS + strawberry
Rscript PCoA_plot_for_target_groups_target_ec.R mice_metatable.txt metagenomics_join_enzyme_cpm.tsv EC_strawberry_both.tsv \
group G2,G3 metaG ./pcoa_group/
    
# DNA DSS + strawberry : D6 vs D18 
Rscript PCoA_plot_for_target_groups_target_ec.R mice_metatable_strawberry.txt metagenomics_join_enzyme_cpm.tsv EC_strawberry.tsv \
day D6,D18 MetaG_straw ./pcoa_day/ 
 
# RNA D6: DSS + Vehicle vs DSS + strawberry
Rscript PCoA_plot_for_target_groups_target_ec.R mice_metatable.txt metatranscriptomics_join_enzyme_cpm.tsv EC_strawberry_both_metat.tsv \
group G2,G3 metaT ./pcoa_group/

# RNA DSS + strawberry : D6 vs D18 
Rscript PCoA_plot_for_target_groups_target_ec.R mice_metatable_strawberry.txt metatranscriptomics_join_enzyme_cpm.tsv EC_strawberry_metat.tsv \
day D6,D18 MetaT_straw ./pcoa_day/ 

