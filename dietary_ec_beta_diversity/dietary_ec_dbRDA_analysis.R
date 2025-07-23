# the aim of this script is to do dbRDA analysis for dietary ec profile 
# author: lu.zhang
# date: 02,March,2024 
 

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
	          stop("Rscript script.R modified_ec_profile output_dir", call.=FALSE)
}


############################################################

# loading package                                          #

############################################################

library(tidyverse)
library(vegan)
library(ggrepel)


############################################################

# prepare input 

############################################################


# dietary ec profile 
modified_ec_profile_beta_input <- readRDS(args[1])

# output_dir 
output_dir <- args[4]
output_dir <- "./dbRDA/"



############################################################

# dbRDA - univariate

############################################################
univariate_dbRDA <- function(rda_input = NA, metadata = NA){

    metadata$depth <- metadata$nonhost
    metadat_sel_specific <- metadata %>% 
        select(age, BMI, continent, gender, depth, country) %>% 
        as.data.frame()


    all <- c()
    for(i in 1:ncol(metadat_sel_specific)){
        print(i)
        cap <- capscale(rda_input ~ metadat_sel_specific[,i], distance="bray", na.action=na.omit) 
        av <- cap %>% 
            anova.cca %>% 
            as.data.frame
        pval <- av$`Pr(>F)`[1]
        Fstat <- av$`F`[1]
        r2 <- RsquareAdj(cap)[[1]]
        adjr2 <- RsquareAdj(cap)[[2]]
        all <- rbind(all, cbind(Fstat,r2,adjr2,pval))
    }

    rownames(all) <- colnames(metadat_sel_specific)
    all <- as.data.frame(cbind(all, padj=p.adjust(all[,"pval"], method="fdr")), stringsAsFactors=F)
    significant <- all[all[,"padj"]<0.05,]
    significant <- na.omit(significant[order(significant$adjr2, decreasing = T),])
    print(paste0(nrow(significant), " variables can individually explain part of the (active) microbiome variation in these samples."))

    uni_result <- list()
    uni_result[["all"]] <- all
    uni_result[["significant"]] <- significant 
    return(uni_result)
}

modified_univariate <- univariate_dbRDA(rda_input = modified_ec_profile_beta_input[[2]], metadata = modified_ec_profile_beta_input[[1]])

############################################################

# dbRDA - multivariate

############################################################

multivariable_model_dbRDA <- function(dbrda_input = NA, metadata = NA){
    metadata$depth <- metadata$nonhost
    metadat_sel_specific <- metadata %>% 
        select(age, BMI, continent, gender, depth, country) %>% # , subject_id, # sequencing_platform, 
        as.data.frame()

    db_rda <- capscale(dbrda_input ~ country + continent + age + gender + BMI + depth, metadat_sel_specific, distance = 'bray') # remove add = TRUE 
    set.seed(1234)
    db_rda_forward_pr <- ordiR2step(capscale(dbrda_input~ 1, metadat_sel_specific, distance = 'bray'), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)
    
    # Print call and anova
    db_rda_forward_pr$call
    #example : capscale(formula = rda_input ~ continent + sequencing_platform + depth + gender + age + BMI, data = metadat_sel_specific, distance = "bray")
    db_rda_forward_pr$anova

    #Rsquare of the model 
    db_rda_forward_pr_Radj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared
    

  
    #Permutation Tests of RDA Results
    db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)

    dbrda_result <- list()
    dbrda_result[["forward"]] <- db_rda_forward_pr
    dbrda_result[["rsquare"]] <- db_rda_forward_pr_Radj
    dbrda_result[["dbrda_p"]] <- db_rda_forward_pr_test

    return(dbrda_result)

}

modified_multivariate <- multivariable_model_dbRDA(dbrda_input = modified_ec_profile_beta_input[[2]], metadata = modified_ec_profile_beta_input[[1]])


############################################################

# dbRDA - ordinate plot 

############################################################

ordination_axis_explanation <- function(db_rda_forward_prv = NA, metadata = NA, model_p = NA, tag = NA){

    metadata$depth <- metadata$nonhost
    metadat_sel_specific <- metadata %>% 
        select(age, BMI, continent, gender, depth, country) %>% # , subject_id, # sequencing_platform, 
        as.data.frame()
    
    #==============PART 1 axis explanation power ===========================================================================

    db_rda_forward_pr.scaling1 <- summary(db_rda_forward_prv, scaling = 1)
    db_rda_forward_pr.site <- data.frame(db_rda_forward_pr.scaling1$sites)[1:2] # extract samples x,y axis. 
    db_rda_forward_pr.env <- data.frame(db_rda_forward_pr.scaling1$biplot)[1:2] # environment catogories : in my case, age,gender,bmi, continent

    db_rda_forward_pr.env$name <- rownames(db_rda_forward_pr.env)
    #db_rda_forward_pr.spe$name <- rownames(db_rda_forward_pr.spe)
    db_rda_forward_pr.site$name <- rownames(db_rda_forward_pr.site)

    db_rda_forward_pr.site$gender <- metadat_sel_specific$gender
    db_rda_forward_pr.site$country <- metadat_sel_specific$country
    db_rda_forward_pr.site$continent <- metadat_sel_specific$continent

    db_rda_forward_pr_Radj <- RsquareAdj(db_rda_forward_prv)$adj.r.squared
    exp_adj <- RsquareAdj(db_rda_forward_prv)$adj.r.squared * db_rda_forward_prv$CCA$eig/sum(db_rda_forward_prv$CCA$eig)
    rda1_exp <- paste('Axis.1:', round(exp_adj[1]*100, 2), '%')
    rda2_exp <- paste('Axis.2:', round(exp_adj[2]*100, 2), '%')
    exp_raw <- db_rda_forward_prv$CCA$eig/sum(db_rda_forward_prv$CCA$eig)
    rda1_raw <- paste('Axis.1:', round(exp_raw[1]*100, 2), '%')
    rda2_raw <- paste('Axis.2:', round(exp_raw[2]*100, 2), '%')


    rda_axis <- c(rda1_exp, rda2_exp, rda1_raw, rda2_raw)
    #==============PART 2 extract the centroid of the factor variables ====================================================

    #extract the centroid of the factor variables 

    # modify the environment variables, some are factors, shouldn't be draw arrow.
    db_rda_forward_pr.env.all <- db_rda_forward_pr.env
    db_rda_forward_pr.env <- db_rda_forward_pr.env.all[c("BMI","depth"),]
    #rownames(db_rda_forward_pr.env)[1]<-"Age"

    # extract the continent & the gender axis 
    db_rda_forward_pr.env.factor <- db_rda_forward_pr.env.all[!db_rda_forward_pr.env.all$name %in%c("depth","BMI"),]
    db_rda_forward_pr.env.factor.modified <- as.data.frame(db_rda_forward_pr.scaling1$centroids)
    rownames(db_rda_forward_pr.env.factor.modified) <- rownames(db_rda_forward_pr.env.factor.modified) %>% 
        str_remove("country") %>% 
        str_remove("gender") 
    
    #==============PART 3 the axis and title setting =======================================================================


    db_rda_forward_pr_test <- model_p

    line_3 = bquote("dbRDA anova:"~"F"~"= "~.(round(db_rda_forward_pr_test$F[1],3))~"P = "~.(db_rda_forward_pr_test$`Pr(>F)`[1]))
    line_4 = bquote("Constrained explained variation:"~ R^2 ~"= "~.(paste(round(db_rda_forward_pr_Radj*100, 2), '%')))
    nonmetric_label = list(line_3,line_4)

    coord_x <- Inf
    coord_x <- c(1.0 * min(db_rda_forward_pr.site$CAP1), 0.9*min(db_rda_forward_pr.site$CAP1))
    coord_y <- c(1.0 * min(db_rda_forward_pr.site$CAP2),
               0.9 *min(db_rda_forward_pr.site$CAP2),
               0.8*max(db_rda_forward_pr.site$CAP2),
               0.7*max(db_rda_forward_pr.site$CAP2))
    title = sprintf("dbRDA (%s) (%s)",tag,"bray curtis distance")

    cbPalette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF")
    envPalette <- "black"

    #==============PART 4 draw ordination plot  ===========================================================================


    options(ggrepel.max.overlaps = 20)

    p2 <- ggplot(db_rda_forward_pr.site, aes(CAP1, CAP2)) +
        geom_point(aes(color = continent, shape = gender),position = "jitter") + 
        scale_color_manual(values = cbPalette[1:4]) +
        labs(x = rda1_exp, y = rda2_exp) +
        geom_hline(yintercept = 0,linetype=1,color="grey")+geom_vline(xintercept = 0,linetype=1,color="grey")+
        scale_shape_manual(values = c(1,16))+
        geom_segment(data = db_rda_forward_pr.env, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.7, color = 'black') +
        geom_text_repel(data=db_rda_forward_pr.env,aes(x = CAP1*1.2, 
                                                    y = CAP2*1, 
                                                    label = rownames(db_rda_forward_pr.env)),size=4, 
                                                    color = "black", fontface = 'bold') + 
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[1], y = db_rda_forward_pr.env.factor.modified$CAP2[1]), color = "black", shape = 15, size = 2) + 
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[2], y = db_rda_forward_pr.env.factor.modified$CAP2[2]), color = "black", shape = 15, size = 2) + 
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[3], y = db_rda_forward_pr.env.factor.modified$CAP2[3]), color = "black", shape = 15, size = 2) + 
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[4], y = db_rda_forward_pr.env.factor.modified$CAP2[4]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[5], y = db_rda_forward_pr.env.factor.modified$CAP2[5]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[6], y = db_rda_forward_pr.env.factor.modified$CAP2[6]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[7], y = db_rda_forward_pr.env.factor.modified$CAP2[7]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[8], y = db_rda_forward_pr.env.factor.modified$CAP2[8]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[9], y = db_rda_forward_pr.env.factor.modified$CAP2[9]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[10], y = db_rda_forward_pr.env.factor.modified$CAP2[10]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[11], y = db_rda_forward_pr.env.factor.modified$CAP2[11]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[12], y = db_rda_forward_pr.env.factor.modified$CAP2[12]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[13], y = db_rda_forward_pr.env.factor.modified$CAP2[13]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[14], y = db_rda_forward_pr.env.factor.modified$CAP2[14]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[15], y = db_rda_forward_pr.env.factor.modified$CAP2[15]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[16], y = db_rda_forward_pr.env.factor.modified$CAP2[16]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[17], y = db_rda_forward_pr.env.factor.modified$CAP2[17]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[18], y = db_rda_forward_pr.env.factor.modified$CAP2[18]), color = "black", shape = 15, size = 2) +
        geom_point(aes(x = db_rda_forward_pr.env.factor.modified$CAP1[19], y = db_rda_forward_pr.env.factor.modified$CAP2[19]), color = "black", shape = 15, size = 2) +
	geom_text_repel(data=db_rda_forward_pr.env.factor.modified,aes(x = CAP1*1, 
                                                        y = CAP2, 
                                                        label = rownames(db_rda_forward_pr.env.factor.modified)),size=4, 
                                                        color = "black", fontface = 'bold') +
        theme_bw() +
        annotate("text", x = -0.40, y = 0.55, label = nonmetric_label[[1]]) + 
        scale_fill_manual(values = envPalette)+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
        theme(
            legend.box.background = element_rect(fill='#fffafa'),
            legend.position = c(.02, .40),
            legend.justification = c("left", "top"),
            legend.margin = margin(6, 6, 6, 6), legend.key.size = unit(4,"mm"),
            legend.background = element_blank()
        )


        p2_raw_axis <- p2 + labs(x = rda1_raw, y = rda2_raw)
        p2.modified <- ggExtra::ggMarginal(p2, type="boxplot", groupColour = TRUE, groupFill = TRUE, size = 8) 
        p2.modified.raw <- ggExtra::ggMarginal(p2_raw_axis, type="boxplot", groupColour = TRUE, groupFill = TRUE, size = 8) 

        result <- list()
        result[["site"]] <- db_rda_forward_pr.site
        result[["factor"]] <- db_rda_forward_pr.env.factor.modified
        result[["env"]] <- db_rda_forward_pr.env
        result[["p"]] <- p2
        result[["p_marginal"]] <- p2.modified
        result[["p_raw"]] <- p2.modified.raw 

        pdf(paste(output_dir, paste0(tag, "_dbrda_country_gender.pdf"), sep = "/"))
        print(p2.modified)
        dev.off()

	return(result)

}

ordination_plot_modified <- ordination_axis_explanation(db_rda_forward_prv = modified_multivariate[["forward"]],
                                                        metadata = modified_ec_profile_beta_input[[1]],
                                                        model_p = modified_multivariate[["dbrda_p"]],
                                                        tag = "modified_ec")

############################################################

# dbRDA - bar plot 

############################################################

draw_barplot_for_dbRDA_result <- function(all_no_platform = NA, db_rda_forward_prv = NA, significant_no_platform = NA, filename = NA){


    # combine the univariate and multivariable results together
    significant_no_platform <- all_no_platform[all_no_platform[,"padj"]<0.05,]
    significant_no_platform <- na.omit(significant_no_platform[order(significant_no_platform$adjr2, decreasing = T),])
    print(paste0(nrow(significant_no_platform), " variables can individually explain part of the (active) microbiome variation in these samples."))


    significant_no_platform <- as_tibble(cbind(significant_no_platform, covariates=rownames(significant_no_platform)))
    anova_table_no_platform <- as_tibble(cbind(db_rda_forward_prv$anova, covariates=gsub(rownames(db_rda_forward_prv$anova),
                                                                pattern="\\+ ",
                                                                replacement="")))

    write_tsv(anova_table_no_platform, paste(output_dir, paste(filename, "multivaraible.tsv", sep = "_"), sep = "/"))
    write_tsv(significant_no_platform, paste(output_dir, paste(filename, "univaraible.tsv", sep = "_"), sep = "/"))

    significant_no_platform <- significant_no_platform %>%
        full_join(anova_table_no_platform, by="covariates")
    significant_no_platform <- significant_no_platform %>%
        filter(covariates!="<All variables>")

    significant_no_platform[is.na(significant_no_platform$R2.adj),"R2.adj"] <- max(significant_no_platform$R2.adj,na.rm=T)
    significantL_no_platform <- gather(significant_no_platform, key = "Class", value="EffectSize", adjr2, R2.adj)
    significantL_no_platform[significantL_no_platform$Class=="adjr2", "Class"] <- "Individual"
    significantL_no_platform[significantL_no_platform$Class=="R2.adj", "Class"] <- "Cumulative"
    significantL_no_platform_sel <- significantL_no_platform %>% filter(Class == "Individual")
    category <- ifelse(significantL_no_platform_sel$covariates == "depth", "Laboratory", "Demographic")
    significantL_no_platform$category <- rep(category,2)

    positions <- rev.default(significant_no_platform$covariates[order(significant_no_platform$R2.adj, significant_no_platform$AIC,
                                                      -significant_no_platform$adjr2,decreasing=F)])
    # positions
    # [1] "continent" "age"       "BMI"       "gender"    "depth"     "country"

    color_list3 <- c("Demographic"="navyblue", "Laboratory"="gold")

    significantL_no_platform$Axis_color <- significantL_no_platform$category

    axis_color2 <- dplyr::distinct(significantL_no_platform, covariates, Axis_color) %>%
        tibble::deframe()

    significantL_2 <- significantL_no_platform
    significantL_2$covariates <- factor(significantL_2$covariates, levels = positions)
    saveRDS(significantL_2, paste(output_dir, paste(filename, "without_sequencingplatform_covariates_barplot_input.rds", sep = "_"), sep = "/"))
    barplot_univariates_multi <-  ggplot(significantL_2, aes(x = covariates, y = EffectSize)) +
            geom_bar(aes(color = NULL, alpha=Class), # fill = category,
            stat = "identity", position = position_dodge()) +
            ggpubr::rotate() +
            theme_bw() +
            #scale_x_discrete(limits = positions) +
            scale_alpha_manual(values=c(0.6,1)) +
                geom_vline(xintercept = nrow(significant_no_platform)-(nrow(anova_table_no_platform)-1.5))+ # anova here contains <all variables>, so needs to + 1 and the line should be in the middle, so should + 0.5
            labs(x="Covariates", y="Effect size")  +
                scale_x_discrete(labels = ~ glue::glue("<span style='color: {color_list3[axis_color2[.x]]}'>{.x}</span>")) +
            ggnewscale::new_scale_fill() +
            geom_col(aes(y = 0, fill = Axis_color)) +
            scale_fill_manual(values = color_list3[unique(significantL_2$Axis_color)], guide = guide_legend(order= 2)) +
            theme(axis.text.y = ggtext::element_markdown())

    ggsave(barplot_univariates_multi, filename = paste(output_dir, paste(filename, "without_sequencingplatform_covariates.pdf", sep = "_"), sep = "/"), device="pdf", width=7, height=6)
    return(barplot_univariates_multi)

}


modified_barplot <- draw_barplot_for_dbRDA_result(all_no_platform = modified_univariate[["all"]],
                                                    db_rda_forward_prv = modified_multivariate[["forward"]],
                                                    significant_no_platform = modified_univariate[["significant"]],
                                                    filename = "modified")



############################################################

# dbRDA - percentage

############################################################

percentage_plot <- function(db_rda_forward_pr_Radj = NA, tag = NA){
        percentage <- tibble(Proportion = c(db_rda_forward_pr_Radj*100, (1- db_rda_forward_pr_Radj)*100), Category = c("Explained", "unknown"), 
                     analysis = "dbRDA")
        per_plot <- ggplot(percentage, aes(x = analysis, y = Proportion, fill = Category)) +
            geom_col(width = 0.2) +
            scale_fill_manual(values=c("red","darkgrey")) +
            geom_text(aes(label = paste0(round(Proportion,3), "%")),
                position = position_stack(vjust = 0.5)) +
            theme(legend.position = "right", 
                legend.title = element_blank()) +
            theme(axis.title.y = element_text(margin = margin(r = 1))) +
            ylab("Microbiome Variation(%)") +
            xlab(NULL) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
            annotate("text", x = c(0.5,0.5), y = c(50, 98), label = c("unknown","explained"), angle = 90) +
            theme(legend.position = "none") +
            theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank())

        pdf(paste(output_dir, paste0(tag, "_percentage.pdf"), sep = "/"))
        print(per_plot)
        dev.off()

        return(per_plot)

}


modifiedec_percentage_plot <- percentage_plot(db_rda_forward_pr_Radj = unmodified_multivariate[["rsquare"]], tag = "unmodified")


############################################################

# combine plot together 

############################################################


# combine plots 


pdf(paste(output_dir, "modified_ec_dbRDA_combine.pdf", sep = "/"), width = 14, height = 8)

ggpubr::ggarrange(modifiedec_percentage_plot, modified_barplot, ordination_plot_modified[["p_marginal"]], widths = c(0.7,3,4), nrow = 1)

dev.off()



############################################################

# saving RData 

############################################################

save.image(paste0(output_dir, "step16_final_dbRDA_manual_check.RData"))


