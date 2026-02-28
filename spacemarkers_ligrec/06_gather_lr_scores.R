library(grid)
library(ComplexHeatmap)

filter_method = "top_n"
top_n = 30

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 1) {
    stop("Usage: Rscript 06_gather_lr_scores.R <spacemarkers_outputs_dir>\n")
} else {
    lr_dir <- args[1]
}


lrs <- list.files(
        lr_dir,
        recursive = TRUE,
        full.names = TRUE,
        pattern = "LRscores.rds"
    )

sample_names <- unlist(regmatches(m=regexec(pattern='HC[0-9]+', text=lrs), x=lrs))
lr_dfs <- lapply(lrs, FUN=function(x) as.data.frame(readRDS(x)))

#get single direction lr scores from df list
get_lr <- function(lrdfs, samples, string){
    res_list <- list()
    for (l in seq_along(lrdfs)) {
        lr <- lrdfs[[l]][,string, drop=FALSE]
        colnames(lr) <- samples[l]
        lr[["index"]] <- rownames(lr)
        rownames(lr) <- NULL
        res_list[[l]] <- lr
    }
    res_df <- Reduce(function(x, y) merge(x, y, by="index", all=TRUE), res_list)
    rownames(res_df) <- res_df[,"index"]
    res_df[,"index"] <- NULL

    res_df
}

lrs_F2P <- get_lr(lr_dfs, sample_names, 'FIBROBLASTS_to_PDAC')
lrs_P2F <- get_lr(lr_dfs, sample_names, 'PDAC_to_FIBROBLASTS')

write.csv(lrs_F2P, file=paste0( "./", "PDAC_near_FIBROBLASTS_significant_spacemarkers.csv"))
write.csv(lrs_P2F, file=paste0("./", "FIBROBLASTS_near_PDAC_significant_spacemarkers.csv"))

n_F2P <- rowSums(!is.na(lrs_F2P))
n_P2F <- rowSums(!is.na(lrs_F2P))

keep_P2F <- names(which(n_P2F==6))
keep_F2P <- names(which(n_F2P==6))

slrs_F2P <- lrs_F2P[keep_F2P,]
slrs_P2F <- lrs_P2F[keep_P2F,]

slrs_F2P[is.na(slrs_F2P)] <- 0
slrs_P2F[is.na(slrs_P2F)] <- 0

colnames(slrs_F2P) <- sample_names
colnames(slrs_P2F) <- sample_names

responders <- c("HC05","HC07","HC08")
nonresponders <- c("HC01","HC03","HC04")


rank_by_effect_size <- function(slrs, group1, group2) {
    # Compute effect sizes for each row
    eff_size <- apply(slrs, 1, function(x) effsize::cohen.d(x[group1], x[group2], pooled=TRUE))

    # Helper to convert effsize object to data frame
    effsize_to_df <- function(effsize) {
    data.frame(cohen.d = effsize$estimate,
                ci_lower = effsize$conf.int[1],
                ci_upper = effsize$conf.int[2],
                var = effsize$var,
                magnitude = effsize$magnitude)
    }

    # Apply helper
    df_eff_size <- do.call(rbind, lapply(eff_size, effsize_to_df))

    res <- merge(slrs, df_eff_size, by="row.names", all=TRUE)
    rownames(res) <- res$Row.names
    res$Row.names <- NULL

    res[order(abs(res$cohen.d), decreasing=TRUE),]
}

cod_F2P <- rank_by_effect_size(slrs_F2P, nonresponders, responders)
write.csv(cod_F2P, file="FIBROBLASTS_near_PDAC_ranked_LR_scores.csv")

plot_F2P <- cod_F2P[1:top_n,c(nonresponders, responders)]

png("LR_scores_F2P_heatmap.png", width=600, height=900, res=150)
plt <-  ComplexHeatmap::Heatmap(
                plot_F2P,
                name = paste("LR Scores Fibroblasts to PDAC"),
                col = circlize::colorRamp2(c(min(plot_F2P,na.rm=TRUE), 0, max(plot_F2P,na.rm=TRUE)), c("blue", "white", "red")),
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "LR Score"),
                column_title = "Sample",
                row_title = "Gene",
                cluster_columns = FALSE,
            )
print(plt)
dev.off()


cod_P2F <- rank_by_effect_size(slrs_P2F, nonresponders, responders)
write.csv(cod_P2F, file="PDAC_near_FIBROBLASTS_ranked_LR_scores.csv")
plot_P2F <- cod_P2F[1:top_n,c(nonresponders, responders)]

png("LR_scores_P2F_heatmap.png", width=600, height=900, res=150)
plt <-  ComplexHeatmap::Heatmap(
                plot_P2F,
                name = paste("LR Scores PDAC to Fibroblasts"),
                col = circlize::colorRamp2(c(min(plot_P2F,na.rm=TRUE), 0, max(plot_P2F,na.rm=TRUE)), c("blue", "white", "red")),
                show_row_names = TRUE,
                show_column_names = TRUE,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 10),
                heatmap_legend_param = list(title = "LR Score"),
                column_title = "Sample",
                row_title = "Gene",
                cluster_columns = FALSE,
            )
print(plt)
dev.off()
