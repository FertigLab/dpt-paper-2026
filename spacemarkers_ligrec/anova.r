library(jsonlite)

datasets <- fromJSON("datasets.json")
sm <- "spacemarkers_ligrec/PDAC_near_FIBROBLASTS_significant_spacemarkers.csv"

df <- read.csv(sm, row.names = 1)
colnames(df) <- paste0(colnames(df), "BTC_visiumHD")

myaov <- function(dfrow, g1=datasets$responders, g2 = datasets$non_responders) {
    df <- data.frame(response = unlist(dfrow[g1]),
                     no_response = unlist(dfrow[g2])) |> stack()

    if (any(is.na(df$values))) {
        res <- NA
    } else {
        fit <- aov(values ~ ind, data = df)
        res <- summary(fit)[[1]][["Pr(>F)"]][1]
    }
    return(res)
}

aov_pvals <- apply(df, 1, myaov)
aov_adj_pvals <- p.adjust(aov_pvals, method = "BH")

res_df <- data.frame(aov_pval = aov_pvals,
                     aov_adj_pval = aov_adj_pvals)

res_df <- res_df[order(res_df$aov_adj_pval), ]
res_df[res_df$aov_adj_pval < 0.05 & !(is.na(res_df$aov_adj_pval)), ]
