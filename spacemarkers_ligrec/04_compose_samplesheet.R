hd_dirs <- list.dirs("~/data/dpt", recursive = TRUE)
hd_dirs <- hd_dirs[grepl("square_016um$", hd_dirs)]

hd_samples <- regmatches(hd_dirs,regexpr("HC[0-9][0-9]",hd_dirs))


rctd_files <- list.files("~/data/dpt", full.names = TRUE, pattern = "^rctd_cell_types.csv", recursive = TRUE)
rctd_samples <- regmatches(rctd_files, regexpr("HC[0-9][0-9]", rctd_files))

#samplesheet 
samplesheet <- merge(data.frame(sample = hd_samples, data_dir = hd_dirs),
      data.frame(sample = rctd_samples, annotation_file = rctd_files),
      by = "sample", all = TRUE)

write.csv(samplesheet, "samplesheet.csv", row.names = FALSE, quote = FALSE)
