temp_asv_names <- rownames(final_filtered_asv_table_t)
to_keep <- temp_asv_names %in% rownames(meta)
final_filtered_asv_table_t <- final_filtered_asv_table_t[to_keep,]
final_filtered_asv_table_t_orig <- final_filtered_asv_table_t
dim(final_filtered_asv_table_t)
final_filtered_asv_table_t <- final_filtered_asv_table_t[,colSums(final_filtered_asv_table_t)> 0]

taxa_silva_full2 <- taxa_silva_full[colnames(final_filtered_asv_table_t),]
dim(taxa_silva_full2)

temp_asv_names <- rownames(final_filtered_asv_table_t)
to_keep <- temp_asv_names %in% rownames(meta)
final_filtered_asv_table_t <- final_filtered_asv_table_t[to_keep,]
final_filtered_asv_table_t <- final_filtered_asv_table_t[,colSums(final_filtered_asv_table_t)> 0]
