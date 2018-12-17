library(cegwas)
library(dplyr)
#setwd("/Users/sding/Documents/cegwas/")
list_of_files<-c('/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_area.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_blob.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_curvature.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_food.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_length.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_speed.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_velocity.csv','/Users/sding/Documents/AggScreening/cegwas/csv_s/swDensityCompare_Tierpsy_width.csv')

for (num_list in seq(1,length(list_of_files))) {
mydata=read.csv(list_of_files[num_list])
strains<-as.character(mydata$strain)
trait_values<-rbind(as.numeric(unlist(mydata[2])),as.numeric(unlist(mydata[3])))
for (num_help in seq(4,length(mydata))) {
  trait_values<-rbind(trait_values,as.numeric(unlist(mydata[num_help])))
}
col_titles<-colnames(mydata)

pheno<-data.frame(col_titles[2:length(mydata)],trait_values)
pheno[,2:ncol(pheno)] <- lapply(pheno[,2:ncol(pheno)], as.numeric)
colnames(pheno) <- c("trait", strains)

processed_phenotypes<-process_pheno(pheno)
mapping_df<-gwas_mappings(processed_phenotypes,cores = 4)
processed_mapping_df<-process_mappings(mapping_df,phenotype_df = processed_phenotypes,CI_size = 50,snp_grouping = 200)

file_name<-sprintf('saved_data_new_features_check%s.Rdata',num_list)
save(pheno,processed_phenotypes,mapping_df,processed_mapping_df, file = file_name)
print(num_list)

}

#Theoretically, 'saved_data_5s_check1.Rdata' is the real values, 5s gap, greedy only, so only 12*3 values are actually tested.
load('/Users/sding/Documents/cegwas/results/saved_data_new_features_checksaved_data_new_features_check_swDensityCompare_Tierpsy_paused.Rdata')
realhits=na.omit(processed_mapping_df)
write.csv(realhits, file = "realhits_new_feat_ran3.csv")

unique(processed_mapping_df$trait)
manplot(processed_mapping_df)

