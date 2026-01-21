#*****************************************************************************#
#
# Step 2: Descriptive statistics: table and plots
#
#*****************************************************************************#
op <- options(warn=1)

logfile=file.path(outdir,"2_Descriptive_statistics.log"); 
if(file.exists(logfile)) file.remove(logfile)

starting.message=paste0("\n2. Obtaining descriptive statistics \n")
cat(starting.message)
cat(starting.message,file=logfile)

# Descriptive statistics ------------------------------------------------------
if(nrow(df) < 30){
  msg=paste0("There were only ", nrow(df), " participants left in ", group,". Skipping analysis of this group.")
  message(msg); 
  cat(msg,file=logfile,append=T)
  next
}

#* remove variables that are 100% missing in neuro or non-neuro datasets
data_missing <- colSums(is.na(df))
if(any(data_missing == nrow(df))){
  print(names(df)[data_missing == nrow(df)])
  missingCol <- names(df)[data_missing == nrow(df)]
  df[, missingCol] <- NULL
}else{
  missingCol = NULL
}

#* record the numbers of missing values for each column
to_write = data.frame(data_missing);to_write = data.frame(vars = rownames(to_write),to_write)
file_to_write=file.path(outdir,paste("table_data_missing.txt",sep=""))
write_tsv(to_write,file_to_write);rm(to_write,file_to_write)


#*****************************************************************************#
# Standardise and remove outliers for cortical volume outside SD of 4 ----
#*****************************************************************************#
#* record number of participants with regional volume with SD >4
ctx_vars.ind = str_detect(names(df),"_volume");sum(ctx_vars.ind)
ctx_vars.ind = ctx_vars.ind | str_detect(names(df),"_thickness");sum(ctx_vars.ind)
ctx_vars.ind = ctx_vars.ind | str_detect(names(df),"_area");sum(ctx_vars.ind)
ctx_vars = names(df)[ctx_vars.ind]
df_tmp = df;count_4SD = c();for(xi in ctx_vars) {
  x =  df_tmp[[xi]]
  ind.4SD = abs(scale(x)[,1])>4
  count_4SD = c(count_4SD,sum(ind.4SD, na.rm=T))
  x[!is.na(ind.4SD) & ind.4SD] <- NA
  df_tmp[[xi]] <- x
}

df.count_4SD = data.frame(vars=ctx_vars,count_4sd = count_4SD)
write_tsv(df.count_4SD,file=file.path(outdir,"df.counts_4sd.txt"))
df = df_tmp;rm(df_tmp)

#*****************************************************************************#
# Summarise neuroimaging data (moved to 2_Descriptive_statistics) ----
#*****************************************************************************#
#* plot histogram
n = nrow(df)
if(!any("BMI"==missingCol)){
  phist_BMI <- plot_hist(dat = df, var = "BMI", col = "#D81B60")
  
}else{
  phist_BMI <- NULL
}
if(!any("waist"==missingCol)){
  phist_waist <- plot_hist(dat = df, var = "waist", col = "#FFC107")
}else{
  phist_waist <- NULL
}
if(!any("WHR"==missingCol)){
  phist_WHR <- plot_hist(dat = df, var = "WHR", col = "#004D40")
}else{
  phist_WHR <- NULL
}
#* arrange plots
p <- ggarrange(phist_BMI,phist_waist,phist_WHR, nrow=1)
#* annotate
p <- annotate_figure(p, 
                     top = text_grob(paste0(opt$cohort_name, " (", group, ")\nN = ",n),
                                     color = "black", 
                                     face = "bold", 
                                     size = 14))
p
#* save plot
p_hist_file=file.path(outdir,"histograms.jpg")
ggsave(p_hist_file, 
       bg = "white",
       plot = p, 
       width = 17, 
       height = 7, 
       units = "cm", 
       dpi = 150)

#* plot heat map
#* suppress warnings because a heatmap always has missing values 
suppressWarnings({
  vnames = unlist(str_split('BMI,waist,WHR,height,total_volume,total_area,mean_thickness,ICV,age_adiposity,age_mri',','))
  vnames = setdiff(vnames, missingCol)
  heat = plot_heatmap(dat = df %>% select(all_of(vnames)),
                      axisNames = vnames)
  
  # annotate
  heat <- annotate_figure(heat, 
                          top = text_grob(paste0(opt$cohort_name, " (", group, ")\nN = ",n),
                                          color = "black", 
                                          face = "bold", 
                                          size = 14))
  
  # save plot
  ggsave(file.path(outdir,"corrplot.jpg"), 
         bg = "white",
         plot = heat, 
         width = 15, 
         height = 15, 
         units = "cm", 
         dpi = 150)
})

#* 
desc_All = df %>% select(-ID,-FID) %>% psych::describe(IQR = T,quant=c(0.25,0.75))
desc_FemalesMales = df %>% select(-ID,-FID) %>% 
  mutate(sex = as.character(sex)) %>% 
  psych::describeBy(group="sex",IQR=T,quant=c(0.25,0.75))
save(desc_All,desc_FemalesMales,file=file.path(outdir,"descriptive_tables.RData"))

# categorical variables ----
descriptive_stats_discrete <- function(x,data){
  ret = table(data[[x]])
  ret = data.frame(ret,prop.table(ret))
  ret = subset(ret,select=-Var1.1)
  names(ret) = c('Value',"Count","Proportion")
  ret
}#descriptive_stats_discrete

cat_vars = attr(unlist(lapply(df, class)),'names')[which(unlist(lapply(df, class))%in%c("factor","character"))]
cat_vars = setdiff(cat_vars,"ID")

desc_All_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(df,select = cat_vars))
names(desc_All_discrete) = cat_vars

desc_Female_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(df,select = cat_vars, sex=="F"))
names(desc_Female_discrete) = cat_vars

desc_Male_discrete = lapply(cat_vars,descriptive_stats_discrete,data=subset(df,select = cat_vars, sex=="M"))
names(desc_Male_discrete) = cat_vars

save(desc_All_discrete,desc_Female_discrete,desc_Male_discrete,
     file=file.path(outdir,"descriptive_discrete_tables.RData"))

# correlation matrices and their plots for 'all' variables ----
corm = cor(subset(df,select=setdiff(names(df),c("ID",cat_vars))),use="p")
cormF = cor(subset(df,sex=="F",select=setdiff(names(df),c("ID",cat_vars))),use="p")
cormM = cor(subset(df,sex=="M",select=setdiff(names(df),c("ID",cat_vars))),use="p")
cor_list = list(corm=corm, cormF=cormF, cormM=cormM)
save(cor_list,file=file.path(outdir,'correlation_matrices_continuous_variables.Rdata'))

# scatter plots with correlation ----------------------------------------------
dA=subset(df,select=c( setdiff(names(non.neuro),c('ID',"FID",missingCol)), "ICV"))
png(file.path(outdir,"correlation_pairwise_scatter.png"),
    width=22,height=22,units="in",res=300)
suppressWarnings({
  pairs.panels(data.frame(dA,stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#BEBEBE33",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)})
dev.off()

png(file.path(outdir,"correlation_pairwise_scatter_Female.png"),
    width=22,height=22,units="in",res=300)
suppressWarnings({
  pairs.panels(data.frame(subset(dA,sex=="F"),stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#8B000033",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)})
dev.off()

png(file.path(outdir,"correlation_pairwise_scatter_Male.png"),
    width=22,height=22,units="in",res=300)
suppressWarnings({
  pairs.panels(data.frame(subset(dA,sex=="M"),stringsAsFactors = F), 
             method = "pearson", # correlation method
             hist.col = "#00008B33",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)})
dev.off()

# characteristic tables -------------------------------------------------------
# 
tab_vars = setdiff(c(names(df)),c("ID","FID"))
tableOne_All <- CreateTableOne(vars = tab_vars,
                               #strata = "trt", 
                               data = df, 
                               factorVars = cat_vars)

tableOne_sex_str <- CreateTableOne(vars = tab_vars,
                                   strata = "sex", 
                                   data = df, 
                                   factorVars = cat_vars)

if(length(table(df$is.underweight, useNA = "ifany"))>1){
  tableOne_bmi_str <- CreateTableOne(vars = tab_vars,
                                     strata = "is.underweight", 
                                     data = df, 
                                     factorVars = cat_vars)
  save(tableOne_bmi_str,
       file=file.path(outdir,
                      "characteristic_table_bmi_stratified.RData"))
}else{
  message("2_Descriptive_statistics: Only one bmi-group exisits.\n")
}

if(length(table(df$age.group, useNA = "ifany"))>1){
  tableOne_age_str <- CreateTableOne(vars = tab_vars,
                                     strata = "age_group", 
                                     data = df, 
                                     factorVars = cat_vars)
  save(tableOne_age_str,
       file=file.path(outdir,"characteristic_table_age_stratified.RData"))
}else{
  message("2_Descriptive_statistics: Only one age-group exisits.\n")
}

save(tableOne_All,
     file=file.path(outdir,
                    "characteristic_table.RData"))
save(tableOne_sex_str,
     file=file.path(outdir,
                    "characteristic_table_sex_stratified.RData"))
# ending ----
ending.message = paste0("\n2_Descriptive_statistics: Finishing obtaining descriptive statistics.\n")
message(ending.message,sep='')
cat(ending.message,sep='',file=logfile,append=T)

options(op)
