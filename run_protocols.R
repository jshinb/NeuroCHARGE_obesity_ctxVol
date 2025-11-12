#*****************************************************************************#
# Load inputs and set working directories ----
#*****************************************************************************#

cohort_input_file = #full path "~/Library/CloudStorage/OneDrive-Personal/neuroCHARGE_obesity_ctxVol_lipids/distribute/cohort_specific_inputs_testing_hypothetical.R"
source(cohort_input_file)
setwd(wd)

opt = data.frame(
  dir = wd,
  input_neuro = input_neuro,# required: name of the file containing the brain variables:
  input_non.neuro = input_non.neuro,# required: name of the file containing the other variables
  cohort_name = cohort_short_name,
  NAsymbol = missing_value_code
)# required:

outdir0=file.path(wd,paste0("outputs_to_send_",opt$cohort_name))
dir.create(outdir0)
file.copy(cohort_input_file ,outdir0,overwrite = T)

# define functions ----
cat("Prep: sourcing functions\n")
source("1_Pheno_preparation.R")

# The next steps will be done in groups (pooled or stratified) ----

# Primary analyses: pooled analyses of non-underweight undividuals at midlife
#* [Pooled analyses of non-underweight individuals at midlife]
pooled_analysis(data=df_all %>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"),
                cov_cohort_specific = cov_cohort_specific, 
                cov_additional = cov_additional, 
                analysis.type = "primary")

# Secondary analyses
#* [APOE4-stratified analyses of non-underweight individuals at midlife]
str_analysis(data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"), 
             str.var = "APOE4",str.var.name =  "APOE4.present",
             cov_cohort_specific=cov_cohort_specific,
             cov_additional = cov_additional)

#* [Age-stratified analyses of non-underweight individuals]
str_analysis(data=df_all %>% filter(is.underweight=="No"), 
             str.var = "Age",str.var.name =  "age.group",
             cov_cohort_specific=cov_cohort_specific,
             cov_additional = cov_additional)

#* [BMI-stratified analyses of individuals at midlife]
str_analysis(data = df_all%>% filter(age.group=="(34,65]_(34,65]"), 
             str.var = "BMI",str.var.name = "is.underweight",
             cov_cohort_specific=cov_cohort_specific,
             cov_additional = cov_additional)

#* [Ethnicity-stratified analyses of non-underweight individuals at midlife] # 
str_analysis(data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"), 
             str.var = "Ethnicity",str.var.name = "ethnicity",
             cov_cohort_specific=cov_cohort_specific,
             cov_additional = cov_additional)

#* [Ancestry-stratified analyses of non-underweight individuals at midlife] # 
str_analysis(data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"), 
             str.var = "Ancestry",str.var.name = "ancestry",
             cov_cohort_specific=cov_cohort_specific,
             cov_additional = cov_additional)

#* [Pooled analyses of non-underweight individuals at midlife adjusting for ethnicity] # 
pooled_analysis(data = df_all%>% filter(is.underweight=="No",age.group=="(34,65]_(34,65]"),
                cov_cohort_specific = c(cov_cohort_specific,"ethnicity"), 
                cov_additional = cov_additional,analysis.type = "adjEthnicity")

#* [Pooled analyses of all individuals]
pooled_analysis(data=df_all,cov_cohort_specific = cov_cohort_specific, 
                cov_additional = cov_additional, 
                analysis.type = "primary_all")
