### Functions written by Anna Furtjes, afurtjes@ed.ac.uk for any questions
### Distributed as part of our lifetime brain atrophy GWAS consortium effort
### Last edited May 2025


#' Inverse normal transformation (INT)
#'
#' https://www.biostars.org/p/80597/
#' See the supplement of Yang et al. Nature 2012. 
#'
#' @example
#' x1 <- 5:1
#' inormal(x1)
#'
#' x2 <- c(NA, 5:1, 5:1, NA)
#' inormal(x2)
inormal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

get_sums_volume <- function(data,ID_col){
  #data -> cortical structure data 
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  totalVol = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"volume",sep="_")
    Rh = paste("rh",roi,"volume",sep="_")
    totalVol_roi = (data[[Lh]] + data[[Rh]])
    totalVol = cbind(totalVol,totalVol_roi) #34 rois
  }
  L_Vol = 'lh_CortexVol_volume'
  R_Vol = 'rh_CortexVol_volume'
  global_totalVol = (data[[L_Vol]] + data[[R_Vol]])/2
  totalVol = data.frame(totalVol, global_totalVol)
  names(totalVol) = c(ID_col,paste(c(rois,"total"),"volume",sep="_"))
  # names(totalVol) = c(ID_col,rois,"total")#ID, roi-volume, and total volume
  totalVol
}

get_sums_area = function(data,ID_col){
  #data -> cortical structure
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  total_area = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"area",sep="_")
    Rh = paste("rh",roi,"area",sep="_")
    total_area_roi = (data[[Lh]] + data[[Rh]])
    total_area = cbind(total_area,total_area_roi) #34 rois
  }
  L_area = 'lh_WhiteSurfArea_area'
  R_area = 'rh_WhiteSurfArea_area'
  global_total_area = (data[[L_area]] + data[[R_area]])/2
  total_area = data.frame(total_area,global_total_area)
  names(total_area) = c(ID_col,paste(c(rois,"total"),"area",sep="_")) #ID, roi-area, and total area
  total_area
}

get_means_thickness <- function(data,ID_col){
  #data -> cortical thickness data 
  rois = fread('freesurfer_34rois.txt',header=F)$V1
  avgTH = subset(data,select=ID_col)
  for(roi in rois){
    Lh = paste("lh",roi,"thickness",sep="_")
    Rh = paste("rh",roi,"thickness",sep="_")
    avgTH_roi = (data[[Lh]] + data[[Rh]])/2
    avgTH = cbind(avgTH,avgTH_roi) #34 rois
  }
  L_TH = 'lh_MeanThickness_thickness'
  R_TH = 'rh_MeanThickness_thickness'
  global_avgTH = (data[[L_TH]] + data[[R_TH]])/2
  avgTH = data.frame(avgTH, global_avgTH)
  names(avgTH) = c(ID_col,paste(c(rois,"mean"),"thickness",sep="_"))
  avgTH
}

successivelyReduceAge = function(data = Share, ageCut = ageCut){
  # object to store results
  storeNames = c("Age cut-off value", "Cor", "p", "ci_l", "ci_u", "Measure")
  store = as.data.frame(matrix(nrow = length(ageCut)*3, ncol = length(storeNames)))
  names(store) = storeNames
  
  # scale data
  data$ratio = data$ratio*(-1)
  data$resid = data$resid*(-1)
  
  # iterate over each age cut-off and calculate scores
  for(i in ageCut){
    # store which age cut off iteration this is
    loc = which(is.na(store$`Age cut-off value`))[1]
    store[loc:(loc+2),"Age cut-off value"] = i
    
    # cut sample
    Youngdata = data[which(data$age <= i),]
    
    # calculate correlations
    ## difference
    store[loc,"Cor"] =
      with(Youngdata, cor.test(age, diff))$estimate
    
    store[loc,"p"] =
      with(Youngdata, cor.test(age, diff))$p.value
    
    store[loc,"ci_l"] =
      with(Youngdata, cor.test(age, diff))$conf.int[1]
    
    store[loc,"ci_u"] =
      with(Youngdata, cor.test(age, diff))$conf.int[2]
    
    store[loc,"Measure"] = "Difference score"
    
    ## ratio
    store[loc+1,"Cor"] =
      with(Youngdata, cor.test(age, ratio))$estimate
    
    store[loc+1,"p"] =
      with(Youngdata, cor.test(age, ratio))$p.value
    
    store[loc+1,"ci_l"] =
      with(Youngdata, cor.test(age, ratio))$conf.int[1]
    
    store[loc+1,"ci_u"] =
      with(Youngdata, cor.test(age, ratio))$conf.int[2]
    
    store[loc+1,"Measure"] = "Ratio score"
    
    ## residual
    store[loc+2,"Cor"] =
      with(Youngdata, cor.test(age, resid))$estimate
    
    store[loc+2,"p"] =
      with(Youngdata, cor.test(age, resid))$p.value
    
    store[loc+2,"ci_l"] =
      with(Youngdata, cor.test(age, resid))$conf.int[1]
    
    store[loc+2,"ci_u"] =
      with(Youngdata, cor.test(age, resid))$conf.int[2]
    
    store[loc+2,"Measure"] = "Residual score"
    
  }
  
  return(store)
  
}

agePlot = function(data = agePlotData){
  # work out nreaks for y-axis
  breaks = round(seq(from = minAge-1, to = maxAge+1, by = 2) )
  # if there are too many and they can't be read, reduce
  if(length(breaks) > 20){
    breaks = breaks[breaks %% 5 ==0]
  }
  
  p=ggplot(data = data)+
    geom_point(aes(x = Cor, y = `Age cut-off value`, colour = Measure), alpha = 0.5)+
    geom_errorbar(aes(y = `Age cut-off value`, xmin = ci_l, xmax = ci_u, colour = Measure), alpha = 0.3)+
    geom_vline(xintercept = 0, color = "grey")+
    xlab("Lifetime brain atrophy\ncorrelation with age")+
    ylab("Maximum age in subset\n(cut-off in age in years)")+
    scale_y_continuous(limits = c(minAge-1, maxAge+1), breaks = breaks)+
    #scale_x_continuous(limits = c(-0.1, 0.2), breaks = seq(from = -0.5, to = 0.35, by = 0.1))+
    scale_color_manual(values = c("#D81B60","#FFC107","#004D40"))+
    ggtitle(paste0(opt$cohort_name, " (N = ", nrow(neuro),")"))+
    theme_bw()+
    theme(panel.border = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(p)
}

plot_hist <- function(dat = neuro, var, col = "grey"){
  
  # make sure input data is data.frame
  dat = as.data.frame(dat)
  # rename for simplicity
  dat$var = dat[,var]
  
  # calculate mean from original variable
  varOr = str_remove(var, "_stand")
  mean = mean(dat[,varOr], na.rm=T)
  sd = sd(dat[,varOr], na.rm=T)
  
  # get title
  if(var == "total"){
    axisName = "total cortical volume"}else
    {
      axisName=var
    }
  
  # binwidth should differ dependeding on sample size
  if(nrow(dat) < 500){width = 0.5}
  if(nrow(dat) >=500){width = 0.1}
  
  # PLOT
  # different output when there is a "sample" column
  plot = ggplot(dat, aes(x = var))+
    geom_histogram(alpha = 0.5, fill = col) + #, binwidth = width)+
    xlab(axisName)+
    ylab("Count")+
    theme_bw()+
    #scale_x_continuous(sec.axis = sec_axis(trans=~.*sd+mean))+
    theme(panel.border = element_blank())
  
  return(plot)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

plot_heatmap = function(dat, axisNames)
{
  # define function to obtain lower triangle
  # get correlation matrix for both samples together
  cor = cor(dat, use = "p")
  cor = get_lower_tri(cor)
  # melt matrix
  melted = reshape2::melt(cor)
  # get rounded value
  melted$value_round = round(melted$value, digit = 2)
  melted$distance0 = abs(melted$value)
  
  # plot
  library(ggplot2)
  
  p = ggplot(data = melted)+
    geom_point(aes(x = Var1, y = Var2, shape = value, fill = value, size = distance0), 
               shape = 21, alpha = 0.7, colour = "white") +
    scale_fill_gradient2(low = "#82A0D8", high = "#8DDFCB", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab" ,name="Correlation", guide = "legend")+
    scale_size_continuous(range = c(1, 15), guide = "none")+
    geom_text(aes(Var1, Var2, label = value_round), color = "black", size = 4)+
    xlab("")+
    ylab("")+
    scale_x_discrete(labels = axisNames)+
    scale_y_discrete(labels = axisNames)+
    guides(fill = "none")+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.text.x = element_text(angle=-45,vjust = 0.5, hjust=0))
  
  return(p)
}

#* will think about the following for outlier detection/removal:
rm.sd4 = function(x){
  x.scaled = scale(x)[,1]
  print(range(x.scaled),na.rm=T)
  ind = !is.na(x.scaled)
  ind = ind & abs(x.scaled)>4
  x[ind] <- NA
  x
}

#* adjust brain outcomes ----
get_adj_ctxV = function(roii,d,covariate_names,remove.ID=T,logfile){
  #age, sex and ICV-corrected
  yname = paste(roii,"volume",sep="_")
  
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2) %>%
    mutate(icv.c = scale(ICV)[,1],icv.c2 = (scale(ICV)/2)^2)
  
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2','icv.c', 'icv.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2','icv.c','icv.c2',setdiff(covariate_names,c('age_mri',"ICV",'sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_volume.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_volume.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_volume.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_volume.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxV: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_volume.adjSex"),paste0(roii,"_volume.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_ctxTH = function(roii,d,covariate_names,remove.ID,logfile){
  #age and sex-corrected
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2)
  yname = paste(roii,"thickness",sep="_")
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2',setdiff(covariate_names,c('age_mri','sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_thickness.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_thickness.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_thickness.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_thickness.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxTH: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_thickness.adjSex"),paste0(roii,"_thickness.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_ctxSA = function(roii,d,covariate_names,remove.ID=T,logfile){
  #* [age and sex-corrected and **ICV**]
  yname = paste(roii,"area",sep="_")
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2) %>%
    mutate(icv.c = scale(ICV)[,1],icv.c2 = (scale(ICV)/2)^2)
  
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2','icv.c', 'icv.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2','icv.c','icv.c2',setdiff(covariate_names,c('age_mri',"ICV",'sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_area.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_area.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_area.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_area.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxSA: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_area.adjSex"),paste0(roii,"_area.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_ctxSA.woICV = function(roii,d,covariate_names,remove.ID=T,logfile){
  #* [age and sex-corrected - not used at the moment]
  #age and sex-corrected (but...what about ICV??? it is affected by ICV)
  d = d %>% 
    # for stability (Gelman)
    mutate(age = age_mri) %>%
    mutate(age.c = scale(age)[,1],age.c2 = (scale(age)/2)^2)
  yname = paste(roii,"area",sep="_")
  analdati = subset(d,select=c("ID",covariate_names,'age.c', 'age.c2',yname)) 
  analdati[['y']] = inormal(analdati[[yname]])
  
  analdati = analdati %>% 
    mutate(y = ifelse(abs(scale(y)[,1]) >4, NA, y))# remove outliers here
  mod0.roii = paste0(c('age.c','age.c2',setdiff(covariate_names,c('age_mri','sex'))),collapse = " + ")
  mod0.roii = paste("y ~ ",mod0.roii,sep="")
  
  # age, sex and other covariates: 
  fit0M.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="M"),na.action=na.exclude)
  fit0F.roii = lm(as.formula(mod0.roii),data=analdati %>% filter(sex=="F"),na.action=na.exclude)
  analdatiM =  analdati%>% filter(sex=="M") %>% select(ID)
  analdatiF =  analdati%>% filter(sex=="F") %>% select(ID)
  
  analdatiM[[paste0(roii,"_area.adjSex")]] = resid(fit0M.roii)
  analdatiM[[paste0(roii,"_area.adj")]] = resid(fit0M.roii)+coef(fit0M.roii)[1]#unadjusted for sex
  analdatiF[[paste0(roii,"_area.adjSex")]] = resid(fit0F.roii)
  analdatiF[[paste0(roii,"_area.adj")]] = resid(fit0F.roii)+coef(fit0F.roii)[1]#unadjusted for sex
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  
  capture.output(paste0("get_adj_ctxSA: ",mod0.roii),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste0(roii,"_area.adjSex"),paste0(roii,"_area.adj")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

#* adjust adiposity measures ----
get_adj_BMI = function(d,covariate_names,remove.ID,logfile){
  analdati = d %>%
    mutate(age = age_adiposity) %>%
    mutate(age.c = scale(age)[,1],)
  
  analdati[['x']] = inormal(analdati[["BMI"]])
  # remove outliers based on BMI
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c',setdiff(covariate_names,c('age_adiposity','sex'))),collapse = " + ")
  mod0 = paste("x ~ ",mod0,sep="")
  
  # age, sex and other covariates: 
  fit0M = lm(as.formula(mod0),data=analdati%>% filter(sex=="M"),na.action=na.exclude)
  fit0F = lm(as.formula(mod0),data=analdati%>% filter(sex=="F"),na.action=na.exclude)
  analdatiM = analdati %>% filter(sex=="M") %>% select(ID) %>% 
    mutate(BMI.adjSex = resid(fit0M)) %>%
    mutate(BMI.adj = resid(fit0M)+coef(fit0M)[1])
  analdatiF = analdati %>% filter(sex=="F") %>% select(ID) %>% 
    mutate(BMI.adjSex = resid(fit0F)) %>%
    mutate(BMI.adj = resid(fit0F)+coef(fit0F)[1])
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))

  capture.output(paste0("get_adj_BMI: ",mod0),file=logfile,append=T)
  ret = subset(analdati,select=c("ID",paste("BMI","adjSex",sep="."),paste("BMI","adj",sep=".")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_WC = function(d,covariate_names,remove.ID=F,logfile){
  analdati = d %>%
    mutate(age = age_adiposity) %>%
    mutate(age.c = scale(age)[,1],)
  
  analdati[['x']] = inormal(analdati[["waist"]])
  # remove outliers based on waist
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c','height',setdiff(covariate_names,c('age_adiposity','sex','height'))),collapse = " + ")
  mod0 = paste("x ~ ",mod0,sep="")
  
  # age, sex and other covariates: 
  fit0M = lm(as.formula(mod0),data=analdati%>% filter(sex=="M"),na.action=na.exclude)
  fit0F = lm(as.formula(mod0),data=analdati%>% filter(sex=="F"),na.action=na.exclude)
  analdatiM = analdati %>% filter(sex=="M") %>% select(ID) %>% 
    mutate(WC.adjSex = resid(fit0M)) %>%
    mutate(WC.adj = resid(fit0M)+coef(fit0M)[1])
  analdatiF = analdati %>% filter(sex=="F") %>% select(ID) %>% 
    mutate(WC.adjSex = resid(fit0F)) %>%
    mutate(WC.adj = resid(fit0F)+coef(fit0F)[1])
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  capture.output(paste0("get_adj_WC: ",mod0),file=logfile,append=T)
  
  ret = subset(analdati,select=c("ID",paste("WC","adjSex",sep="."),paste("WC","adj",sep=".")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

get_adj_WHR = function(d,covariate_names,remove.ID=F,logfile){
  analdati = d %>%
    mutate(age = age_adiposity) %>%
    mutate(age.c = scale(age)[,1],)
  
  analdati[['x']] = inormal(analdati[["WHR"]])
  # remove outliers based on waist
  analdati = analdati %>% mutate(y = ifelse(abs(scale(x)[,1]) >4, NA, x))
  mod0 = paste0(c('age.c',setdiff(covariate_names,c('age_adiposity','sex'))),collapse = " + ")
  mod0 = paste("x ~ ",mod0,sep="")
  
  # age, sex and other covariates: 
  fit0M = lm(as.formula(mod0),data=analdati%>% filter(sex=="M"),na.action=na.exclude)
  fit0F = lm(as.formula(mod0),data=analdati%>% filter(sex=="F"),na.action=na.exclude)
  
  analdatiM = analdati %>% filter(sex=="M") %>% select(ID) %>% 
    mutate(WHR.adjSex = resid(fit0M)) %>%
    mutate(WHR.adj = resid(fit0M)+coef(fit0M)[1])
  analdatiF = analdati %>% filter(sex=="F") %>% select(ID) %>% 
    mutate(WHR.adjSex = resid(fit0F)) %>%
    mutate(WHR.adj = resid(fit0F)+coef(fit0F)[1])
  analdati = analdati %>% left_join(rbind(analdatiM,analdatiF),join_by(ID))
  capture.output(paste0("get_adj_WHR: ",mod0),file=logfile,append=T)
  
  ret = subset(analdati,select=c("ID",paste("WHR","adjSex",sep="."),paste("WHR","adj",sep=".")))
  if(remove.ID) ret = ret %>% select(-ID)
  ret
}

#* check sample sizes for discrete covariste (if 1, cannot run regression adjusting for that var)----
check.sample.size <- function(anal_outdir,disc_cov_names){
  load(file.path(anal_outdir,"descriptive_discrete_tables.RData"))
  sm.F = desc_Female_discrete
  sm.M = desc_Male_discrete
  disc_cov_names = intersect(names(sm.F),disc_cov_names)
  
  if(length(disc_cov_names)>0){
    minFreq.F <- minFreq.M <- minLevel.F <- minLevel.M <- c(); for(i in disc_cov_names){
      min.ind = which.min(sm.F[[i]] %>% pull(Count))
      minFreq.F <- c(minFreq.F, sm.F[[i]][["Count"]][min.ind])
      minLevel.F <- c(minLevel.F,as.character(sm.F[[i]][['Value']])[min.ind])
      
      min.ind = which.min(sm.M[[i]] %>% pull(Count))
      minFreq.M <- c(minFreq.M, sm.M[[i]][["Count"]][min.ind])
      minLevel.M <- c(minLevel.M,as.character(sm.M[[i]][['Value']])[min.ind])
    }
    df_minFreq = data.frame(vars=disc_cov_names,minLevel.F,minFreq.F,minLevel.M,minFreq.M)
    lessthan_5 = min(minFreq.F,minFreq.M)<5
  }else{#length(disc_cov_names)=0
    df_minFreq = NULL
    lessthan_5 = NULL
  }
  ret = list(df_minFreq=df_minFreq,lessthan_5 = lessthan_5)
  ret
}
#* stratified or pooled analysis functions ----
str_analysis = function(data,str.var,str.var.name,cov_cohort_specific,cov_additional,is.family.data){
  groups = unique(na.omit(data[[str.var.name]]))
  if(length(groups)>1){
    # 2. descriptive stats
    for(group in groups){
      outdir = file.path(outdir0,paste0(str.var,"_",str.var.name,"_",group))
      dir.create(outdir)
      df = data %>% filter(.data[[str.var.name]] == group)
      source("2_Descriptive_statistics.R",local = T)
      # check sample sizes for each categorical variable:
      small.sample_base <- check.sample.size(outdir,cov_cohort_specific)
      if(!is.null(small.sample_base$lessthan_5)){
        # discrete cohort-specific covs exist and can do sex-specific adjustment --> proceed to fitting
        run_base = !small.sample_base$lessthan_5
      }else{#is.null(small.sample_base$lessthan_5
        run_base = TRUE# discrete cohort-specific covs don't exist --> proceed to fitting
      }
      if(run_base) {
        if(is.family.data){
          source("3a_roiAssociation_statistics_baseModel_FAM.R",local = T)
        }else{
          source("3a_roiAssociation_statistics_baseModel.R",local = T)
        }
      }else{
        msg = paste0('\nBefore 3a: ',str.var,"_",group," (n=",nrow(df),
                     '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                     'we will not run association tests. \n')
        message(msg)
        print(small.sample_base$df_minFreq)
        cat(msg,file=logfile,append = T)
        capture.output(print(small.sample_base$df_minFreq),file = logfile,append=T)
      }#else for if(run_base)
      small.sample_full <- check.sample.size(outdir,c(cohort_specific_covs,cov_additional))
      run_full = !small.sample_full$lessthan_5
      if(run_full){
        if(is.family.data){
          source("3b_roiAssociation_statistics_fullModel_FAM.R",local = T)
        }else{
          source("3b_roiAssociation_statistics_fullModel.R",local = T)
        }
        source("3c_plot_roiAssociation_results_base_vs_fullModels.R",local=T)
      }else{
        msg = paste0('\nBefore 3b: ',str.var,"_",group," (n=",nrow(df),
                     '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                     'we will not run association tests. \n')
        message(msg)
        print(small.sample_full$df_minFreq)
        cat(msg,file=logfile,append = T)
        capture.output(print(small.sample_full$df_minFreq),file = logfile,append=T)
      }#else for if(run_full)
    }#for(group in groups)
  }else{#length(groups)=1 or 0
    group = groups[1]
    outdir = file.path(outdir0,paste0(str.var,"_",str.var.name,"_",group))
    dir.create(outdir)
    logfile=file.path(outdir,"strtified_analyses.log")
    df = data
    msg = paste0("\n There is only one group:", str.var,"_",group," (n=",nrow(df),"): No need to run ",str.var,"-stratified analyses.\n")
    message(msg)
    cat(msg,file=logfile,append = T)
  }
}

pooled_analysis = function(
    data,
    cov_cohort_specific,
    cov_additional,
    analysis.type,# c("primary","primary_all",'adjEthnicity')
    is.family.data)
  {
  group = "pooled"
  # 2. descriptive stats
  outdir = file.path(outdir0,paste0(group,"_",analysis.type))
  dir.create(outdir)
  df = data
  if(str_detect(analysis.type, 'primary')){
    do.run = TRUE
  }else{
    do.run = length(unique(df[['ethnicity']])) >1 #run if more than 1 ethnicity
  }
  if(do.run){
    source("2_Descriptive_statistics.R",local = T)
    # check sample sizes for each categorical variable:
    small.sample_base <- check.sample.size(outdir,cov_cohort_specific)
    if(!is.null(small.sample_base$lessthan_5)){
      # discrete cohort-specific covs exist and can do sex-specific adjustment --> proceed to fitting
      run_base = !small.sample_base$lessthan_5
    }else{#is.null(small.sample_base$lessthan_5
      run_base = TRUE# discrete cohort-specific covs don't exist --> proceed to fitting
    }
    if(run_base) {
      if(is.family.data){
        source("3a_roiAssociation_statistics_baseModel_FAM.R",local = T)
      }else{
        source("3a_roiAssociation_statistics_baseModel.R",local = T)
      }
    }else{
      msg = paste0('\nBefore 3a: ',group," (n=",nrow(df),
                   '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                   'we will not run association tests. \n')
      message(msg)
      print(small.sample_base$df_minFreq)
      cat(msg,file=logfile,append = T)
      capture.output(print(small.sample_base$df_minFreq),file = logfile,append=T)
    }#else for if(run_base)
    small.sample_full <- check.sample.size(outdir,c(cohort_specific_covs,cov_additional))
    run_full = !small.sample_full$lessthan_5
    if(run_full){
      if(is.family.data){
        source("3b_roiAssociation_statistics_fullModel_FAM.R",local = T)
      }else{
        source("3b_roiAssociation_statistics_fullModel.R",local = T)
      }
      source("3c_plot_roiAssociation_results_base_vs_fullModels.R",local=T)
    }else{
      msg = paste0('\nBefore 3b: ',group," (n=",nrow(df),
                   '): Because sex-specific sample sizes for some discreate covariates are less than 5, ',
                   'we will not run association tests. \n')
      message(msg)
      print(small.sample_full$df_minFreq)
      cat(msg,file=logfile,append = T)
      capture.output(print(small.sample_full$df_minFreq),file = logfile,append=T)
    }#else for if(run_full)
  }else{#!do.run
    logfile=file.path(outdir,"pooled_analyses.log")
    msg = paste0("\n There is only one ethnic group, \'",unique(df$ethnicity),"\' (n=",nrow(df),"): We will not run ethnicity-adjusted analyses.\n")
    message(msg)
    capture.output(msg,file=logfile,append=T)
  }
}
