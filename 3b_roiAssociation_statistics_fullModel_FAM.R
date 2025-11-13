#*****************************************************************************#
#
# 3b: Examine associations between adiposoty vs. ctx structures across the 34 regions
# adjusting for age, sex, ICV (for volume), and genotype PCs, additional 
# covariates and any cohort-specific covariates.
#
#*****************************************************************************#
op <- options(warn=1)
logfile=file.path(outdir,"3b_roiAssociation_statistics_fullModel_FAM.log"); 
if(file.exists(logfile)) file.remove(logfile)

starting.message=paste0("\n3b. Start obtaining association estimates for adiposity vs. roi-ctx for the full model (family data).\n")
cat(starting.message)
cat(starting.message,file=logfile,append=F)

# starting ----
roi.34 = fread("freesurfer_34rois.txt",header=F) %>% pull(V1)

##* identify columns with all missing or with SD=0 (no variation)
missing <- colSums(is.na(df))
missingCol <- names(df)[missing == nrow(df)]
col.SD = apply(df %>% select(where(is.numeric)),2,sd,na.rm=T)
mono = apply( df %>% select(where(is.character)), 2, function(x){ all(length(unique(na.omit(x))) == 1) } )
monoVar = attr(mono,"names")[mono]
rm.var = c(names(col.SD)[col.SD==0],monoVar,missingCol)
cat('\nRemoved variables due to 100% missing or no variations: ', paste0(rm.var,collapse = ","),".","\n",sep='',
    file=logfile,append=T)

##* remove the identified columns from the analyses
cov_ctxV = setdiff(c('age_mri','age_diff','sex','ICV',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_ctxSA = setdiff(c('age_mri','age_diff','sex','ICV',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_ctxTH = setdiff(c('age_mri','age_diff','sex',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_BMI = setdiff(c('age_adiposity','age_diff','sex',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_WC = setdiff(c('age_adiposity','age_diff','sex','height',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)
cov_WHR = setdiff(c('age_adiposity','age_diff','sex',cov_genoPC,cov_cohort_specific,cov_additional),rm.var)

# adjust Ctx outcomes ----
adj_ctxV <- adj_ctxSA <- adj_ctxTH <-df %>% select(ID)
for (roii in c(roi.34,'global')){ 
  adj_ctxV = adj_ctxV %>% 
    left_join(get_adj_ctxV(
      roii=roii,
      d=df%>% rename(global_volume = total_volume),
      covariate_names = cov_ctxV,remove.ID=F,logfile=logfile),join_by(ID))
  adj_ctxSA = adj_ctxSA %>% 
    left_join(get_adj_ctxSA(
      roii=roii,d=df%>% rename(global_area = total_area),
      covariate_names = cov_ctxSA,remove.ID=F,logfile=logfile),join_by(ID))
  adj_ctxTH = adj_ctxTH %>% 
    left_join(get_adj_ctxTH(
      roii=roii,d=df%>% rename(global_thickness = mean_thickness),
      covariate_names = cov_ctxTH,remove.ID=F,logfile=logfile),join_by(ID))
}

# adjust adiposity measures ----
adj_BMI = get_adj_BMI(d=df,covariate_names = cov_BMI,remove.ID = F,logfile=logfile)
adj_WC = get_adj_WC(d=df,covariate_names = cov_WC,remove.ID = F,logfile=logfile)
adj_WHR = get_adj_WHR(d=df,covariate_names = cov_WHR,remove.ID = F,logfile=logfile)

# association tests ----
res_colnames = c("term", "Estimate", "Std..Error", "df", "Pr...t..", "N", 
  "sex", "ctx.pheno", "adiposity", "mod")
for(adiposity in c("BMI","WC","WHR")){
  
  if(adiposity == "BMI"){
    adj_adiposity = adj_BMI;cov_adiposity = cov_BMI
  }else if (adiposity == "WC") {
    adj_adiposity = adj_WC;cov_adiposity = cov_WC
  }else{
    adj_adiposity = adj_WHR;cov_adiposity = cov_WHR
  }
  
  for(pheno in c("volume","thickness","area")){
    if(pheno=='volume'){
      adj_ctx = adj_ctxV;cov_ctx = cov_ctxV
    }else if(pheno =='thickness'){
      adj_ctx = adj_ctxTH;cov_ctx = cov_ctxTH
    }else{
      adj_ctx = adj_ctxSA;cov_ctx = cov_ctxSA
    }
    names(adj_ctx) = str_remove(names(adj_ctx),paste0("_",pheno))
    
    #* approach 1 ----
    assoc_res_adj = c();for(roii in c(roi.34,'global') ){
      
      analdati = df %>% dplyr::select(ID,FID,sex) %>%
        left_join(adj_adiposity, join_by(ID)) %>% 
        left_join(subset(adj_ctx, select = c('ID',paste0(roii,".adjSex"),paste0(roii,".adj"))),
                  join_by(ID))
      
      analdati[['y']] <- inormal(analdati[[paste0(roii,".adj")]])
      analdati[['x']] <- inormal(analdati[[paste0(adiposity,".adj")]])
      
      # sex-specific adiposity  effect on cortical structure
      ## female
      if(max(table(analdati %>% filter(sex=="F") %>% pull(FID)))>1){
        mod.F = "y~x+(1|FID)"
        fit.assoc.F = lmer(as.formula(mod.F),data=analdati,subset=sex=="F", na.action = na.exclude)
      }else{
        mod.F = "y~x"
        fit.assoc.F = lm(as.formula(mod.F),data=analdati,subset=sex=="F", na.action = na.exclude)
        message("fit.assoc.F: will fit \'lm\' instead of \'lmer\'.")
      }
      excl_col = which(colnames(summary(fit.assoc.F)$coef)=="t value")
      assoc_res_adj_roii.F = summary(fit.assoc.F)$coef["x",-excl_col,drop=F]
      terms = rownames(assoc_res_adj_roii.F)
      rownames(assoc_res_adj_roii.F)=NULL
      assoc_res_adj_roii.F = data.frame(term=terms,assoc_res_adj_roii.F)
      assoc_res_adj_roii.F$N = nobs(fit.assoc.F)
      assoc_res_adj_roii.F$sex ="F"
      assoc_res_adj_roii.F$ctx.pheno = pheno
      assoc_res_adj_roii.F$adiposity = adiposity
      assoc_res_adj_roii.F$mod = mod.F
      if(class(fit.assoc.F)=="lm"){
        assoc_res_adj_roii.F$df = nobs(fit.assoc.F)
        assoc_res_adj_roii.F= subset(assoc_res_adj_roii.F,select=res_colnames)
      }
      ## male
      if(max(table(analdati %>% filter(sex=="M") %>% pull(FID)))>1){
        mod.M = "y~x+(1|FID)"
        fit.assoc.M = lmer(as.formula(mod.M),data=analdati,subset=sex=="M", na.action = na.exclude)
      }else{
        mod.M = "y~x"
        fit.assoc.M = lm(as.formula(mod.M),data=analdati,subset=sex=="M", na.action = na.exclude)
        message("fit.assoc.M: will fit \'lm\' instead of \'lmer\'.")
      }
      excl_col = which(colnames(summary(fit.assoc.M)$coef)=="t value")
      assoc_res_adj_roii.M = summary(fit.assoc.M)$coef["x",-excl_col,drop=F]
      terms = rownames(assoc_res_adj_roii.M)
      rownames(assoc_res_adj_roii.M)=NULL
      assoc_res_adj_roii.M = data.frame(term=terms,assoc_res_adj_roii.M)
      assoc_res_adj_roii.M$N = nobs(fit.assoc.M)
      assoc_res_adj_roii.M$sex ="M"
      assoc_res_adj_roii.M$ctx.pheno = pheno
      assoc_res_adj_roii.M$adiposity = adiposity
      assoc_res_adj_roii.M$mod = mod.M
      if(class(fit.assoc.M)=="lm"){
        assoc_res_adj_roii.M$df = nobs(fit.assoc.M)
        assoc_res_adj_roii.M= subset(assoc_res_adj_roii.M, select=res_colnames)
      }
      ## sex-combined and sex-by-interaction models
      if(max(table(analdati%>% pull(FID)))>1){
        ## sex-combined
        mod.all = "y~x+sex+(1|FID)"
        fit.assoc.sex_combined = lmer(as.formula(mod.all),data=analdati, na.action = na.exclude)
        ## sex-by-adiposity interaction
        mod.interaction = "y~sex*x+(1|FID)"
        fit.assoc = lmer(as.formula(mod.interaction),data=analdati,na.action = na.exclude)        
      }else{
        ## sex-combined
        mod.all = "y~x+sex"
        fit.assoc.sex_combined = lm(as.formula(mod.all),data=analdati, na.action = na.exclude)
        ## sex-by-adiposity interaction
        mod.interaction = "y~sex*x"
        fit.assoc = lmer(as.formula(mod.interaction),data=analdati,na.action = na.exclude)        
        message("fit.assoc.sex_combined, fit.assoc: will fit \'lm\' instead of \'lmer\'.")
      }
      ## sex-combined
      excl_col = which(colnames(summary(fit.assoc.sex_combined)$coef)=="t value")
      assoc_res_adj_roii.sex_combined = summary(fit.assoc.sex_combined)$coef["x",-excl_col,drop=F]
      terms = rownames(assoc_res_adj_roii.sex_combined)
      rownames(assoc_res_adj_roii.sex_combined)=NULL
      assoc_res_adj_roii.sex_combined = data.frame(term=terms,assoc_res_adj_roii.sex_combined)
      assoc_res_adj_roii.sex_combined$N = nobs(fit.assoc.sex_combined)
      assoc_res_adj_roii.sex_combined$sex ="sex-combined"
      assoc_res_adj_roii.sex_combined$ctx.pheno = pheno
      assoc_res_adj_roii.sex_combined$adiposity = adiposity
      assoc_res_adj_roii.sex_combined$mod = mod.all
      if(class(fit.assoc.sex_combined)=="lm"){
        assoc_res_adj_roii.sex_combined$df = nobs(fit.assoc.sex_combined)
        assoc_res_adj_roii.sex_combined=subset(assoc_res_adj_roii.sex_combined,
          select=res_colnames)
      }
      # sex-by-adiposity interaction effect
      excl_col = which(colnames(summary(fit.assoc)$coef)=="t value")
      assoc_res_adj_roii = summary(fit.assoc)$coef[,-excl_col,drop=F]
      terms = rownames(assoc_res_adj_roii)
      rownames(assoc_res_adj_roii)=NULL
      assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii)
      assoc_res_adj_roii$N = nobs(fit.assoc)
      assoc_res_adj_roii$sex = 'sex*adiposity interaction'
      assoc_res_adj_roii$ctx.pheno = pheno
      assoc_res_adj_roii$adiposity = adiposity
      assoc_res_adj_roii$mod = "y~sex*x+(1|FID)"
      if(class(fit.assoc)=="lm"){
        assoc_res_adj_roii$df = nobs(fit.assoc)
        assoc_res_adj_roii= subset(assoc_res_adj_roii,select=res_colnames)
        }
      
      # sex-combined adiposity-effect
      res.roii = rbind(assoc_res_adj_roii.F,
                       assoc_res_adj_roii.M,
                       assoc_res_adj_roii.sex_combined,
                       assoc_res_adj_roii)
      res.roii$roi = roii
      assoc_res_adj = rbind(assoc_res_adj,res.roii)
    }
    
    #* approach 2 ----
    cov_required = unique(c(cov_ctx,cov_adiposity))
    assoc_res_cov = c();for(roii in c(roi.34,'global')){
      yname = paste0(roii,".adj")
      xname = paste0(adiposity,".adj")
      select_vars=c('age.c','age.c2')
      if(pheno!="thickness"){#* [adjust SA and Vol for ICV]
        select_vars = c(select_vars,'icv.c','icv.c2')
      }
      ## data to fit
      analdati = subset(
        df %>% 
          mutate(age.c = scale(age_mri)[,1],
                 age.c2 = (scale(age_mri)/2)^2,
                 icv.c = scale(ICV)[,1],
                 icv.c2 = (scale(ICV)/2)^2),# for stability (Gelman)
        select=c("ID","FID",select_vars,
                 setdiff(cov_required,c("age_mri","ICV")))) %>% 
        left_join(subset(adj_ctx, select=c("ID",yname)),join_by(ID)) %>%
        left_join(subset(adj_adiposity, select=c("ID",xname)),join_by(ID))
      analdati[['y']] <- inormal(analdati[[yname]])
      analdati[['x']] <- inormal(analdati[[xname]])
      
      ## scale continuous covariates for better stability
      cont.covs = sapply(subset(df,select=setdiff(cov_required,c('age_mri','age_adiposity',"ICV",'sex'))), is.numeric)
      cont.covs = names(cont.covs)[cont.covs]
      analdati <- analdati %>% mutate(across(cont.covs, scale))
      
      ## define regression model
      mod0 = paste0(
        c('x',select_vars,
          setdiff(cov_required,c('age_mri','age_adiposity',"ICV",'sex'))),
        collapse = " + ")
      if(max(table(analdati%>% pull(FID)))>1){
        mod0 = paste0("y ~ sex*(",mod0,") + (1|FID)")
        fit.assoc = lmer(as.formula(mod0),data=analdati,na.action = na.exclude)
      }else{
        ## sex-by-adiposity interaction
        mod0 = paste0("y ~ sex*(",mod0,")")
        fit.assoc = lm(as.formula(mod0),data=analdati,na.action = na.exclude)
        message("fit.assoc (approach 2): will fit \'lm\' instead of \'lmer\'.")
      }
      assoc_res_adj_roii = summary(fit.assoc)$coef
      terms = rownames(assoc_res_adj_roii)
      excl_col = which(colnames(assoc_res_adj_roii)=="t value")
      assoc_res_adj_roii = data.frame(term=terms,assoc_res_adj_roii[,-excl_col])
      assoc_res_adj_roii = subset(assoc_res_adj_roii, term %in% c("x","sexM:x"))
      assoc_res_adj_roii$N = nobs(fit.assoc)
      assoc_res_adj_roii$roi = roii
      assoc_res_adj_roii$sex = 'sex*adiposoty interaction'
      assoc_res_adj_roii$ctx.pheno = pheno
      assoc_res_adj_roii$adiposity = adiposity
      assoc_res_adj_roii$mod = mod0
      if(class(fit.assoc)=="lm"){
        assoc_res_adj_roii$df = nobs(fit.assoc)
        assoc_res_adj_roii= subset(assoc_res_adj_roii,select=res_colnames)
      }
      rownames(assoc_res_adj_roii)=NULL
      assoc_res_cov = rbind(assoc_res_cov,assoc_res_adj_roii)
    }
    
    #* formatting results ----
    res1 = subset(assoc_res_adj,term %in% c("x","sexM:x"))
    res1 = res1 %>% rename(SE=Std..Error,P=Pr...t..)
    res1$method = 'm1'
    res1$term = str_remove(res1$term,"[.]adj")
    
    res2 = assoc_res_cov
    res2 = res2 %>% rename(SE=Std..Error,P=Pr...t..)
    res2$method = 'm2'
    
    #* creating outputs ----
    # table containing both sets of the results
    res_two_methods = res1 %>% filter(sex=="sex*adiposity interaction") %>%
      left_join(res2, join_by(roi,term,ctx.pheno,adiposity))
    
    roi_assoc_res_two_methods = res_two_methods
    save(res1,res2,roi_assoc_res_two_methods,
         file=file.path(outdir,
                        paste0("roi_assoc_res_two_methods_fullModel_",adiposity,"_ctx_",pheno,".Rdata"))
    )
    
    # plot
    beta.p = res_two_methods %>% 
      ggplot(aes(x=Estimate.x,y=Estimate.y,color=term)) +
      geom_point() + geom_abline(intercept = 0, slope = 1) + 
      theme_bw() +
      theme(legend.position = 'none')
    SE.p = res_two_methods %>% ggplot(aes(x=SE.x,y=SE.y,color=term)) +
      geom_point() + geom_abline(intercept = 0, slope = 1)+ 
      theme_bw() +
      theme(legend.position = 'none')
    Pval.p = res_two_methods %>% ggplot(aes(x=-log10(P.x),y=-log10(P.y),color=term)) +
      geom_point() +   theme_bw() +
      geom_abline(intercept = 0, slope = 1)
    p_method_comp = beta.p + SE.p + Pval.p
    p_method_comp +  plot_annotation(
      title = 'Comparing two approaches to adjusting for covariates',
      subtitle = 'Method 1 (x-axis): correct first -> examine associations vs. Method2 (y-axis): correct covariates while fitting.'
    )
    ggsave(file.path(outdir,
                     paste0('roi_assoc_res_two_methods_fullModel_',adiposity,"_ctx_",pheno,'.png')),
           width=8.5,height = 3, units="in", dpi=300)
    
    roi.34.ordered = subset(res1,sex=="sex-combined") %>% arrange(Estimate)
    roi.34.ordered = roi.34.ordered$roi
    res1 = res1 %>% mutate(roi=factor(roi,levels=rev(roi.34.ordered))) %>% 
      mutate(U95CI = Estimate + 1.96*SE,L95CI = Estimate - 1.96*SE) %>%
      arrange(roi)
    
    pos <- position_dodge(width=0.75)#
    p <- subset(res1,sex!="sex*adiposity interaction") %>% 
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,color=sex,group=sex)) + 
      geom_point(size=1, stroke=0.5, position=pos,alpha=0.5) +
      geom_errorbar(position=pos, width=0.25,size=0.5,alpha=0.5) + 
      geom_hline(yintercept=0, size=0.4) +
      geom_path(aes(group=sex),position=pos,alpha=0.5) +
      geom_hline(yintercept = 0)+
      xlab(NULL) + 
      coord_flip() +
      ylim((range(res1$L95CI,res1$U95CI)))+#only when coord_flip
      ggtitle(paste0(adiposity," vs. cortical ",pheno," associations") ) + 
      guides(col = guide_legend(reverse=T))
    
    p + theme_bw()
    ggsave(file.path(outdir,
                     paste0("forest_plot_assoc_",adiposity,"_roi_ctx_",pheno,"_fullModel.png")),
           width=6,height=6.5,units='in',dpi = 300)
  }
  #rm(adj_adiposity,adiposity)
}
rm(p,res1,res2,res_two_methods,adj_ctx,adj_ctxV,adj_ctxTH,adj_ctxSA,adj_BMI,adj_WC,adj_WHR)

#ending ----
ending.message = paste0("\nFinished obtaining association estimates for ",adiposity," vs. roi-ctx_",pheno," under full model in ",group,".\n")
cat(ending.message,sep='')
cat(ending.message,sep='',file=logfile,append=T)

options(op)
