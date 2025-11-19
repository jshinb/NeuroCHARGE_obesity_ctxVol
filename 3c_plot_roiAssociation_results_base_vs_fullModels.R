#*****************************************************************************#
#
# 3. compare the association results between base vs. full models
#
#*****************************************************************************#
op <- options(warn=1)

cat("\n3c. Loading association estimates for adiposity vs. roi-ctx for both the base and full models.\n")
logfile=file.path(outdir,"3c_plot_roiAssociation_results_base_vs_fullModels.log"); 
if(file.exists(logfile)) file.remove(logfile)

starting.message=paste0("\n3c. Start comparing association estimates between base vs. full models.\n")
cat(starting.message)
cat(starting.message,file=logfile,append=F)


# starting --------------------------------------------------------------------
for(adiposity in setdiff(c("BMI","waist","WHR"),rm.var)){
  for(pheno in c('volume','thickness','area')){
    baseres = paste0("roi_assoc_res_two_methods_baseModel_",adiposity,"_ctx_",pheno,".Rdata")
    fullres = paste0("roi_assoc_res_two_methods_fullModel_",adiposity,"_ctx_",pheno,".Rdata")
    load(file.path(outdir,baseres))
    res1_base = res1
    res2_base = res2 
    roi_assoc_res_two_methods_base = roi_assoc_res_two_methods
    rm(res1,res2,roi_assoc_res_two_methods)
    
    load(file.path(outdir,fullres))
    res1_full = res1
    res2_full = res2 
    roi_assoc_res_two_methods_full = roi_assoc_res_two_methods
    rm(res1,res2,roi_assoc_res_two_methods)
    
    # forest plot
    roi.34.ordered = subset(res1_base,sex=="sex-combined") %>% arrange(Estimate)
    roi.34.ordered = roi.34.ordered$roi
    res1 = c()
    res1 = rbind(res1,data.frame(res1_base,model="base"))
    res1 = rbind(res1,data.frame(res1_full,model="full"))
    
    res1 = res1 %>% 
      mutate(roi = factor(roi,levels=rev(roi.34.ordered)),
             model = factor(model, levels=rev(c('base','full'))),
             sex = factor(sex, levels=c('sex-combined','M','F')),
             U95CI = Estimate + 1.96*SE,
             L95CI = Estimate - 1.96*SE) %>% 
      arrange(roi)
    alpha.scale = 0.9
    pos <- position_dodge(width=0.75)#
    p <- subset(res1,sex!="sex*adiposity interaction") %>% 
      ggplot(aes(x=roi,y=Estimate,ymin=L95CI,ymax=U95CI,linetype=model,color=sex, shape=model)) + 
      geom_point(size=2, stroke=0.5, position=pos,alpha=alpha.scale) +
      geom_errorbar(position=pos, width=0.25,size=0.5,alpha=alpha.scale) + 
      geom_hline(yintercept=0, linewidth=0.4) +
      geom_path(aes(group=model),position=pos,alpha=alpha.scale) +
      scale_linetype_manual(values=c('base'="dashed",'full'='solid')) + 
      scale_shape_manual(values=c('base'=19,'full'=17)) + 
      scale_color_manual(values=c("sex-combined"="#619CFF","M"="#00BA38","F"="#F8766D")) +
      geom_hline(yintercept = 0)+
      xlab(NULL) + 
      facet_grid(cols = vars(sex)) + 
      coord_flip() +
      ylim((range(res1$L95CI,res1$U95CI)))+#only when coord_flip
      ggtitle(paste0(adiposity," vs. cortical ",pheno," associations")) + 
      theme_bw() +
      guides(col = guide_legend(reverse=F),
             linetype = guide_legend(reverse=F),
             shape = guide_legend(reverse=F))
    p 
    pfile = file.path(outdir,paste0("forest_plot_assoc_",adiposity,"_ctx_",pheno,"_base_fullModel.png"))
    ggsave(pfile,width=15,height=8,units = 'in',dpi=300)
  }
}

rm(res1_full,res2_full,res1_base,res2_base,
   roi_assoc_res_two_methods_base, roi_assoc_res_two_methods_full,
   p,pfile,pos,alpha.scale,
   roi.34.ordered,res1)#res under approach2 were not used for plotting

#ending ----
ending.message = paste0("\nFinished forest plots of association estimates for adiposity vs. roi-ctx for both the base and full models.\n")
cat(ending.message,sep='')
cat(ending.message,sep='',file=logfile,append=T)

options(op)
