library(ggplot2)
library(dplyr)
library(lme4)
library(effects)
library(effectsize)
library(tidyverse) 
library(RColorBrewer) 
library(lmerTest)
library(margins)
library(vctrs)
library(sjPlot)
library(sjmisc)


setwd('D:/[seted_path]')
set.seed(1)
dt = read.csv('tt_comment_case1_formlm.csv')

dt = dt %>% mutate(cm_knowledge = (cm_knowledge - min(cm_knowledge)) / (max(cm_knowledge) - min(cm_knowledge)))
dt <- within(dt, treatment <- relevel(as.factor(treatment), ref='Control'))
dt <- within(dt, assignment <- relevel(as.factor(assignment), ref='Control'))
dt = dt %>% mutate(prepost_b_da = ifelse(prepost=='Before', 'Pre', 'Post'),
                   prepost_bd_a = ifelse(prepost=='After', 'Post', 'Pre'))
dt <- within(dt, prepost_b_da <- relevel(as.factor(prepost_b_da), ref='Pre'))
dt <- within(dt, prepost_bd_a <- relevel(as.factor(prepost_bd_a), ref='Pre'))


vary = 'cm_knowledge_bi'
varx = 'treatment'
varm = 'sessions_bi'
vary_l = c('cm_knowledge','cm_knowledge_bi', 'cm_express_appreciation',
           'cm_knowledge_agreement', 'cm_knowledge_application', 
           'cm_knowledge_clarification', 'cm_knowledge_disagreement',
           'cm_knowledge_reconceptualization', 'cm_knowledge_reflection',
           'cm_relevant_to_video', 'cm_seek_professional_help', 'cm_seek_wb_info',
           'cm_self_disclose_mh', 'cm_wb_coping')



get_mlm_base_res <- function(vary, varx, lm_res_df, dt, need_two_level=F){
  col_select_l = c(vary, varx, 'prepost_b_da', 'TTID', 'VID')
  df = dt %>% dplyr::select(all_of(col_select_l))
  colnames(df) = c('y', 'intervention',  'prepost_b_da', 'TTID', 'VID')
  
  ## three level
  model = lmer(formula = y ~ 1 + intervention*prepost_b_da + (1|TTID) + (1|TTID:VID), data = df)
  #plot_model(model, type = "int", terms = c('intervention', 'prepost_b_da'))
  
  #margin = margins(model, variables = "intervention", at=list(prepost_b_da = c('Pre', 'Post')))
  margin = margins(model, variables = "prepost_b_da", at=list(intervention = c('Control', 'C+M', 'MO')))
  
  std_m = effectsize::standardize_parameters(model)
  lm_df = as.data.frame(summary(model)$coefficients)
  lm_df$var = row.names(lm_df)
  lm_df = cbind(lm_df, data.frame(std_m)[,c('Std_Coefficient', 'CI_low', 'CI_high')])
  lm_df['std_SE'] = attr(std_m, 'standard_error')
  
  lm_df$between_var = as.data.frame(VarCorr(model))['vcov'][1,]
  lm_df$within_var = as.data.frame(VarCorr(model))['vcov'][2,]
  lm_df$between_var_std = as.data.frame(VarCorr(model))['sdcor'][1,]
  lm_df$within_var_std = as.data.frame(VarCorr(model))['sdcor'][2,]
  lm_df$n_group = summary(model)$ngrps[1]
  lm_df$n_all = summary(model)$devcomp$dims['N'][1]
  # also calculate vpc: variance partitioning coefficient
  lm_df$vpc = as.data.frame(VarCorr(model))['vcov'][1,] / (as.data.frame(VarCorr(model))['vcov'][1,]+as.data.frame(VarCorr(model))['vcov'][2,])
  lm_df$type = 'comment|video|creator'
  lm_df$x = varx
  lm_df$y = vary
  lm_df$model_spe = 'base'
  
  ## add the marginal effect
  mm_df = data.frame(summary(margin))
  mm_df[setdiff(names(lm_df), names(mm_df))] <- NA
  lm_df[setdiff(names(mm_df), names(lm_df))] <- NA
  lm_df = rbind(lm_df, mm_df)
  
  lm_res_df = rbind(lm_res_df, lm_df)
  
  ## two level
  if(need_two_level){
    
    model = lmer(formula = y ~ 1 + intervention*prepost_b_da + (1 |TTID), data = df)
    
    std_m = effectsize::standardize_parameters(model)
    lm_df = as.data.frame(summary(model)$coefficients)
    lm_df = cbind(lm_df, data.frame(std_m)[,c('Std_Coefficient', 'CI_low', 'CI_high')])
    lm_df['std_SE'] = attr(std_m, 'standard_error')
    lm_df$var = row.names(lm_df)
    lm_df$between_var = as.data.frame(VarCorr(model))['vcov'][1,]
    lm_df$within_var = as.data.frame(VarCorr(model))['vcov'][2,]
    lm_df$between_var_std = as.data.frame(VarCorr(model))['sdcor'][1,]
    lm_df$within_var_std = as.data.frame(VarCorr(model))['sdcor'][2,]
    lm_df$n_group = summary(model)$ngrps[1]
    lm_df$n_all = summary(model)$devcomp$dims['N'][1]
    # also calculate vpc: variance partitioning coefficient
    lm_df$vpc = as.data.frame(VarCorr(model))['vcov'][1,] / (as.data.frame(VarCorr(model))['vcov'][1,]+as.data.frame(VarCorr(model))['vcov'][2,])
    lm_df$type = 'comment|creator'
    lm_df$x = varx
    lm_df$y = vary
    lm_df$model_spe = 'base'
    lm_res_df = rbind(lm_res_df, lm_df)
  }
  
  print(paste0(varx, '-', vary, ' done'))
  
  return(lm_res_df)
}


lm_res_df = data.frame()
for (vary in vary_l){
  for (varx in c('treatment')){
    lm_res_df = get_mlm_base_res(vary, varx, lm_res_df, dt)
  }
}
write.csv(lm_res_df, 'tt_case1_mlm_res_0824.csv', row.names = F, na='')


## version 2
dt_v2 = dt %>% filter(prepost != 'During')
lm_res_df2 = data.frame()
for (vary in vary_l){
  for (varx in c('treatment')){
    lm_res_df2 = get_mlm_base_res(vary, varx, lm_res_df2, dt_v2)
  }
}
write.csv(lm_res_df2, 'tt_case1_mlm_res_aprfuzzy_0824.csv', row.names = F, na='')


## engagement
dt = dt %>% mutate(engagement_large = ifelse(engagement_bi=='big_than10', 'Yes', 'No'))

## LARGE FOLLOWER
dt = dt %>% mutate(largefollow = ifelse(ttfollow>750000, 'Yes', 'No'))
dt <- within(dt, largefollow <- as.factor(largefollow))
dt = dt %>% mutate(ttf = as.integer(ttfollow/500000))
dt <- within(dt, ttf <- as.factor(ttf))
unique(dt$ttf)


## licensed, coaching
dt <- within(dt, licensed <- relevel(as.factor(licensed), ref='No'))
dt <- within(dt, coaching <- relevel(as.factor(coaching), ref='No'))
dt <- within(dt, engagement_large <- relevel(as.factor(engagement_large), ref='No'))
dt <- within(dt, largefollow  <- relevel(as.factor(largefollow), ref='No'))


set.seed(1)
get_mlm_moderation_res <- function(vary, varx, varm, lm_res_df, dt, need_two_level=F){
  col_select_l = c(vary, varx, varm,'prepost_b_da', 'TTID', 'VID')
  df = dt %>% dplyr::select(all_of(col_select_l))
  colnames(df) = c('y', 'intervention', 'moderation',  'prepost_b_da', 'TTID', 'VID')
  
  ## three level
  model = lmer(formula = y ~ 1 + intervention*prepost_b_da*moderation + (1 |TTID) + (1|TTID:VID), data = df)
  margin_no = margins(model, variables = "prepost_b_da", at=list(intervention = c('Control', 'C+M', 'MO'),
                                                                 moderation = c('No')))
  margin_yes = margins(model, variables = "prepost_b_da", at=list(intervention = c('Control', 'C+M', 'MO'),
                                                                  moderation = c('Yes')))
  
  std_m = effectsize::standardize_parameters(model)
  lm_df = as.data.frame(summary(model)$coefficients)
  lm_df = cbind(lm_df, data.frame(std_m)[,c('Std_Coefficient', 'CI_low', 'CI_high')])
  lm_df['std_SE'] = attr(std_m, 'standard_error')
  lm_df$var = row.names(lm_df)
  lm_df$between_var = as.data.frame(VarCorr(model))['vcov'][1,]
  lm_df$within_var = as.data.frame(VarCorr(model))['vcov'][2,]
  lm_df$between_var_std = as.data.frame(VarCorr(model))['sdcor'][1,]
  lm_df$within_var_std = as.data.frame(VarCorr(model))['sdcor'][2,]
  lm_df$n_group = summary(model)$ngrps[1]
  lm_df$n_all = summary(model)$devcomp$dims['N'][1]
  # also calculate vpc: variance partitioning coefficient
  lm_df$vpc = as.data.frame(VarCorr(model))['vcov'][1,] / (as.data.frame(VarCorr(model))['vcov'][1,]+as.data.frame(VarCorr(model))['vcov'][2,])
  lm_df$type = 'comment|video|creator'
  lm_df$x = varx
  lm_df$y = vary
  lm_df$m = varm
  lm_df$model_spe = 'moderation'
  
  ## add the marginal effect
  mm_df = data.frame(summary(margin_no))
  colnames(mm_df) = c('factor_no', 'intervention_no', 'moderation_no', 'ame_no',
                      'SE_no', 'z_no', 'p_no', 'lower_no', 'upper_no')
  mm_df[setdiff(names(lm_df), names(mm_df))] <- NA
  lm_df[setdiff(names(mm_df), names(lm_df))] <- NA
  lm_df = rbind(lm_df, mm_df)
  
  mm_df = data.frame(summary(margin_yes))
  colnames(mm_df) = c('factor_yes', 'intervention_yes', 'moderation_yes', 'ame_yes',
                      'SE_yes', 'z_yes', 'p_yes', 'lower_yes', 'upper_yes')
  mm_df[setdiff(names(lm_df), names(mm_df))] <- NA
  lm_df[setdiff(names(mm_df), names(lm_df))] <- NA
  lm_df = rbind(lm_df, mm_df)
  
  lm_res_df = rbind(lm_res_df, lm_df)
  
  
  ## two level
  if(need_two_level){
    
    model = lmer(formula = y ~ 1 + intervention*prepost_b_da*moderation + (1 |TTID), data = df)
    
    std_m = effectsize::standardize_parameters(model)
    lm_df = as.data.frame(summary(model)$coefficients)
    lm_df = cbind(lm_df, data.frame(std_m)[,c('Std_Coefficient', 'CI_low', 'CI_high')])
    lm_df['std_SE'] = attr(std_m, 'standard_error')
    lm_df$var = row.names(lm_df)
    lm_df$between_var = as.data.frame(VarCorr(model))['vcov'][1,]
    lm_df$within_var = as.data.frame(VarCorr(model))['vcov'][2,]
    lm_df$between_var_std = as.data.frame(VarCorr(model))['sdcor'][1,]
    lm_df$within_var_std = as.data.frame(VarCorr(model))['sdcor'][2,]
    lm_df$n_group = summary(model)$ngrps[1]
    lm_df$n_all = summary(model)$devcomp$dims['N'][1]
    # also calculate vpc: variance partitioning coefficient
    lm_df$vpc = as.data.frame(VarCorr(model))['vcov'][1,] / (as.data.frame(VarCorr(model))['vcov'][1,]+as.data.frame(VarCorr(model))['vcov'][2,])
    lm_df$type = 'comment|creator'
    lm_df$x = varx
    lm_df$y = vary
    lm_df$m = varm
    lm_df$model_spe = 'moderation'
    lm_res_df = rbind(lm_res_df, lm_df)
  }
  print(paste0(varx, '-', vary,'-',varm, ' done'))
  
  return(lm_res_df)
}

lm_res_df2 = data.frame()
for (vary in vary_l){
  for (varx in c('treatment')){
    for (varm in c('largefollow', 'engagement_large', 'lgbtq', 'licensed', 'coaching')){
      lm_res_df2 = get_mlm_moderation_res(vary, varx, varm, lm_res_df2, dt)
    }
  }
  write.csv(lm_res_df2, 'tt_case1_mlm_res_moderation_tmpsave_0824.csv', row.names = F, na='')
}
write.csv(lm_res_df2, 'tt_case1_mlm_res_moderation_0824.csv', row.names = F, na='')
