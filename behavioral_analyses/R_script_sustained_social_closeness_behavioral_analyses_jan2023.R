# R script for the behavioral analyses reported in the manuscript:
# "Connected in bad times and in good times: Empathy induces stable social closeness"

# author: Anne Saulin
# email: anne.saulin@yahoo.com

################################################################################
# load required packages

library(dplyr)
library(car)
library(psych)
library(lme4)
library(tidyr)
library(emmeans)
library(MASS)
library(tibble)
library(brms)
library(BayesFactor)



######## set working directory to where you downloaded the files to
# setwd("~/happy_data_place")

#################################################################################
## read in data

# questionnaire data
qa_data               <- read.csv("questionnaire_scores.csv")

# emotion_ratings
emotion_ratings_fMRI  <- read.csv("emotion_ratings_fMRI.csv")
emotion_ratings_repl  <- read.csv("emotion_ratings_repl.csv")
emotion_ratings_contr <- read.csv("emotion_ratings_contr.csv")

# closeness ratings, trait empathy (empathic concern & perspective-taking), and emotion ratings
behavior_fMRI_repl <- read.csv("behavior_fMRI_repl.csv")
behavior_contr     <- read.csv("behavior_contr.csv")

# closeness ratings only 
# (containing only the variables required for mixed models analyses for closeness ratings)
closeness_fMRI  <- read.csv("closeness_fMRI.csv")
closeness_repl  <- read.csv("closeness_repl.csv")
closeness_contr <- read.csv("closeness_contr.csv")

# dataframe for post-hoc analyses on extracted betas from left IFG/AI and bil STS/TPJ
behavior_and_neural_betas <- read.csv("behavior_and_neural_betas.csv")

##################################################################################
## analyses

###############################################
## fMRI study and replication study
###############################################

### differences in questionnaire scores across studies?

lm_PR_studies <- lm(unlist(PR) ~ study, qa_data)
summary(lm_PR_studies)
lm_PT_studies <- lm(unlist(PT) ~ study, qa_data)
summary(lm_PT_studies)
lm_EC_studies <- lm(EC.value ~ study, qa_data)
summary(lm_EC_studies)


### manipulation checks

## differential response to trial type (reinforced vs. non-reinforced) 


lmer_emotion_ratings_fMRI_repl <- lmer(emotion_rating ~ as.factor(trial_type)*as.factor(block_no)*condition*sample + (1 | id),
				behavior_fMRI_repl)
summary(lmer_emotion_ratings_fMRI_repl)
Anova(lmer_emotion_ratings_fMRI_repl, type = 3)

## emotional response to other's pain ~ trait empathy

# first aggregate the data
agg_emo_facet_emp_agg <- aggregate(emotion_rating ~ EC + PT + sample + trial_type + id, behavior_fMRI_repl, mean)%>%
mutate(., PT = scale(PT))%>%
pivot_longer(cols = 1:2, names_to= "trait", values_to = "score" )

# conduct mixed model analysis
lmer_emo_facet_emp_agg <- lmer(emotion_rating ~ trait*score*as.factor(trial_type)*sample+
(1|id),
agg_emo_facet_emp_agg)
Anova(lmer_emo_facet_emp_agg, type = 3)
summary(lmer_emo_facet_emp_agg)



############################################
## behavioral results

### fMRI study

# closeness over time (due to the mirrored presentation in the fMRI scanner, ratings have to be mirrored, to)
lmer_closeness_fMRI <- lmer(scale(closeness*-1+100) ~ condition*scale(trial_no)*as.factor(block_no) + 
                                                       (1|id) + (0 + scale(trial_no)|id), 
                                                       closeness_fMRI)
Anova(lmer_closeness_fMRI, type = 3)
summary(lmer_closeness_fMRI)

# use bayes to test for null efect of three-way interaction (block*trial number*condition)
pr = prior(normal(0, 1), class = 'b')
lmer_bayes_closeness_additive_fMRI <- brms::brm(scale(closeness*-1+100) ~ scale(trial_no) + condition*as.factor(block_no) + 
                                                  (1|id) + (0 + scale(trial_no)|id), 
                                                closeness_fMRI, prior = pr, cores = 4, save_pars = save_pars(all = TRUE), iter = 9000)

lmer_bayes_closeness_interaction_fMRI <- brms::brm(scale(closeness*-1+100) ~ scale(trial_no)*condition*as.factor(block_no) + 
                                                     (1|id) + (0 + scale(trial_no)|id), 
                                                   closeness_fMRI, prior = pr, cores = 4, save_pars = save_pars(all = TRUE), iter = 9000)

# interaction effect with trial_no?
bayesfactor_closeness_fMRI <- bayes_factor(lmer_bayes_closeness_interaction_fMRI,lmer_bayes_closeness_additive_fMRI)


# post-doc contrast
df_last5_last5_fMRI <- aggregate(closeness ~ condition + block_no + id, closeness_fMRI[closeness_fMRI$trial_no>=21,], mean)
t.test(df_last5_last5_fMRI$closeness[df_last5_last5_fMRI$condition=="decay" & df_last5_last5_fMRI$block_no=="2"],
       df_last5_last5_fMRI$closeness[df_last5_last5_fMRI$condition=="decay" & df_last5_last5_fMRI$block_no=="4"],
       paired = TRUE)
library(BayesFactor)
ttestBF(df_last5_last5_fMRI$closeness[df_last5_last5_fMRI$condition=="decay" & df_last5_last5_fMRI$block_no=="4"],
        df_last5_last5_fMRI$closeness[df_last5_last5_fMRI$condition=="decay" & df_last5_last5_fMRI$block_no=="2"],
        paired = TRUE) 

### replication study
lmer_closeness_repl <- lmer(scale(closeness) ~ condition*scale(trial_no)*as.factor(block_no) + 
                                                     (1|id) + (0 + scale(trial_no)|id), 
                                                     closeness_repl)
Anova(lmer_closeness_repl, type = 3)
summary(lmer_closeness_repl)

# bayes
lmer_bayes_closeness_additive_repl <- brms::brm(scale(closeness) ~ scale(trial_no) + condition*as.factor(block_no) + 
                                                  (1|id) + (0 + scale(trial_no)|id), 
                                                closeness_repl, prior = pr, cores = 4, save_pars = save_pars(all = TRUE), iter = 9000)

lmer_bayes_closeness_interaction_repl <- brms::brm(scale(closeness) ~ scale(trial_no)*condition*as.factor(block_no) + 
                                                     (1|id) + (0 + scale(trial_no)|id), 
                                                   closeness_repl, prior = pr, cores = 4, save_pars = save_pars(all = TRUE), iter = 9000)

# interaction effect with trial_no?
bayesfactor_closeness_repl <- bayes_factor(lmer_bayes_closeness_interaction_repl,lmer_bayes_closeness_additive_repl)


# compare last five trials of block 1 treatment("decay") and last five trials of block 2 of the treatment condition
agg_last5_last5_repl <- aggregate(closeness ~ condition + block_no + id, closeness_repl[closeness_repl$trial_no>=21,], mean)
t.test(agg_last5_last5_repl$closeness[agg_last5_last5_repl$condition=="empathy_decay" & agg_last5_last5_repl$block_no=="2"],
       agg_last5_last5_repl$closeness[agg_last5_last5_repl$condition=="empathy_decay" & agg_last5_last5_repl$block_no=="4"],
       paired = TRUE)

ttestBF(agg_last5_last5_repl$closeness[agg_last5_last5_repl$condition=="empathy_decay" & agg_last5_last5_repl$block_no=="2"],
        agg_last5_last5_repl$closeness[agg_last5_last5_repl$condition=="empathy_decay" & agg_last5_last5_repl$block_no=="4"],
        paired = TRUE)

###########################
## compare samples directly
lmer_closeness_fMRI_repl <- lmer(closeness ~ condition*scale(trial_no)*block_no*sample + 
                                           (1|id) + (0 + scale(trial_no)|id), 
                                           behavior_fMRI_repl)
Anova(lmer_closeness_fMRI_repl, type = 3)
summary(lmer_closeness_fMRI_repl)


#############################################
## control study (reciprocity motive)
#############################################

##### manipulation check
## emotion ratings sensitive to trial type?
lmer_emotion_ratings_contr <- lmer(scale(choice_val.mean) ~ trial_type*condition + 
                                                         (1|id), emotion_ratings_contr)
Anova(lmer_emotion_ratings_contr, type = 3)
summary(lmer_emotion_ratings_contr)

## emotion ratings influenced by trait positive reciprocity (PNR)?
agg_emo_PNR_cond_contr   <- aggregate(emotion_rating ~ pos_PNR + trial_type + id + condition + block_no, 
                                 behavior_contr, mean)
lmer_emo_PNR_cond_contr <- lmer(emotion_rating ~ pos_PNR*condition*as.factor(block_no)*as.factor(trial_type)+
                             (1|id),
                         agg_emo_PNR_cond_contr)
Anova(lmer_emo_PNR_cond_contr, type = 3)
summary(lmer_emo_PNR_cond_contr)

#### behavioral results

lmer_closeness_reciprocity <- lmer(scale(closeness) ~ condition*scale(trial_no)*as.factor(block_no) + 
                               (1|id) +  (0 + scale(trial_no)|id), 
                               closeness_contr)
Anova(lmer_closeness_reciprocity, type = 3)
summary(lmer_closeness_reciprocity)
emtrends(lmer_closeness_reciprocity, c("block_no", "condition"), var = "trial_no")

# bayes
lmer_bayes_closeness_additive_repl <- brms::brm(scale(closeness) ~ scale(trial_no) + condition*as.factor(block_no) + 
                                                  (1|id) + (0 + scale(trial_no)|id), 
                                                closeness_repl, prior = pr, cores = 4, save_pars = save_pars(all = TRUE), iter = 9000)

lmer_bayes_closeness_interaction_repl <- brms::brm(scale(closeness) ~ scale(trial_no)*condition*as.factor(block_no) + 
                                                     (1|id) + (0 + scale(trial_no)|id), 
                                                   closeness_repl, prior = pr, cores = 4, save_pars = save_pars(all = TRUE), iter = 9000)

# interaction effect with trial_no?
bayesfactor_closeness_repl <- bayes_factor(lmer_bayes_closeness_interaction_repl,lmer_bayes_closeness_additive_repl)


############################################################################################
#### analyses with neural betas ############################################################

# neural sensitivity to other's pain is differentially predictive of social closeness
# depending on trial type (observed pain vs. observed non-pain)

# IFG/AI
lmer_IFG_closeness_trait_trial_type <- lmer(closeness ~ betas_PM_wtreatMINcon_IFG*block_no*trial_type + facet*score + 
                                              (1|id) + (1|trial_no), behavior_and_neural_betas)
Anova(lmer_IFG_closeness_trait_trial_type, type = 3)
summary(lmer_IFG_closeness_trait_trial_type)

# TPJ/STS
lmer_TPJ_closeness_trait_trial_type <- lmer(closeness ~ betas_PM_wtreatMINcon_TPJ*block_no*trial_type + facet*score+ 
                                              (1|id) + (1|trial_no), behavior_and_neural_betas)
Anova(lmer_TPJ_closeness_trait_trial_type, type = 3)
summary(lmer_TPJ_closeness_trait_trial_type)


