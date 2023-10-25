# R-script for producing the figures for the article
# - Cycling fall threshold based on perturbed cyclist experiments -- version 2"
# "Marco Reijne, Frank van der Meulen, Arend Schwaab and Frans van der Helm"

library(tidyverse)
library(DataExplorer)
library(brms)
library(rstanarm)
library(ggdist)  
library(caret)
library(rpart.plot)
library(gridExtra)
library(lme4)
library(BAS)
library(latex2exp)

# Set theme for ggplot
mytheme = theme_bw()
theme_set(mytheme)  

setwd("~/Sync/DOCUMENTS/onderzoek/fietsvalproject/analyse/supplementary")

# Read the data 
d <-  read_csv("PerturbedCyclistExperimentalDataTUDelft2022.csv", 
        col_types = cols(participant_ID = col_character(), 
        outcome = col_factor()))

# We reorder `participant_ID` according to `age`, filter out initial search data and remove some columns that we don't further need in the analysis.
d <- d %>% 
  mutate(participant_ID = fct_reorder(participant_ID, age)) %>%
  filter(initial_search==FALSE, exclude=="NO") %>% 
  select(-c(pull_force, comments, initial_search, exclude)) 

# Exploratory plot
p_exp <- d %>% mutate(velocity=as.factor(velocity)) %>%
  ggplot(aes(x=angular_momentum, y = outcome, colour=velocity)) + 
  geom_jitter(height=0.1,size=0.5) + 
  facet_wrap(~participant_ID,nrow=4) + 
  labs(x="angular momentum", y="y")+
  theme(legend.position = "top")

pdf("exploratory.pdf",width=7, height=4)
p_exp
dev.off()

# Random forest
d_dropna  <- d %>% drop_na()
X <- d_dropna %>% dplyr::select(-outcome)
Y <- d_dropna %>% dplyr::select(outcome) %>% deframe()
set.seed(14)
tr_rf2 <- train(X,Y, method='rf',
            trControl=trainControl("cv"), importance = TRUE, tuneLength=10)
p2 <- ggplot(varImp(tr_rf2))

# Refit, but now further excluding `participant_ID`.
Xsub <- X %>% dplyr::select(-participant_ID)
set.seed(14)
tr_rf3 <- train(Xsub, Y, method='rf',
             trControl=trainControl("cv"), importance = TRUE, tuneLength=10)
 p3 <- ggplot(varImp(tr_rf3))

# Put the two figures side by side
pdf('fig_randomforest.pdf', width=7, height=3.5)
  grid.arrange(p2,p3, ncol=2)
dev.off()

# Bayesian Model Averaging

# First, we standardise (mean centered, scaled to have unit standard deviation) the numerical variables. This helps to interpret th magnitude of the estimated coefficients.
preProcValues <- preProcess(d, method = c("center", "scale"))
dTransformed <- predict(preProcValues, d) %>% drop_na()

# Bar-chart of number of experiments for each participant
pdf("n_experiments_participants.pdf",width=7, height=3)
  dTransformed %>% group_by(participant_ID) %>%  count() %>%
    ggplot(aes(x=participant_ID,y=n)) + geom_bar(stat='identity')
dev.off()

# We use Bayesian model averaging using the `BAS`-library to gain insight into which variables are important to be included. We force categorical variables to be fully included or not (hence we don't allow certain factor levels to be included, while other levels of the same factors not). We insist on a model that contains `angular_momentum`, for else we cannot determine the required threshold value from the inferred model. 
# We use a uniform prior over all possible models and use the default mixture of Zellner-g-prior on the coefficients. 

IT <- 20000 # nr of MCMC iterations
set.seed(15)
bma0 = bas.glm(outcome ~ ., data=dTransformed, method="MCMC", 
               MCMC.iterations=IT, family=binomial(), 
               include.always = ~ angular_momentum,
               modelprior=uniform(), force.heredity=TRUE)
print(summary(bma0),digits=2)
diagnostics(bma0) # should be on straight line (assesses MCMC sampling)

pheat <- image(bma0, rotate=F)
pdf('heat_bma.pdf', width=7, height=10)
  image(bma0, rotate=F)
dev.off()

# make figure of posterior inclusion probabilities
aa <- tibble(probne0=bma0$probne0, names=bma0$namesx) 
partic_p <- aa %>% filter(names=="participant_ID1") %>% select(probne0)

p_inclusion <- aa %>% 
  filter(!str_detect(names,"participant_ID")) %>%
  filter(!str_detect(names,"Interc")) %>%
  add_row(names="participant_ID", probne0=as.numeric(partic_p)) %>% 
  mutate(names2= reorder(names, probne0, median))  %>% 
  ggplot(aes(x=probne0, y=names2)) + 
  geom_point() + geom_vline(xintercept=0.5, colour="red") +
  labs(x="posterior inclusion probability", y="")

pdf('prob_inclusion_bma.pdf', width=3.5, height=3)
  p_inclusion
dev.off()


# Multilevel approach

# We first set the prior distributions on the coefficients to be zero-mean Normally distributed with standard-deviation of 2.

prior1 <- set_prior("normal(0, 2)", class = "b")

set.seed(16)
# We fit a multilevel version based on the median probability model appearing from `bma0`. 
fit0 <- brm(outcome ~ 0 + Intercept + phi_zero + pull_direction +
            skill_performance + delta_dot_zero + velocity +
            pull_ID + angular_momentum +  (1 | participant_ID), 
            data=dTransformed,
            family = bernoulli(link="logit"),
            prior=prior1,
            warmup = 5000, 
            iter = IT, 
            chains = 4, 
            seed = 123)

fit0 %>% summary()

prior_summary(fit0)

# The PSIS-LOO score has been advocated for Bayesian model selection:
l0 <- loo(fit0)

# We also consider the model under `fit0`, but dropping delta_dot_zero and phi_zero
fit2 <- brm(outcome ~ 0 + Intercept + pull_direction +
              skill_performance + velocity +
              pull_ID + angular_momentum +  (1 | participant_ID), 
            data=dTransformed,
            family = bernoulli(link="logit"),
            prior=prior1,
            warmup = 5000, 
            iter = IT, 
            chains = 4, 
            seed = 123)
l2 <- loo(fit2)


loo_compare(l0, l2) # stick to fit0


# Compare Bayesian with frequentist results for model 0
dTransformed_reformat <- dTransformed %>%
                  mutate(outcomenum = as.numeric(outcome)-1)
glmer_fit0 <- glmer(outcomenum ~  phi_zero + pull_direction +
                      skill_performance + delta_dot_zero + velocity +
                      pull_ID + angular_momentum +  (1 | participant_ID), 
                     data=dTransformed_reformat, family = binomial)

# To facilitate comparison, we show both the Bayesian and frequentist fitted models:
sum0 <- fit0 %>% summary() %>% unclass() 
sum0$fixed
sum0$random
# to be compared to
glmer0 <- glmer_fit0 %>% summary() %>% unclass() 
coefficients(glmer0)
glmer0$varcor

# The posterior means from `brm` and estimates from `glmer` are very similar. 
# Additionally, we obtain significance of all variables in the model. 





### Posterior predictive check

checks2 <- function(d, fit) # fraction falling by participant and factor(velocity)
{
  pred <- posterior_predict(fit, ndraws=1000)
  fracfalling1 <- d %>% 
    mutate(velocity=as.factor(round(velocity,1))) %>%     
    mutate(predicted_average=colMeans(pred)) %>% 
    group_by(velocity, participant_ID) %>%
    summarise(n=n(), observed=mean(outcome=="1"),
              predicted=mean(predicted_average)) %>%
    pivot_longer(cols=c(predicted, observed), names_to="type", values_to="fraction")
  
  p <- fracfalling1 %>% 
    ggplot() + 
    geom_point(aes(x=velocity, y = fraction,colour=type)) +
    labs(x="velocity (standardized)", y="fraction of falls of the total number of perturbations") 
  p + facet_wrap(~participant_ID) 
}

pdf("predictive_check.pdf", width=7, height=7)
  checks2(dTransformed, fit0)
dev.off()


# We extract `IT`  posterior samples from  the fitted model. 
ps <- brms::as_draws_df(fit0) %>% tibble() %>% sample_n(IT)

# We plot all parameters, except for participant effects

ln <- c("b_Intercept", "b_skill_performance", "b_velocity", "b_phi_zero",
        "b_pull_ID", "b_pull_directionCW",  "b_angular_momentum",
         "b_delta_dot_zero","sd_participant_ID__Intercept")

ps_coef <- ps %>% dplyr::select(-contains("r_particip")) %>%
  dplyr::select(contains("b_") | contains("sd_")) %>%
  tidyr::pivot_longer(cols=ln,names_to="coefficient", values_to="value") %>%
  mutate(coefficient=str_remove(coefficient,"b_")) %>% 
  mutate(coefficient = recode(coefficient, `sd_participant_ID__Intercept` = "sd_participant"))

pdf("posterior_coefs.pdf",width=7, height=4.5)
  ps_coef %>% ggplot(aes(x=value, y = reorder(coefficient, value,median))) +
    stat_halfeye( point_size=1, interval_size=0.5)+
    labs(x="coefficient", y = "")
dev.off()

# Similarly for  participant effects
ps_sub <- ps %>% dplyr::select(contains("r_particip"))
ps_sublong <-  ps_sub %>% tidyr::pivot_longer(cols=contains("r_particip"),
                names_to="participant", values_to="value") %>%
                mutate(participant_ID=as.factor(readr::parse_number(participant)))
ps_complete <- left_join(ps_sublong, d, by ="participant_ID")

pdf("posterior_participants.pdf",width=7, height=7)
  ps_complete %>% ggplot(aes(x=value, 
            y = reorder(participant_ID, value,median),fill=age )) +
            stat_halfeye(point_colour="orange", point_size=1, interval_size=0.5)+
         #   labs(x="coefficient", y = "participant_ID")
    labs(x=TeX(r'($\beta_i$)'), y="participant_ID")
dev.off()



rm(ps_coef, ps_sub, ps_sublong, ps_complete) # free memory


## Computing perturbation threshold for a few participants in the study.

# 
# Consider the results for following persons:
#   
# - proefpersoon 1 bij 12 km/h en 6 km/h
# - proefpersoon 9 bij 12 km/h en 6 km/h
# - proefpersoon 11 bij 12 km/h, bij 6 km/h en bij 18 km/h
# - eventueel proefpersoon 12 bij 12 km/h en bij 6 km/h
# 
#  phi_zero = 0 (dit betekent dat de fiets rechtovereind staat en niet al naar een kant leunt)
#  pull_directionCW = 1
#  delta_dot_zero = 0
#  pull_id = 1


pIDs <- c(1,1,9,9,11,11,11,12,12) # arbitrarily choose a few participants

d1 <- d %>% filter(participant_ID==as.character(1)) %>% slice(1)
d9 <- d %>% filter(participant_ID==as.character(9)) %>% slice(1)
d11 <- d %>% filter(participant_ID==as.character(11)) %>% slice(1)
d12 <- d %>% filter(participant_ID==as.character(12)) %>% slice(1)
d_existing <- bind_rows(d1,d1,d9,d9,d11,d11,d11,d12,d12) %>% 
  mutate(participant_ID2=c("P1_v6","P1_v12","P9_v6", "P9_v12", "P11_v6","P11_v12","P11_v18", "P12_v6", "P12_v12" )) %>%
  mutate(velocity=c(6,12,6,12,6,12,18,6,12)) %>%
  mutate(phi_zero=rep(0,9), pull_direction=rep("CW",9), delta_dot_zero=rep(0,9), pull_ID=rep(1,9))

d_existingTransformed <- predict(preProcValues, d_existing)


nps <- nrow(ps) # nr of posterior samples
zz = rep(0,nps)
z <- data.frame("P1_v6" = zz,"P1_v12" = zz,"P9_v6" = zz, "P9_v12" = zz, "P11_v6" = zz,
                "P11_v12" = zz,"P11_v18" = zz, "P12_v6" = zz, "P12_v12" = zz)
head(z)

# Function to compute threshold function
angular_threshold <- function(participant_ch, postdraw, participant_postdraw)
{
    xx <- participant_postdraw +
        postdraw$b_Intercept + 
        postdraw$b_skill_performance * participant_ch$skill_performance +
        postdraw$b_velocity * participant_ch$velocity+
        postdraw$b_pull_ID * participant_ch$pull_ID +
        postdraw$b_pull_directionCW * (participant_ch$pull_direction=="CW") + 
      #  postdraw$b_phi_zero * participant_ch$phi_zero +
        postdraw$b_delta_dot_zero * participant_ch$delta_dot_zero 
    out <- -xx/(postdraw$b_angular_momentum) 
    out
}  

for (j in 1:9)
{
  ID <- pIDs[j]
  # select for the participant posterior draws from its intercept:
  lab <- paste0("r_participant_ID[",as.character(ID),",Intercept]")
  raneff <- ps[,lab] %>% deframe()
  #  choose participant experimental setting:
  participant_ch <- d_existingTransformed[j,] 

  for (i in 1:nps)  # for each posterior sample
  { 
    postdraw <- ps %>% slice(i) 
    z[i, j] <- angular_threshold(participant_ch, postdraw, raneff[i])
  }
}
m_ang <- mean(na.omit(d)$angular_momentum)
sd_ang <- sd(na.omit(d)$angular_momentum)

angular_threshold_df <- apply(z,2, function(x) {  m_ang + sd_ang * x  }  )  # rescale each column
angular_threshold_tidy <- as_tibble(angular_threshold_df) %>%
                          pivot_longer(cols=1:9, names_to="participant", values_to="eta") 


TeX(r'($\beta_i$)')

# Plotting the results 
pdf("angular_thresholds_instudy.pdf",width=7, height=4.5)
angular_threshold_tidy %>% 
  #ggplot(aes(x=eta, y=reorder(participant, eta, median))) +
  ggplot(aes(x=eta, y=participant)) +
  stat_halfeye(fill="lightblue") +
  labs(x=TeX(r'($\Delta L_{50\%}$)'), y ="") +
  ggtitle("Threshold value for angular momentum")
dev.off()




## Computing perturbation threshold for six new fictive participants.

# We manually construct "new" participants, for which we will compute the perturbation threshold (for the variables in model 5 we provide values, while all remaining columns are just filled with NA-values).


# Nieuwe fictieve personen
# - skill performance 40, velocity 12 km/h
# - skill performance 90, velocity 12 km/h
# - skill performance 40, velocity 6 km/h
# - skill performance 90, velocity 6 km/h 
# - skill performance 40, velocity 18 km/h 
# - skill performance 40, velocity 25 km/h (snelheid buiten bereik experiment)

dpred <- tibble(participant_ID=c("new_v12_sp40", "new_v12_sp90", "new_v6_sp40", "new_v6_sp90","new_v18_sp40","new_v25_sp40"),
                velocity=c(12,12,6,6,18,25), age = rep(NA,6),
                skill_performance=c(40,90,40,90,40,40), pull_ID=rep(1,6),
                pull_direction=rep("CW",6), weight=rep(NA,6), length=rep(NA,6),
                skill_effort = rep(NA,6), initial_search=rep(NA,6), phi_zero=rep(0,6),
                phi_dot_zero=rep(NA,6), delta_zero=rep(NA,6), delta_dot_zero=rep(0,6),
                psi_zero=rep(NA,6), Q_zero=rep(NA,6), comments=rep(NA,6),
                exclude=rep(NA,6), reaction_time=rep(NA,6), pull_force=rep(NA,6),
                angular_momentum=rep(NA,6))

# Now we will predict for these new participants:

dpred_tf <- predict(preProcValues, dpred)  # tf=transformed as we did for the model
dpred_tf$delta_dot_zero <- rep(0,6)

nparticipants = nrow(dpred)  # nr of participants for which we wish to predict
zz = rep(0,nps)
z <- data.frame(new_participant1 = zz,  new_participant2 = zz,
                new_participant3 = zz,new_participant4 = zz,
                new_participant5 = zz, new_participant6 = zz)

for (j in 1:nparticipants)  # for each participant we wish to find threshold value
{
  for (i in 1:nps)  # for each participant
  { 
    particip =  rnorm(1,sd=ps$sd_participant_ID__Intercept[i])
    z[i, j] <- angular_threshold(dpred_tf[j,], ps[i,], particip)
  }
}  
# rescale each column:
angular_threshold_df_new <- apply(z,2, function(x) { m_ang + sd_ang * x})  
angular_threshold_tidy_new <- as_tibble(angular_threshold_df_new) %>%
      pivot_longer(cols=contains("new"), names_to="participant", values_to="eta") %>%
      mutate(participant=fct_recode(participant, 
                                    s40v12="new_participant1",
                                    s90v12="new_participant2",
                                    s40v6="new_participant3",
                                    s90v6="new_participant4",
                                    s40v18="new_participant5",
                                    s40v25="new_participant6"))



# Plotting the results 
pdf("angular_thresholds_fictive.pdf",width=7, height=4.5)
angular_threshold_tidy_new %>% 
#  ggplot(aes(x=eta, y=reorder(participant, eta, median))) +
  ggplot(aes(x=eta, y=participant)) +
  stat_halfeye(fill="lightblue") +
  labs(x=TeX(r'($\Delta L_{50\%}$)'), y ="") +
  ggtitle("Threshold value for angular momentum")
dev.off()

# We can also make a combined plot for participants which were part of the study and "new" participants:

nr <- nrow(angular_threshold_tidy)
nr_new <- nrow(angular_threshold_tidy_new)
factor_new <- c(rep("no",nr), rep("yes",nr_new))
pdf("angular_thresholds.pdf", width=7,height=7)
  bind_rows(angular_threshold_tidy, angular_threshold_tidy_new) %>% 
    mutate(new=factor_new) %>%  
    ggplot(aes(x=eta, y=reorder(participant, eta, median), fill=new)) +
    stat_halfeye(fill="lightblue") + labs(x="eta", y = "participant") + facet_wrap(~new)
dev.off()


