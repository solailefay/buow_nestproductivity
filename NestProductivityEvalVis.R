#######################################################################
#######################################################################
#######################################################################
#########       This script was created by Solai Le Fay         #######
#########       to evaluate and visualize model results         #######
#########      for Burrowing Owl nest productivity model        #######
######### 
##                                                                    
# We import results for Bayesian analysis of Burrowing Owl nest prod
# which were modeled using a zero-inflated model in a Bayesian framework 
# We visualize model results and assess model fit in this script                 
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse ) 
library( patchwork ) 
library ( jagsUI )
library( ggdist )

################################################################################
#### Load workspace with model results -----------------------------------------
load( "ProdMod_ZIP.RData" )

#rename model
zm <- zip.fullab

###############################################################################
####################     model evaluation  ####################################
###############################################################################

#####################    trace plots    #############################
par( mfrow = c( 2, 2 ), ask = F, mar = c(3,4,2,2) )

#intercept 
traceplot( zm, parameters = c( 'int.lam') )
#random intercepts
traceplot( zm, parameters = c( "sigma.j", "sigma.s") )
#fixed effects
traceplot( zm, parameters = c( 'beta') )



###################    MAE & Bayesian R^2   ##########################
# MAE with 95% CI by calculating per posterior draw [5000 values]
mae_draws <- apply(zm$sims.list$y_hat, 1, function(draw) {
  mean(abs(y_obs - draw))
})

# Summarize
data.frame(
  Estimate = median(mae_draws),
  Q2.5     = quantile(mae_draws, 0.025),
  Q97.5    = quantile(mae_draws, 0.975)
)

# Calculate bayesian r2 (Gelman et al. 2019)
# Observed y
y <- moddf[, yid]

# Extract posterior predicted matrix [5000 draws x n_obs]
mu_mat <- zm$sims.list$y_hat

# Compute Bayes R2 per draw (Gelman et al. 2019)
# = variance of predicted values over the variance of predicted values plus variance of residuals
# mu = each row of mu_mat as apply loops through them. 
# For each of the 5000 draws, mu is a vector of length n_obs (581)
r2_draws <- apply(mu_mat, 1, function(mu) {
  var_fit <- var(mu)
  var_res <- var(y - mu)
  var_fit / (var_fit + var_res)
})

# Summarize
data.frame(
  Estimate = median(r2_draws),
  Q2.5     = quantile(r2_draws, 0.025),
  Q97.5    = quantile(r2_draws, 0.975)
)


### Pull and visualize observed and predicted values
y_obs  <- win.data$y_obs
y_hatz  <- zm$sims.list$y_hat

# calculate posterior mean prediction per observation
y_pred_meanz <- colMeans(y_hatz)

# combine into df
evalplotdfz <- data.frame(
  observed = y_obs,
  predicted = y_pred_meanz,
  year = moddf$year,
  site = moddf$site
)

# Plot observed vs mean predicted 
zeval <- ggplot(evalplotdfz, aes(x = observed, y = predicted, color = factor(year))) +
  geom_jitter(width = 0.2, height = 0.05, size = 2, alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(breaks = 0:max(evalplotdfz$observed)) +
  scale_y_continuous(
    limits = c(0, 12),
    breaks = 0:12
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Observed chick count",
    y = "Predicted mean chick count",
    color = "Year",
    title = "ZIP Full Productivity Model: Observed vs Predicted"
  )
zeval

######### end model evaluation ######################################

###############################################################################
#################### viewing model results ####################################
###############################################################################
# predictor names/order
colnames(X)

#beta summaries
beta_summary <- zm$summary[grep("beta", rownames(zm$summary)), c("mean","2.5%","97.5%")]
#add predictor names
beta_summary <- data.frame(
  predictor = colnames(X),
  beta_summary
)
beta_summary

##### random intercepts 
# year random effect
year_re <- zm$sims.list$eps.j
year_re_df <- data.frame(
  level = 1:ncol(year_re),
  mean = apply(year_re, 2, mean),
  lower = apply(year_re, 2, quantile, 0.025),
  upper = apply(year_re, 2, quantile, 0.975)
)
#variance
zm$mean$sigma.j^2

# site random effect
site_re <- zm$sims.list$eps.s
site_re_df <- data.frame(
  level = 1:ncol(site_re),
  mean = apply(site_re, 2, mean),
  lower = apply(site_re, 2, quantile, 0.025),
  upper = apply(site_re, 2, quantile, 0.975)
)
#variance
zm$mean$sigma.s^2


### PROBABILITY OF DIRECTION
# Probability of direction is the proportiton of the posterior distribution
# that is of the median's sign. I.E. what percent of the posterior is positve
# or negative. It is directly interpretable as the probability that the effect
# is directional. Ref: https://joss.theoj.org/papers/10.21105/joss.01541
# Function to calculate pd
prob_direction <- function( samples ) {
  if( median( samples ) > 0 ) {
    return( mean( samples > 0 ))
  } else {
    return( mean( samples < 0 ))
  }
}

##Coefficients
pd_beta <- apply( zm$sims.list$beta, 2, prob_direction )
pd_beta

########################   Posterior plot   ####################################
# Manuscript Figure 3
#To plot our model coefficient distributions 
#start with extracting relevant fixed effects from model
beta.matz <- zm$sims.list$beta

#rename predictors
colnames(beta.matz) <- c("Neighbor density", "Neighbor density²", 
                         "Distance to\nagriculture (m)","Shrub cover (%)",
                         "Perennial biomass\n(kg/ha)", "Annual biomass\n(kg/ha)", 
                         "Male age", "Female age", "Temperature\n(days)")

#set colors manually
colors <- c(
  "Neighbor density"              = "#7B9E87",
  "Neighbor density²"             = "#7B9E87",
  "Distance to\nagriculture (m)"  = "#C4A882",
  "Shrub cover (%)"               = "#D5B9B2",
  "Perennial biomass\n(kg/ha)"    = "#7A8C99",
  "Annual biomass\n(kg/ha)"       = "#B5836A",
  "Male age"                      = "#8DA9B5",
  "Female age"                    = "#A07A8A",
  "Temperature\n(days)"           = "#A8B87C"
)

#convert beta matrix to df
beta.mat.df <- as.data.frame(beta.matz)
# Convert to long 
beta_long <- pivot_longer(beta.mat.df, cols=everything(), names_to = "parameter", values_to = "value")

# Set order from colors vector
beta_long$parameter <- factor(beta_long$parameter, levels = rev(names(colors)))

## Compute min/max of each posterior for bottom line of posterior in plot
baseline_df <- beta_long %>%
  group_by(parameter) %>%
  summarise(lo = min(value), hi = max(value))

# Figure 3:
#plot
posteriordistplot <- ggplot( beta_long, aes( x =parameter, y = value, fill = parameter )) +
  stat_halfeye(
    trim = FALSE,
    .width = 1,        # shades 95% CI
    point_interval = mean_qi,
    slab_color = "black", # outline color
    slab_size = 0.2, #thin outline
    interval_color = NA,   # removes interval line
    point_color = NA,      # removes point estimate dot
    show.legend = FALSE
  ) +
  #add bottom lines
  geom_segment(
    data = baseline_df,
    aes(x = parameter, xend = parameter, y = lo, yend = hi),
    color = "black",
    linewidth = 0.2,      
    inherit.aes = FALSE
  )+
  scale_fill_manual(values = colors) +
  theme_classic( base_size = 9 ) +
  coord_flip() +
  geom_hline( yintercept = 0, linewidth = 0.5, color = 'black', linetype = "dashed" ) +
  labs(
    y = "Standardized effect size",
    x = ""
  ) +
  theme(
    plot.title = element_text( hjust = -0.06, face = "bold" ),
    axis.ticks = element_blank(),
    axis.text.y  = element_text(size = 8),   # predictor names on y axis
    axis.text.x  = element_text(size = 8),   # numbers on x axis
    axis.title.x = element_text(size = 8)    # "Standardized Effect Size" label
  ) +
  theme( #cut down on empty space 
    plot.margin = margin(
      t = 2,
      r = 2,
      b = 2,
      l = 0,
      unit = "mm"
    )
  )+ #add probability of direction
  annotate( "text", x = 9.13, y = -1.0, label = "98.72%", color = "black", size = 2.5 ) + # density
  annotate( "text", x = 8.13, y = 0.78, label = "90.24%", color = "black", size = 2.5) + # density quad
  annotate( "text", x = 7.13, y = -0.4, label = "64.04%", color = "black", size = 2.5 ) + # dist
  annotate( "text", x = 6.13, y = -0.45, label = "64.38%", color = "black", size = 2.5) + # shrub
  annotate( "text", x = 5.13, y = -1.1, label = "100%", color = "black", size = 2.5 ) + # peren
  annotate( "text", x = 4.13, y = -1, label = "100%", color = "black", size = 2.5 ) + # annual
  annotate( "text", x = 3.13, y = 0.38, label = "99.60%", color = "black", size = 2.5) + # male
  annotate( "text", x = 2.13, y = 0.35, label = "95.84%", color = "black", size = 2.5 ) + # female
  annotate( "text", x = 1.13, y = 0.5, label = "52.70%", color = "black", size = 2.5)  # temp

posteriordistplot

###############################################################################
###################      PARTIAL PREDICTIONS PLOTS      #######################
###############################################################################
# Manuscript Figure 4
# Estimate partial prediction plots (marginal effect plots) for predictors 
# with 95% CIs not overlapping zero: density, annual, perennial, male, female
# predictors with 85% CI not overlapping zero: density2

# Start by creating our datasets to predict over
# how many values do we use:
n <- 100
#define a vector of ones for intercept
int <- rep( 1, n )

#check order of predictors
colnames(X)

#convert biomass predictors from lbs/acre to metric (kg/ha)
moddf$perennial_bio_metric <- moddf$perennial * 1.12085
moddf$annual_bio_metric <- moddf$annual * 1.12085

##### PERENNIAL ---------------------------------------
# Use the observed values to define range of predictor:
peren <- seq( min( moddf[,"perennial_bio_metric"]),max( moddf[,"perennial_bio_metric"]),
                  length.out = n )
#standardize predictors:
peren.std <- scale2sd( peren)

#extract relevant fixed coefficient from model results
fixedperen <- cbind( zm$sims.list$int.lam, zm$sims.list$beta[,5] )

#estimate predicted values
predperen <- exp( fixedperen %*% t( cbind( int, peren.std) ) )
#calculate mean values
mperen <- apply( predperen, MARGIN = 2, FUN = mean )
#calculate 95% credible intervals
CIperen <- apply( predperen, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
perendf <- data.frame( mperen, t(CIperen),
                       peren.std, peren )
#view
head( perendf); dim( perendf)
#rename columns
colnames(perendf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalized effects 
p_peren <- ggplot(perendf, aes(x = peren, y = Mean)) +
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha = 0.3, fill = "#7A8C99") +
  geom_line(size = 0.4, color = "#7A8C99") +
  scale_y_continuous(
    limits = c(0,10),
    breaks = seq(0,10,2)
  ) +
  theme_classic(base_size = 9) +
  ylab("Nest productivity") +
  xlab("Perennial biomass (kg/ha)") +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
  )+
  annotate( "text", x = 695, y = 7.5, label = "A", color = "black", size = 4) 

p_peren

##### ANNUAL ---------------------------------------
# Use the observed values to define range of predictor:
annual <- seq( min( moddf[,"annual_bio_metric"]),max( moddf[,"annual_bio_metric"]),
               length.out = n )
#standardize predictors:
annual.std <- scale2sd( annual )

#extract relevant fixed coefficient from model results
fixedannual <- cbind( zm$sims.list$int.lam, zm$sims.list$beta[,6] )

#estimate predicted productivity
predannual <- exp( fixedannual %*% t( cbind( int, annual.std) ) )
#calculate mean productivity
mannual <- apply( predannual, MARGIN = 2, FUN = mean )
#calculate 95% credible intervals
CIannual <- apply( predannual, MARGIN = 2, FUN = quantile, 
                   probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
annualdf <- data.frame( mannual, t(CIannual),
                        annual.std, annual )
#view
head( annualdf); dim( annualdf)
#rename columns
colnames(annualdf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalized effects 
p_annual <- ggplot(annualdf, aes(x = annual, y = Mean)) +
  # ylab("Nest productivity") +
  xlab("Annual biomass (kg/ha)") +
  ylab(NULL) +
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha = 0.3, fill = "#B5836A") +
  geom_line(size = 0.4, color = "#B5836A") +
  scale_x_continuous(breaks = seq(400, 1400, by = 300)) +
  scale_y_continuous(
    limits = c(0,10),
    breaks = seq(0,10,2)
  )+
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  annotate( "text", x = 1233, y = 7.5, label = "B", color = "black", size = 4)

p_annual

##### DENSITY (QUADRATIC) ---------------------------------------
# Conspecific neighbor density
# Use the observed values to define range of predictor:
density <- seq( min( moddf[,"density"]),max( moddf[,"density"]),
               length.out = n )
#standardize predictors:
density.std <- scale2sd( density )
density2.std <- scale2sd( density^2 )

#extract relevant fixed coefficient from model results
fixeddensity <- cbind( zm$sims.list$int.lam, zm$sims.list$beta[,1] ,zm$sims.list$beta[,2])

#estimate predicted productivity
preddensity <- exp( fixeddensity %*% t( cbind( int, density.std, density2.std) ) )
#calculate mean productivity
mdensity <- apply( preddensity, MARGIN = 2, FUN = mean )
#calculate 95% credible intervals 
CIdensity <- apply( preddensity, MARGIN = 2, FUN = quantile, 
                   probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
densitydf <- data.frame( mdensity, t(CIdensity),
                        density.std, density )
#view
head( densitydf); dim( densitydf)
#rename columns
colnames(densitydf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalized effects 
p_density <- ggplot(densitydf, aes(x = density, y = Mean)) +
  ylab(NULL) +
  xlab("Neighbor density") +
  geom_ribbon(aes(ymin = lowCI, ymax = highCI), alpha = 0.3, fill = "#7B9E87") +
  geom_line(size = 0.4, color = "#7B9E87") +
  scale_y_continuous(
    limits = c(0,10),
    breaks = seq(0,10,2)
  )+
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  annotate("text", x = 58, y = 7.5, label = "C", color = "black", size = 4)

p_density

#### Figure 4
#### combine into panel
margeff_panel <-  p_peren + p_annual + p_density


##### AGE (MALE & FEMALE) ---------------------------------------
##### Figure 5
male <- c(-1, -0.25, 0, 0.25, 1)
female <- c(-1, -0.25, 0, 0.25, 1)
ageint <- rep( 1, 5 )

#extract relevant fixed coefficient from model results
fixedmale <- cbind(zm$sims.list$int.lam, zm$sims.list$beta[,7]) 
fixedfemale <- cbind(zm$sims.list$int.lam, zm$sims.list$beta[,8]) 

#estimate predicted productivity
predmale <- exp( fixedmale %*% t( cbind( ageint, male) ) )
predfemale <- exp( fixedfemale %*% t( cbind( ageint, female) ) )

#calculate mean productivity
mmale <- apply( predmale, MARGIN = 2, FUN = mean )
mfemale <- apply( predfemale, MARGIN = 2, FUN = mean )
#calculate 95% credible intervals
CImale<- apply( predmale, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )
CIfemale<- apply( predfemale, MARGIN = 2, FUN = quantile, 
                probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
maledf <- data.frame( mmale, t(CImale),
                        male = c("SY", "AHY", "TY", "ASY", "ATY" ))
femaledf <- data.frame( mfemale, t(CIfemale),
                      female = c("SY", "AHY", "TY", "ASY", "ATY" ))
#view
head( maledf); dim( maledf)
head( femaledf); dim( femaledf)
#rename columns
colnames(maledf )[1:3] <- c(  "Mean", "lowCI", "highCI" )
colnames(femaledf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#add sex column
maledf$sex <- "Male"
femaledf$sex <- "Female"

#change age column names to "age" instead of "male" and "female"
colnames(maledf)[colnames(maledf) == "male"] <- "Age"
colnames(femaledf)[colnames(femaledf) == "female"] <- "Age"
#merge
agedf <- rbind(maledf, femaledf)

#change age classes to ordinal factor
agedf$Age <- factor(agedf$Age,
                    levels = c("SY", "AHY", "TY", "ASY", "ATY"))

# Figure 5
# Plot both sexes together
margeff_sex <- ggplot(agedf, aes(x = Age, y = Mean, color = sex)) +
  # Points
  geom_point(position = position_dodge(width = 0.4),
             size = 1.4) +
  # Error bars
  geom_errorbar(aes(ymin = lowCI, ymax = highCI),
                position = position_dodge(width = 0.4),
                width = 0.25,
                size = 0.6) +
  # Labels
  labs(
    x = "Age class",
    y = "Nest productivity",
    color = NULL
  ) +
  # Colors
  scale_color_manual(
    values = c(
      "Male" = "#8DA9B5", 
      "Female" = "#A07A8A"  
    )
  ) +
  scale_y_continuous(
    limits = c(3.5, 6),
    breaks = seq(4, 6, 1)
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.6, "cm"),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

margeff_sex

##########      end of marginal effect plots      ##############################
################################################################################
##########      end of script                     ##############################