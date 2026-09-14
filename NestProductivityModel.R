#######################################################################
#######################################################################
#######################################################################
#########       This script was created by Solai Le Fay         #######
#########         to fit a productivity model for BUOW          #######
######### 
##
## This script was adapted by Solai from script written by Dr. Jen Cruz 
##                                                                    
## Here we import our cleaned data, containing each nest attempt 
## with a count of chicks survived to near-fledging, and all ecological
## predictors: annual and perennial biomass, shrub % cover,
## distance to agriculture, hot days, male and female age
## with random effects of year and siteID
##
## The model is hierarchical    
## with : (1) an ecological submodel linking productivity to             
## environmental predictors for each nest attempt; 
## (2) a zero-inflated submodel           
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse ) 
library( jagsUI )
###################################################################
#### Load data ----------------------------------------------------

#load relevant workspace---------------
load( "ProdMod_ZIP.RData" )

# set directory where your data are:
datapathclean <- "Z:Common/BurrowingOwls/CleanData/"


# load unscaled model predictor df 
moddf <- read.csv( file = paste( datapathclean, "prodmoddf_ab.csv", sep = ""),
                   header = TRUE )
#581 nest attempts

#view
head( moddf ); dim( moddf ) 


#### End of data load -------------
####################################################################
##### Ready data for analysis --------------

# extract broad parameters of interest
#number of nest attempts:
I <- length( moddf$nestID )
# No detection model/repeated observations to loop through

# create index for chick count observations - response
yid <- grep( "count", colnames( moddf ), value = FALSE)

### Prepare fixed effects
# create dataframe for fixed predictors and scale numeric by 2 SD 
# (for comparisons with categorical predictors)
# Scaling function for 2 sd
scale2sd <- function(x){
  (x - mean(x))/(sd(x)*2)
}

## add quadratic term for delta before scaling
moddf$delta2 <- (moddf$delta)^2

# scale and rename predictors
XI <- moddf %>%
  dplyr::select(delta, delta2, edge = dist_to_edge_meters, shrub, peren = perennial_bio,
                annual = annual_bio, male = age_male_5C, female = age_female_5C, temp = tnz_high) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), scale2sd))


#view
head( XI); dim( XI ) 

# JAGS cannot take predictors as factors
# Space factor levels across numeric range from -1 to 1
X <- XI %>%
  mutate(
    male   = as.character(male),
    female = as.character(female)
  ) %>%
  mutate(
    male = recode(male,
                  "SY"  = -1,
                  "AHY" = -0.5,
                  "TY"  = 0,
                  "ASY" = 0.5,
                  "ATY" = 1),
    female = recode(female,
                    "SY"  = -1,
                    "AHY" = -0.5,
                    "TY"  = 0,
                    "ASY" = 0.5,
                    "ATY" = 1)
  )


head(X)


# RANDOM INTERCEPTS
# We will include year and siteID as random effects - turn both numeric
# Year random intercept
moddf$year <- as.numeric( as.factor(moddf$year ))

#number of years
J <- length(unique(moddf$year)) 

# Site ID random intercept
moddf_sub$site <- as.numeric( as.factor(moddf_sub$siteID ))

#number of sites
S <- length(unique(moddf_sub$site)) #93


#-------------------------------------------------------------------------
##########################################################################
####### Zero inflated poisson - nest productivity model
### Year and siteID are random intercepts
### Fixed effects include weather, vegetation, and social factors:
#- Weather: tnz_high (number of days that season where max temp was above the max thermoneutral zone)
#- Vegetation:
# perennial_bio, annual_bio, shrub (perennial and annual mean biomass within 1400m radius of burrow, % shrub cover)
# dist_to_edge_meters (nearest distance of nest to agriculture in meters)
#- Social: neighbor density (delta), male age, female age
# --- Delta as quadratic, male and female age as categorical factors with 5 levels each
############################################################################
############## Specify model in bugs language:  #####################
sink( "zip.fullab.txt" )
cat( "
     model{
     
      ### PRIORS
      #random intercept for year
      for( j in 1:J ){
        eps.j[j] ~ dnorm( 0, pres.j ) T(-7, 7)
      }
      #associated variance of random intercept:     
      pres.j <- 1/ ( sigma.j * sigma.j )
      #sigma prior specified as a student t half-normal:
      #SD can not have negative values, trunkate to make it positive
      sigma.j ~ dt( 0, 2.5, 7 ) T( 0, )
      
      #random intercept for site
      for( s in 1:S ){
        eps.s[s] ~ dnorm( 0, pres.s ) T(-7, 7)
      }
      #associated variance of random intercepts:     
      pres.s <- 1/ ( sigma.s * sigma.s )
      #sigma prior specified as a student t half-normal:
      #SD can not have negative values, trunkate to make it positive
      sigma.s ~ dt( 0, 2.5, 7 ) T( 0, )
  
      #priors for fixed coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[ b ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      # Zero-inflated prior
      omega ~ dbeta( 4, 4 )
      
      #prior for ecological model intercept 
      int.lam ~ dnorm( 0, 0.01 )
      
    
    # ecological model - nest productivity
    # loop through each nest attempt i
    for( i in 1:I ){
      #latent suitability state 
      #omega is not related to something else (only on right side of equations)
      #so it was given a prior above
      # Z = zero inflation, whether a nest attempt was successful (at least one chick) or failed (0 chicks)
      # distributed as a bernoulli conditional on omega
      z[ i ] ~ dbern( omega )
      
      # true productivity now conditional on z (nest success)
      # z is 0 or 1, so either multiplying by 0 or just lamda 
      # observed count of chicks for each nest attempt is distributed as a poisson
      # conditional on lambda and z (for zero inflation, nest success)
      y_obs[ i ] ~ dpois( lambda[ i ] * z[ i ] )
      
      #mean relative productivity related to ecological predictors
      log( lambda[ i ] ) <- int.lam + 
                        inprod( beta, X[ i, ]  ) +
                        #random intercepts
                        eps.j[year[i]] 
                        + eps.s[site[i]]
                        
            
        # Model evaluation
        #for model evaluation we estimate expected productivity yhat
        y_hat[ i ] ~ dpois( lambda[ i ] * z[ i ] )
        
    } #close I

    } #model close
     
     ", fill = TRUE )

sink()


################ end of model specification  #####################################     
modelname <- "zip.fullab.txt"
#parameters monitored
params <- c(
  'int.lam' #intercept for lamda
  #, 'sigma.j' #random intercept for year
  #, 'sigma.s' #random intercept for site
  , 'sigma.j' #error for random intercept for year
  , 'sigma.s' #error for random intercept for site
  , 'omega' #suitability parameter
  , 'beta' #abundance coefficients
  , 'y_hat' #predicted observations
)

#how many ecological predictors that are fixed effects
B <- dim(X)[2]

#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           z = rep(1, I))}


#define data that will go in the model
str( win.data <- list( y_obs = moddf[ ,yid] ,
                       #number of nest attempts, years, sites, and fixed predictors
                       I = I, J = J, B = B, S = S,
                       #ecological predictors
                       X = X,
                       #random effects
                       year=moddf$year, site=moddf$site
                       
) )


#call JAGS and summarize posteriors:
zip.fullab <- autojags( win.data, inits = inits, params, modelname, #
                        n.chains = 5, n.thin = 10, n.burnin = 20000,
                        iter.increment = 10000, max.iter = 500000, 
                        Rhat.limit = 1.1,
                        save.all.iter = FALSE, parallel = TRUE ) 


###### end zip.fullab model***** ########----------------------------------------