#' ---
#' title: "Integrated SDM - possum data"
#' author: "Saras Windecker & David Uribe"
#' date: "14 May, 2023"
#' output: html_document
#' ---

#' packages
library(raster)
library(jagsUI)

#' covariate for distribution of possums
tree_cover <- readRDS('data/nveg_ras.rds')

# extract center and scale for full region
# this is key! you need to apply the same center
# and scaling to future extractions of the raster.
tree_cover_scaled <- scale(getValues(tree_cover))
scale <- attr(tree_cover_scaled, 'scaled:scale')
center <- attr(tree_cover_scaled, 'scaled:center')

#' COUNT DATA
counts <- read.csv('data/counts.csv')
counts_sp <- SpatialPoints(counts[, c('lon', 'lat')])

# visualise counts
plot(tree_cover)
points(counts_sp, pch = 20)

# area of count data is search radius 50m. convert to sqkm here
area_counts <- (pi * 50 ^ 2) / 1e6

tree_cover_counts <- as.vector(scale(
  extract(tree_cover, 
          counts_sp), 
  center = center,
  scale = scale))

#' DETECTION DATA
det <- read.csv('data/det.csv')
det_sp  <- SpatialPoints(det[, c('lon', 'lat')])

# visualise det data
plot(tree_cover)
points(det_sp, pch = 20)

# search radius 10m to sqkm
area_det <- (pi * 10 ^ 2) / 1e6

tree_cover_det <- as.vector(scale(
  extract(tree_cover, 
          det_sp), 
  center = center,
  scale = scale))

#' PRESENCE-ONLY DATA
po_raw <- read.csv('data/po.csv')

po_ras <- rasterize(po_raw, tree_cover, 
                    fun = 'count', background = 0)
po <- getValues(po_ras)
table(po)

# visualize data
plot(tree_cover)
points(SpatialPoints(po_raw), pch = 20, cex = .4)

# area for po data is sqkm
area_ras <- area(tree_cover)
area_po <- getValues(area_ras)

# covariates for po data
tree_cover_po <- as.vector(tree_cover_scaled)

access_ras <- readRDS('data/access_ras.rds')
access <- as.vector(scale(getValues(access_ras)))

#' predictions across whole region
tree_cover_pred <- as.vector(tree_cover_scaled)
access_pred <- access

#' create jags data
win.data <- list(
  
  # spotlight count data
  n_counts = nrow(counts),
  counts = counts$count,
  tree_cover_counts = tree_cover_counts,
  area_counts = area_counts,
  
  # det data 
  n_det = nrow(det),
  tree_cover_det = tree_cover_det,
  area_det = area_det,
  
  # presence-only data
  n_po = length(po), 
  po = po,
  tree_cover_po = tree_cover_po, 
  area_po = area_po,
  access = access,
  
  # prediction
  n_pred = length(tree_cover_pred),
  tree_cover_pred = tree_cover_pred,
  access_pred = access_pred
)

sink("jags_mod.txt")
cat("

# Define model and write into R working directory

model {

  ## Priors ##
  # nb. these have not been thought about at all! They are generic. 
  
  # Priors for abundance of possums
  alpha ~ dnorm(0, 0.1)

  # Prior for slope of abundance predictor
  beta ~ dnorm(0, 0.1) 

  # Prior for possum spatial sampling bias
  alpha_bias ~ dnorm(0, 0.1)

  # Prior for slope in bias predictors
  beta_bias ~ dnorm(0, 0.1)

  ## Likelihood ##

  # count data of possums
  for (i in 1:n_counts){

    counts[i] ~ dpois(lambda_1[i]*area_counts)
    log(lambda_1[i]) <- alpha + beta*tree_cover_counts[i]

  }

  # detection data
  for (k in 1:n_det) { 

    det[k] ~ dbern(prob_det[k])
    
    # prob_det is prob osberving 1 or more possums given expected abund possums
    prob_det[k] <- 1 - exp(-abund_possum[k]) 
    abund_possum[k] <- lambda_3[k]*area_det
    log(lambda_3[k]) <- alpha + beta*tree_cover_det[k]
    
  }
  
  # presence-only possum data
  for (j in 1:n_po){

    po[j] ~ dpois(lambda_2[j]*bias[j]*area_po[j])
    
    log(lambda_2[j]) <- alpha + beta*tree_cover_po[j] 
    log(bias[j]) <- alpha_bias + beta_bias*access[j]
   
  }

  ## Derived quantities ##
  
  for (i in 1:n_pred){

    log(lambda_pred[i]) <- alpha + beta*tree_cover_pred[i]
    log(bias_pred[i]) <- alpha_bias + beta_bias*access_pred[i]

  }

}

", fill = TRUE)
sink()

#' MCMC
# Run JAGS, check convergence and summarize posteriors:
out <- jagsUI::jags(data = win.data, inits = NULL,
                    parameters.to.save = c("beta", "alpha",
                                           "alpha_bias", "beta_bias",
                                           "lambda_pred"),
                    model.file = "jags_mod.txt",
                    n.chains = 3, n.thin = 1, n.iter = 5000, n.burnin = 2500,
                    codaOnly = TRUE, parallel = FALSE)
# saveRDS(out, "jags_mod.rds")
# out <- readRDS('jags_mod.rds')

#' prediction map
pred_map <- tree_cover
values(pred_map) <- as.vector(out$mean$lambda_pred)
pred_map <- mask(pred_map, tree_cover)
plot(pred_map)

