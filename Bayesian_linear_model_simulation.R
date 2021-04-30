## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulate our data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# A function that generates 1000 (x,x*,y,z) data points
simdat <- function(beta, gamma, taux, tauy, tauxstr,
                   alpha, zmean, zvar, n1=1000, n2=100) {
  z      = rnorm(n1,zmean, zvar)
  x.true =  alpha[1] + alpha[2]*z + rnorm(n1, 0,1/taux)
  x.partial = c(x.true[1:n2], rep(NA, n1-n2))
  y.true = beta[1] + beta[2]*x.true + beta[3]*z
  
  # Observed y-values
  y  = y.true + rnorm(n1, 0, 1/tauy)
  
  ### Measurement model
  xstr  = x.true + gamma + rnorm(n1,0,1/tauxstr)
  
  output <- list(z, x.true, x.partial, xstr, y.true, y)
  names(output) <- c('z', 'x.true', 'x.partial', 'xstr', 'y.true', 'y')
  return(output)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Bayesian model
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(rjags)
myburnin=500
mysample=2000

### Set up data, initializations, and model
bayes.estimate <- function(s, n1=1000,
                           f="C:/Users/Shannon/Documents/STAT_536D/simulation_2.bug") {
  # Set the parameter values
  jags.data <- list('y'=s$y, 'xstr'=s$xstr, 'x'=s$x.partial, 'z'=s$z,
                    'n1'=n1,
                    'betamean'=c(0,0,0),
                    'betaprec'=diag(rep(0.0001,3)),
                    'gammamean'=c(0,0),
                    'gammaprec'=diag(c(0.0001,0.0001)),
                    'alphamean'=c(0,0),
                    'alphaprec'=diag(c(0.0001,0.0001)))
  
  # Choose which parameters to monitor
  jags.params = c('gamma', 'taux', 'tauxstr', 'tauy', 'beta', 'alpha')
  
  # Keep initializations relatively reasonable
  # Theoretically, we should have good identifiability, so we could use
  # more extreme values for inits but then we'd have to wait longer (more iterations)!
  inits <- list(list(beta=c(0,1,2), alpha=c(2,3), gamma=c(1,-1), taux=3, tauy=3, tauxstr=1),
                list(beta=c(-1,2,-1), alpha=c(-1,0), gamma=c(0,3), taux=1, tauy=2, tauxstr=2),
                list(beta=c(3,3,-3), alpha=c(0,1), gamma=c(3,0), taux=1, tauy=1, tauxstr=1))
  
  m <- jags.model(f, 
                  data=jags.data, n.chains=3, inits=inits)
  
  ### Make a preliminary run of 10000 iterations, with monitoring
  # X <- coda.samples(m, variable.names=jags.params, 
  #                   n.iter=20000, thin=100)
  # This lets us figure out an appropriate burn-in and number of simulations to run
  # Seems like mixing usually occurs by about 500
  # Let's run another 1500 after that for each chain
  
  # Now let's re-run this entire code-chunk (reset the model) to run with burn-in
  burnin     <- coda.samples(m, jags.params, n.iter=myburnin)
  mainsample <- coda.samples(m, jags.params, n.iter=mysample)
  return(mainsample)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Some code so we can nicely display our findings!
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## A function to get the beta-related quantiles
beta.q <- function(mcmc.s, model='m') {
  mcmc.quantiles <- as.data.frame(summary(mcmc.s)[[2]])
  mcmc.quantiles$var <- rownames(mcmc.quantiles)
  mcmc.quantiles$model <- model
  mcmc.quantiles <- mcmc.quantiles[startsWith(mcmc.quantiles$var, 'beta'),c(1,3,5,6,7)]
  names(mcmc.quantiles) <- c('low', 'est', 'high', 'var', 'model')
  return(mcmc.quantiles)
}

## Now some code to output the beta-coefficient estimates from an OLS model
# b1 + b2*x + b3*z
ols.ests <- function(ols.model, model){do.call(rbind, lapply(1:3, function(i){
  
  # Run the actual linear model itself
  ols.summary <- summary(ols.model)$coef
  
  # Extract the three different covariate values
  high <- ols.summary[i,1] + 1.96*ols.summary[i,2]
  low <- ols.summary[i,1] - 1.96*ols.summary[i,2]
  est <- ols.summary[i,1]
  return(data.frame(low=low, est=est, high=high,
                    var=c('beta[1]', 'beta[2]', 'beta[3]')[i],
                    model=model))
}))}

library(ggplot2)
library(reshape2)
# Plot the estimates given a table with format:
# low   est   high   var   model
# "vars" is a list of variable parameters: e.g. 
plot.ests <- function(ests, vars) {
  
  # First, define the label:
  ests$var <- gsub('beta\\[1\\]', 'Intercept', ests$var)
  ests$var <- gsub('beta\\[2\\]', 'beta_X', ests$var)
  ests$var <- gsub('beta\\[3\\]', 'beta_Z', ests$var)
  p <- ggplot(ests, aes(group=paste(var,model), col=model)) +
    
        # Draw a line at zero
        geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
        # Draw a line at one (our true value)
        geom_hline(yintercept = 1, colour = gray(1/2), lty = 2) +
    
        geom_point(aes(x = var, y = est), 
                   position = position_dodge(width = 1/2),
                   size=2) + 
        geom_linerange(aes(x = var, 
                     ymin = low, ymax = high),
                     lwd = 1, position = position_dodge(width = 1/2)) +
        ylab("Coefficient value") + xlab("Coefficient") +
        coord_flip() + theme_bw() +
        theme(legend.position='none')
  return(p)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Define parameters for which we want to simulate data and run the models
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

beta=list(c(1,0,1), c(1,1,1))

# Measurement bias
gamma=list(0,1)

# Measurement error precision
tauxstr=list(0.5,2)

# X=a1+a2*Z
alpha=list(c(1,0), c(1,1))

# Let's get a cross of all possible combinations of our variables
vars=expand.grid(gamma, tauxstr, alpha, beta)
names(vars) <- c("gamma", "tauxstr", "alpha", "beta")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Let's go ahead and run our simulation
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ests.2 <- vector("list", dim(vars)[1])
for (i in 1:dim(vars)[1]) {
  
  # Generate a simulated sample with the desired coefficients
  s <- simdat(beta=vars$beta[[i]], # outcome model:  y=b0+b1*X+b2*Z
              # Measurement model: additive & multiplicative measurement bias
              gamma=vars$gamma[[i]],
              ### Variances (precisions) of true X and Y
              taux=1, tauy=1, 
              # Precision of measurement (accounts for random measurement error)
              tauxstr=vars$tauxstr[[i]],
              ### Exposure model: relating Z to X
              alpha=vars$alpha[[i]],
              zmean=1, zvar=1)
  
  # Run Bayesian models
  bayes_bias <- bayes.estimate(s, n1=1000,
                           f="bayesmod_bias.bug")
  bayes_nobias <- bayes.estimate(s, n1=1000,
                           f="bayesmod_nobias.bug")
  # Run the ordinary least squares
  ols <- lm(s$y~s$xstr + s$z)
  # Run the OLS on the small sample with true values for X
  ols.small <- lm(s$y[1:100]~s$x.partial[1:100] + s$z[1:100])
  
  # Gather the estimates and save them as a table in our ests list!
  ests.2[[i]] <- rbind(ols.ests(ols, 'OLS using X* (n=1000)'), 
              ols.ests(ols.small, 'OLS using true X (n=100)'),
              beta.q(bayes_bias, model='Bayes, measurement bias'),
              beta.q(bayes_nobias, model='Bayes, no measurement bias'))

}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot our estimated coefficients!
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


plots <- lapply(1:16, function(i){plot.ests(ests.2[[i]])})
for (i in 1:16) {
  png(paste0("C:/Users/Shannon/Documents/STAT_536D/coefplots/plot_", i, ".png"))
  
  # Add a quick description of the variables we are working with!
  est.lab <- paste0('X*|X~N(',vars2$alpha[[i]][1],'+',vars2$alpha[[i]][2],'Z,',vars2$tauxstr[[i]],')',
              '\n X|Z~N(X+',vars2$gamma[[i]], ',1)',
              '\n beta_X=', vars2$beta[[i]][2])
  print(plots[[i]] + annotate('text', x=3.3, y=0.25,label=est.lab, size=3))
  dev.off()
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## An aside: don't forget to check traceplots
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(MCMCvis)
MCMCtrace(bayes_bias)

