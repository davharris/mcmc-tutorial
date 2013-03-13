## @knitr seed
set.seed(1)


## @knitr show-plot
f = function(x){
  (sin(.2 * x^3) + sin(x^2)/5 + 2.5 + sin(10 * x)/2) * dnorm(x, mean = 1)
}

xs = seq(-6, 6, length = 500)

plot(
  xs, 
  f(xs), 
  type = "n",
  xlim = range(xs),
  ylim = c(0, max(f(xs) + 1/2)),
  xlab = "x",
  ylab = "f(x)",
  yaxs = "i",
  yaxt = "n"
)

polygon(
  c(min(xs), xs, max(xs)), 
  c(0, f(xs), 0), 
  density = 100, 
  lwd = 3, 
  col = "darkblue"
)



## @knitr rejection-sampler
N = 10000000
ceiling = max(f(xs) + 1)

# Step 1: Throw N darts (uniformly distributed x-y pairs) at the plot above.
test.xs = runif(N, min = min(xs), max = max(xs))
test.ys = runif(N, min = 0, max = ceiling)


# Step 2: Throw away any samples that hit white areas (above the curve)
accepted.samples = test.xs[test.ys < f(test.xs)]

# Step 3: Study the samples however you want
hist(accepted.samples, breaks = 150, xlim = range(xs), freq = FALSE)

mean(accepted.samples)

quantile(accepted.samples, probs = seq(0, 1, length = 25))



## @knitr Gaussian
library(MASS)
library(ggplot2)
lik = function(x, y){dnorm(x-3) * dnorm(y - x + 2)}

grid = expand.grid(x = xs, y = xs)
z = lik(grid$x, grid$y)

gaussian.plot = ggplot(data = grid, aes(x = x, y = y)) + geom_raster(aes(fill = z)) + scale_fill_gradient2() + coord_equal()

gaussian.plot


## @knitr MH
maxit = 50

samples = matrix(NA, nrow = maxit, ncol = 2, dimnames = list(NULL, c("x", "y")))
samples[1, ] = c(0,0) # start at 0,0

for(i in 2:maxit){
  
  # propose a new sample point
  proposal = samples[i - 1, ] + rnorm(2, mean = 0, sd = 1)

  # Compare its likelihood with the current position
  old.lik = lik(samples[i - 1, "x"], samples[i - 1, "y"])
  new.lik = lik(proposal["x"], proposal["y"])
  
  ratio = new.lik / old.lik
  
  # flip a coin and accept the new proposal with probability min(ratio, 1)
  if(rbinom(1, size = 1, prob = min(ratio, 1))){
    samples[i, ] = proposal
  }else{
    # If you don't accept the proposal, just keep what you had in the last time step
    samples[i, ] = samples[i - 1, ]
  }
    
}


gaussian.plot + geom_path(mapping = aes(x = samples[,"x"], y = samples[,"y"]), color = "orange") + geom_point(mapping = aes(x = samples[,"x"], y = samples[,"y"]))


## @knitr MH2
maxit = 50

samples = matrix(NA, nrow = maxit, ncol = 2, dimnames = list(NULL, c("x", "y")))
samples[1, ] = c(0,0) # start at 0,0

for(i in 2:maxit){
  
  # propose a new sample point
  proposal = samples[i - 1, ] + rnorm(2, mean = 0, sd = 1)

  # Compare its likelihood with the current position
  old.lik = lik(samples[i - 1, "x"], samples[i - 1, "y"])
  new.lik = lik(proposal["x"], proposal["y"])
  
  ratio = new.lik / old.lik
  
  # flip a coin and accept the new proposal with probability min(ratio, 1)
  if(rbinom(1, size = 1, prob = min(ratio, 1))){
    samples[i, ] = proposal
  }else{
    samples[i, ] = samples[i - 1, ]
  }
    
}


gaussian.plot + geom_path(mapping = aes(x = samples[,"x"], y = samples[,"y"]), color = "orange") + geom_point(mapping = aes(x = samples[,"x"], y = samples[,"y"]))


## @knitr MH3
maxit = 10000

samples = matrix(NA, nrow = maxit, ncol = 2, dimnames = list(NULL, c("x", "y")))
samples[1, ] = c(0,0) # start at 0,0

for(i in 2:maxit){
  
  # propose a new sample point
  proposal = samples[i - 1, ] + rnorm(2, mean = 0, sd = 1)

  # Compare its likelihood with the current position
  old.lik = lik(samples[i - 1, "x"], samples[i - 1, "y"])
  new.lik = lik(proposal["x"], proposal["y"])
  
  ratio = new.lik / old.lik
  
  # flip a coin and accept the new proposal with probability min(ratio, 1)
  if(rbinom(1, size = 1, prob = min(ratio, 1))){
    samples[i, ] = proposal
  }else{
    samples[i, ] = samples[i - 1, ]
  }
    
}


ggplot(data = NULL, mapping = aes(x = samples[,"x"], y = samples[,"y"])) + stat_density2d(geom = "tile", aes(fill = ..density..), contour = FALSE) + scale_fill_gradient2() + xlim(-6,6) + ylim(-6,6) + coord_equal()


