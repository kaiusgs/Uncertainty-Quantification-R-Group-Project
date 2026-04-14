# ******************************
# This should be a single R script which reproduces the best fit parameters
# and all plots in the poster, and must run without any changes needed. 
# 
# The code should be tidy and well commented on what each aspect is implementing. 
# 
# Clearly indicate which element of the poster or which figure each part of 
# the code is producing. 
# 
# All code should be understood by all group members.
# 
# Group Number  |  Model  |  High/low density  
#      30       |   III   |        low
# ******************************


# ******************************
# Existing Data
# ******************************

# estimated experimental values of data points in Figure 2
time_data <- c(0, 2, 4, 5, 6, 7, 12) # Time(day-degrees/100)
S_data    <- c(0.0, 5.4, 3.8, 3.9, 3.9, 3.5, 3.3) # Susceptible stems
I_data    <- c(0.0, 0.1, 0.3, 0.3, 0.5, 0.5, 0.2) # Infected stems

# best-fit parameters for Model III Low density in Table 2
b       <- 0.0005
kappa   <- 11.725
lambda0 <- 0.064
mu      <- 0.140
alpha1  <- 0.721
alpha2  <- 0.261
alpha3  <- 0.0003
d       <- 0.215
# corresponding performance 
RSS     <- 10.73
df      <- 9


# ******************************
# Task 1 
# For your prescribed Model of f(I, S), solve the ODEs numerically in R with
# zero as initial conditions using the best-fit values for the parameters given
# in Table 2 of the article.
# ******************************

# library for solving differential equations
library(deSolve) 

# define the ODE system function
model3 = function(t, y, pars){
    with(as.list(c(y, pars)), {
        # Model III System 
        dS = b * (kappa - (S + I + R)) - L * S - 
            (alpha1 * I^2 - alpha2 * I) / (alpha3 + I^2)
        dI = L * S - d * I 
        dR = d * I 
        dL = - mu * L 

        return(list(c(dS, dI, dR, dL)))
    })
}

# define the parameter using the best-fit values
pars = c(
    b       = b       ,
    kappa   = kappa   ,
    lambda0 = lambda0 ,
    mu      = mu      ,
    alpha1  = alpha1  ,
    alpha2  = alpha2  ,
    alpha3  = alpha3  ,
    d       = d       
)

# set the initial condition for y 
yini = c(
    S = 0,
    I = 0,
    R = 0,
    L = lambda0 # lambda
)

# define the times to store the solution
times <- seq(0, 12, by = 0.01)

# Call the ODE solver and print the summary of the solution
out1  <- ode(yini, times, model3, pars)
summary(out1)

# plot numerical solution of the ODE system 
# plot the predicted susceptible stems
plot(
  out1[,1], out1[,2], type = "l", 
  xlim = c(0, 12), ylim = c(0, 6), 
  xlab = "Time(day-degrees/100)", ylab = "Susceptible stems", 
  main = "Model III with Best-fit Parameters"
)
points(time_data, S_data, pch = 17, cex = 2)
# plot the predicted infected stems
plot(
  out1[,1], out1[,3], type = "l", 
  xlim = c(0, 12), ylim = c(0, 6), 
  xlab = "Time(day-degrees/100)", ylab = "Infected stems", 
  main = "Model III with Best-fit Parameters"
)
points(time_data, I_data, pch = 17, cex = 2)


# ******************************
# Task 2 
# With your estimates of the data points, determine the best-fit for the 
# parameters with f(I, S) for your Model. 
# Compare your best-fit parameters with Table 2.
# Include the RSS errors for your fits (high and low density for both S and I). 
# Plot the solution of the ODEs with (simultaneously) your best-fit values, 
# the article’s best-fit values, and also your data points.
# ******************************

# library for least squares fit using levenberg-marquart algorithm
library(minpack.lm) 
# library for reshaping data (tall-narrow <-> short-wide)
library(reshape2) 

# preprocess for experimental data 
experiment_df = data.frame(time_data, S_data, I_data)
names(experiment_df) = c("time", "S", "I")

ssq = function(pars){
    # define the parameter, fixing b kappa mu & alpha3
    parameters = c(
        b       = b                      ,
        kappa   = kappa                  ,
        lambda0 = unname(pars['lambda0']),
        mu      = mu                     ,
        alpha1  = unname(pars['alpha1']) ,
        alpha2  = unname(pars['alpha2']) ,
        alpha3  = unname(pars['alpha3']) ,
        d       = unname(pars['d'])      
    )
  
    # initial values for y 
    yini = c(
        S = 0,
        I = 0,
        R = 0,
        L = unname(pars['lambda0'])
    )

    # define the times which have already included the experimental time points
    # times <- c(seq(0, 12, by = 0.01), experiment_df$time)
    # times = sort(unique(times)) 
    times <- seq(0, 12, by = 0.01)

    # solve ODE for a given set of parameters
    out  <- ode(yini, times, model3, parameters)

    # Filter data that contains experimental time points
    out_df = data.frame(out)
    out_df = out_df[out_df$time %in% time_data, c("time", "S", "I")]
    
    # Evaluate predicted vs experimental residual
    pred_df = melt(out_df, id.var="time", variable.name="y", value.name="stems")
    exp_df = melt(experiment_df, id.var="time", variable.name="y",value.name="stems")
    ssqres = (pred_df$stems - exp_df$stems)

    # return predicted vs experimental residual
    return(ssqres)
    
}

# parameter fitting using levenberg marquart algorithm
# initial guess for parameters
init_pars = c(
    lambda0 = 0.05,
    alpha1  = 0.8,
    alpha2  = 0.3,
    alpha3  = 0.001,
    d       = 0.2
)
# fit and print the summary
fit_val = nls.lm(par = init_pars, fn = ssq)
summary(fit_val)

# print fitted parameter from estimated data points
# ⁠Best-fit values: λ₀, α₁, α₂, d
par_est = as.list(coef(fit_val))
print(par_est)

# setting the fitted parameters and initial condition of y 
esti_pars = c(
    b       = b              ,
    kappa   = kappa          ,
    lambda0 = par_est$lambda0,
    mu      = mu             ,
    alpha1  = par_est$alpha1 ,
    alpha2  = par_est$alpha2 ,
    alpha3  = par_est$alpha3 ,
    d       = par_est$d       
)
esti_yini = c(
    S = 0,
    I = 0,
    R = 0,
    L = par_est$lambda0
)

# define the times 
times <- seq(0, 12, by = 0.01)

# solve ODE using the fitted parameters and print the summary
out2  <- ode(esti_yini, times, model3, esti_pars)
summary(out2)

# Calculate the Total RSS and ⁠RSS for S and I separately
# Filter data that contains experimental time points
out2_df = data.frame(out2)
out2_df = out2_df[out2_df$time %in% time_data, c("time", "S", "I")]
# Evaluate predicted vs experimental residual
pred_df = melt(out2_df, id.var="time", variable.name="y", value.name="stems")
exp_df = melt(experiment_df, id.var="time", variable.name="y",value.name="stems")
print(paste("Total RSS:",  sum((pred_df$stems - exp_df$stems)^2) ))
print(paste("RSS for S:",  sum((out2_df$S - experiment_df$S)^2) ))
print(paste("RSS for I:",  sum((out2_df$I - experiment_df$I)^2) ))

# ⁠Predicted S and I at each time point (0,2,4,5,6,7,12)
print(out2_df)

# Plot the solution of the ODEs with our best-fit values, 
# the article’s best-fit values, and also your data points.
# plot the predicted susceptible stems
plot(
  out1[,1], out1[,2], type = "l", col = "blue", lwd = 2.5,
  xlim = c(0, 12), ylim = c(0, 6), 
  xlab = "Time(day-degrees/100)", ylab = "Susceptible stems", 
  main = "Model III"
)
lines(out2[,1], out2[,2], type = "l", col = "green", lwd = 2.5)
points(time_data, S_data, pch = 17, cex = 2)
legend("topright", legend = c( "papaer's best-fit","our best-fit"),
    col = c("blue", "green"), lwd = 2.5)
# plot the predicted infected stems
plot(
  out1[,1], out1[,3], type = "l", col = "blue", lwd = 2.5,
  xlim = c(0, 12), ylim = c(0, 6), 
  xlab = "Time(day-degrees/100)", ylab = "Infected stems", 
  main = "Model III"
)
lines(out2[,1], out2[,3], type = "l", col = "green", lwd = 2.5)
points(time_data, I_data, pch = 17, cex = 2)
legend("topright", legend = c( "papaer's best-fit","our best-fit"),
    col = c("blue", "green"), lwd = 2.5)


# ******************************
# Task 3 
# Sensitivity Analysis
# ******************************


