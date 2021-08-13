genSeason <- function(temps, mrange = c(8, 6), years) {
  stopifnot(length(mrange) == 2)
  N <- length(years)
  stopifnot(N >= 1)
  Season <- list()
  if(inherits(temps, "list"))  temperatures<-temps[[1]] else temperatures<-temps
  for(j in c(1 : N)) {
    ## season starts in November of the preceding year
    ## the winter season from August to June
    ## we add the Year as the last element to ship it to sapply
    Season[[j]] <- c(which(((temperatures$Month %in% c(mrange[1] : 12)) &
                              (temperatures$Year == years[j] - 1))),
                     which(((temperatures$Month %in% c(1 : mrange[2])) &
                              (temperatures$Year == years[j]))))
  }
  return(invisible(Season))
}


# Helper to wrap the DynModel_driver function from chillR

apply_const_temp <- function(temp, A0, A1, E0, E1, Tf, slope, portions = 1200, deg_celsius = TRUE)
  {
    temp_vector <- rep(temp, times = portions)
    res <- chillR::DynModel_driver(temp = temp_vector,
                                   A0 = A0,
                                   A1 = A1,
                                   E0 = E0,
                                   E1 = E1,
                                   Tf = Tf,
                                   slope = slope,
                                   deg_celsius = deg_celsius)
    return(invisible(res$y[length(res$y)]))
}

# Helper to generate a bell shape based on chill response

gen_bell <- function(par, temp_values = seq(-5, 20, 0.1)) {
  E0 <- par[5]
  E1 <- par[6]
  A0 <- par[7]
  A1 <- par[8]
  Tf <- par[9]
  slope <- par[12]
  
  y <- c()
  for(i in seq_along(temp_values)) {
    y[i] <- apply_const_temp(temp = temp_values[i],
                             A0 = A0,
                             A1 = A1,
                             E0 = E0,
                             E1 = E1,
                             Tf = Tf,
                             slope = slope)
  }
  return(invisible(y))
}


# Helper to generate the response to GDH accumulation

GDH_response <- function(par, T){
  
  Tb <- par[11]
  Tu <- par[4]
  Tc <- par[10]
  
  GDH_weight <- rep(0, length(T))
  
  GDH_weight[which(T >= Tb & T <= Tu)] <- 1/2 * (1 + cos(pi + pi * (T[which(T >= Tb & T <= Tu)] - Tb)/(Tu - Tb)))

  GDH_weight[which(T > Tu & T <= Tc)] <- (1 + cos(pi/2 + pi/2 * (T[which(T >  Tu & T <= Tc)] -Tu)/(Tc - Tu)))
  
  return(GDH_weight)
}

