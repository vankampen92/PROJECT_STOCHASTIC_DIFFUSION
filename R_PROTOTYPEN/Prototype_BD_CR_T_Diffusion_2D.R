#############################################################################################
#### 01-12-2020 Blanes
#### second prototype to model diffusion of individuals on a lattice
#### with birth death and immigration
#### and consumer-resource dynamics
#############################################################################################


#############################################################################################
#### set language seed and other settings
#### load packages
#############################################################################################

if ( !require("rlist") )   { install.packages("rlist"); library("rlist") }
if ( !require("deSolve") ) { install.packages("deSolve"); library("deSolve") }

rm(list=ls())
Sys.setenv(LANG = "en")
set.seed(12345)

PLOT_IN_TERMINAL <- F
COMPUTE_AVERAGES <- F
PLOT_STEADY      <- T
#############################################################################################
#### set spatial lattice, cells and individuaks
#############################################################################################

#### size of the lattice
L <- 1
#### number of cells
M <- L*L
#### initial number of resource individuals
N_R  <- 0
#### initial number of consumer individuals
N_A  <- 0
#### initial number of compounds
N_C  <- 0
#### initial number of triplets
N_T  <- 0

#############################################################################################
#### model parameters: Rates of the process 
#############################################################################################

#### rate of movement of resources
mu_R     <- 0
#### rate of immigration of resources
lambda_R <- 50
#### death rate of resources
delta_R  <- 1
####lamda_R/delta_R is the average numebr of resource units
#### rate of movement of consumer A
mu_A     <- 0
#### rate of immigration of consumer A
lambda_A <- 1
#### death rate of consumer A
delta_A  <- 1
#### encounter rate between consumer and resource
alpha    <- 1
#### system's Volume
V <- 10
#### handling time
h <- 1
#### "decomposition rate" of consumer-resource pairs (compounds)
nu <- 1/h
#### rate of triplet formation
chi <- 1
#### rate of triplet decomposition
eta <- 100

#############################################################################################
#### stochastic simulation settings
#############################################################################################

#### number of events
NEV    <- 7000
### number of replicates
NREP   <- 40

#############################################################################################
#### start dynamics
#############################################################################################

#### definelist to put results of simulations
list_simu <- list()
####loop in replicates
kk <- 1
for(kk in 1:NREP){
  list_simu[[kk]] <- list()
  #### initial time
  tt <- 0
  #### define a matrix to describe resources in the lattice (with all zeros)
  spatial_grid_R <-  matrix(nrow = L, ncol = L, 0)
  #### define a matrix to describe consumers in the lattice (with all zeros)
  spatial_grid_A <-  matrix(nrow = L, ncol = L, 0)
  #### define a matrix to describe compounds in the lattice (with all zeros)
  spatial_grid_C <-  matrix(nrow = L, ncol = L, 0)
  #### define a matrix to describe triplets in the lattice (with all zeros)
  spatial_grid_T <-  matrix(nrow = L, ncol = L, 0)
  
  #### set initial conditions
  #### put all individuals at the center of the lattice
  spatial_grid_R[ceiling(L/2),ceiling(L/2)] <- N_R
  spatial_grid_A[ceiling(L/2),ceiling(L/2)] <- N_A
  spatial_grid_C[ceiling(L/2),ceiling(L/2)] <- N_C
  spatial_grid_T[ceiling(L/2),ceiling(L/2)] <- N_T

  ##### transition rate at beginning
  # R_movement   <- (mu_R+delta_R)*sum(spatial_grid_R)+(mu_A+delta_A)*sum(spatial_grid_A)
  # R_immigration<- (lambda_R+lambda_A)*M
  # R_feeding    <- (alpha/V)*sum(spatial_grid_A*spatial_grid_R) + nu*sum(spatial_grid_C)
  # 
  # R <- R_movement+R_immigration+R_feeding
  
  
  #### loop in events
  jj <- 1
  for(jj in 1:NEV){
    #### calculate the rate of events
    R_movement   <- (mu_R+delta_R)*sum(spatial_grid_R)+(mu_A+delta_A)*sum(spatial_grid_A)
    R_immigration<- (lambda_R+lambda_A)*M
    R_feeding    <- (alpha/V)*sum(spatial_grid_A*spatial_grid_R) + nu*sum(spatial_grid_C)
    R_feeding    <- R_feeding + chi*sum(spatial_grid_A*spatial_grid_C) + eta*sum(spatial_grid_T)
      
    R <- R_movement+R_immigration+R_feeding
    
    #### draw time step for next event
    pp      <- runif(1,0,1)
    delta_t <- -log(pp)/R
    tt      <- tt + delta_t
    
    #### calculate matrix of rates
    spatial_grid_rates <- (mu_R+delta_R)*spatial_grid_R+(mu_A+delta_A)*spatial_grid_A+lambda_R+lambda_A
    spatial_grid_rates <- spatial_grid_rates + (alpha/V)*spatial_grid_A*spatial_grid_R + nu*spatial_grid_C
    spatial_grid_rates <- spatial_grid_rates + chi*spatial_grid_A*spatial_grid_C + eta*spatial_grid_T

    #### normalize to obtain the matrix of probabilities
    spatial_grid_pp    <- spatial_grid_rates/sum(spatial_grid_rates)
    
    #### pick a site where the event happens
    repeat{
      rand_col <- sample(1:ncol(spatial_grid_pp),1)
      rand_row <- sample(1:nrow(spatial_grid_pp),1)
      pp <- runif(1,0,1)
      if(spatial_grid_pp[rand_row,rand_col] > pp)break
    }
    #### check the same thing with sample (probabilities from matrix)
    #### site has been selected
    
    #### probabilities of the 10 events 
    #### movement
    p_m_R    <- mu_R*spatial_grid_R[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_m_A    <- mu_A*spatial_grid_A[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    #### immigration
    p_i_R    <- lambda_R/spatial_grid_rates[rand_row,rand_col]
    p_i_A    <- lambda_A/spatial_grid_rates[rand_row,rand_col]
    #### deaths
    p_d_R    <- delta_R*spatial_grid_R[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_d_A    <- delta_A*spatial_grid_A[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    #### feeding
    p_e      <- (alpha/V)*spatial_grid_R[rand_row,rand_col]*spatial_grid_A[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_d      <- nu*spatial_grid_C[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    
    p_t_f    <- chi*spatial_grid_A[rand_row,rand_col]*spatial_grid_C[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]
    p_t_d    <- eta*spatial_grid_T[rand_row,rand_col]/spatial_grid_rates[rand_row,rand_col]

    #### check normalization to 1
    p_m_R+p_m_A+p_i_R+p_i_A+p_d_R+p_d_A+p_e+p_d+p_t_f+p_t_d 
    
    #### now select the event that happens in the site
    pp <- runif(1,0,1)
    
    #### immigration of resources
    p1 <- p_i_R
    if(pp<p1){
      spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]+1
    }
    #### immigration of consumers
    p2 <- p1 + p_i_A
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]+1
    }
    #### death of resources
    p1 <- p2
    p2 <- p1 + p_d_R
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]-1
    }
    #### death of consumer
    p1 <- p2
    p2 <- p1 + p_d_A
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]-1
    }
    #### formation of compund through encounter
    p1 <- p2
    p2 <- p1 + p_e
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_R[rand_row,rand_col]<- spatial_grid_R[rand_row,rand_col]-1
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]-1
      spatial_grid_C[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]+1
    }
    #### decomposition or handling
    p1 <- p2
    p2 <- p1 + p_d
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_C[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]-1
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]+2
    }
    #### formation of triplets
    p1 <- p2
    p2 <- p1 + p_t_f
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_C[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]-1
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]-1
      spatial_grid_T[rand_row,rand_col]<- spatial_grid_T[rand_row,rand_col]+1
    }
    #### decompostion of triplets
    p1 <- p2
    p2 <- p1 + p_t_d
    if((pp>=p1)&&(pp<p2)){
      spatial_grid_T[rand_row,rand_col]<- spatial_grid_T[rand_row,rand_col]-1
      spatial_grid_C[rand_row,rand_col]<- spatial_grid_C[rand_row,rand_col]+1
      spatial_grid_A[rand_row,rand_col]<- spatial_grid_A[rand_row,rand_col]+1
    }
    p1 <- p2
    p2 <- p1 + p_m_R
    if((pp>=p1)&&(pp<p2)){      
      
      ##### remove one individual from the selected site
      spatial_grid_R[rand_row,rand_col] <- spatial_grid_R[rand_row,rand_col]-1
      
      #### select a new site where the individual moves
      pp <- runif(1,0,1)
      
      #### go down
      if(pp<0.25){
        next_row <- rand_row - 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
      #### go up
      if((pp>=0.25)&&(pp<0.5)){
        next_row <- rand_row + 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
      #### go left
      if((pp>=0.5)&&(pp<0.75)){
        next_row <- rand_row
        next_col <- rand_col - 1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
      #### go right
      if(pp>=0.75){
        next_row <- rand_row
        next_col <- rand_col+1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_R[next_row,next_col] <- spatial_grid_R[next_row,next_col]+1
      }
    }
    #### movement of consumers
    if((pp>=p2)){
      
      ##### remove one individual from the selected site
      spatial_grid_A[rand_row,rand_col] <- spatial_grid_A[rand_row,rand_col]-1
      
      #### select a new site where the individual moves
      pp <- runif(1,0,1)
      
      #### go down
      if(pp<0.25){
        next_row <- rand_row - 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
      #### go up
      if((pp>=0.25)&&(pp<0.5)){
        next_row <- rand_row + 1
        next_col <- rand_col
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
      #### go left
      if((pp>=0.5)&&(pp<0.75)){
        next_row <- rand_row
        next_col <- rand_col - 1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
      #### go right
      if(pp>=0.75){
        next_row <- rand_row
        next_col <- rand_col+1
        
        ## periodic boundary conditions
        if(next_col == L+1)next_col <- 1
        if(next_col == 0)  next_col <- L
        if(next_row == L+1)next_row <- 1
        if(next_row == 0)  next_row <- L
        
        spatial_grid_A[next_row,next_col] <- spatial_grid_A[next_row,next_col]+1
      }
    }
    
    #############################################################################################
    #### save results in a list
    list_simu[[kk]][[jj]] <- list(tt,spatial_grid_R,spatial_grid_A,
                                     spatial_grid_C,spatial_grid_T)
    
    #### plot results in terminal if needed
    if(PLOT_IN_TERMINAL){
      print(paste("replica ",kk,sep=""))
      print(paste("event "  ,jj,sep=""))
      print(paste("time "   ,tt,sep=""))
      print("Population Matrix Resources")
      print(spatial_grid_R)
      print(paste("Total number of resources ", sum(spatial_grid_R),sep=""))
      print("Population Matrix Consumers")
      print(spatial_grid_A)
      print(paste("Total number of consumers ", sum(spatial_grid_A),sep=""))
      print("Population Matrix Compounds")
      print(spatial_grid_C)
      print(paste("Total number of compounds ", sum(spatial_grid_C),sep=""))
      print("Population Matrix Triples")
      print(spatial_grid_T)
      print(paste("Total number of triplets ", sum(spatial_grid_T),sep=""))
    }
  }###   close loop in events
  print(kk)
}### close loop in replicates

#### save list
list.save(list_simu,file="list_simu.rdata")

#############################################################################################
### sampling of the results of stochastic simulation 
#############################################################################################

### time step
delta_t <- 0.01

### load list with simulation results
list_results <- list.load("list_simu.rdata")
## select site
rand_row <- ceiling(L/2)
rand_col <- ceiling(L/2)

### data frame to put elements from the replicates
dd_R           <- data.frame(matrix(nrow = NEV,ncol = 2*NREP,0))
dd_A           <- data.frame(matrix(nrow = NEV,ncol = 2*NREP,0))
dd_C           <- data.frame(matrix(nrow = NEV,ncol = 2*NREP,0))
dd_T           <- data.frame(matrix(nrow = NEV,ncol = 2*NREP,0))

colnames(dd_R) <- c(rep(paste("time_",seq(1,NREP,1),sep="")),
                    rep(paste("pop_", seq(1,NREP,1),sep="")))
colnames(dd_A) <- colnames(dd_R)
colnames(dd_C) <- colnames(dd_R)
colnames(dd_T) <- colnames(dd_R)

### extract time series from one site and put them in the matrix
kk <- 1; jj <- 1
for(kk in 1:NREP){for(jj in 1:NEV){
  dd_R[jj,paste("time_",kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[1]][1])
  dd_R[jj,paste("pop_", kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[2]][rand_row,rand_col])
  dd_A[jj,paste("time_",kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[1]][1])
  dd_A[jj,paste("pop_", kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[3]][rand_row,rand_col])
  dd_C[jj,paste("time_",kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[1]][1])
  dd_C[jj,paste("pop_", kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[4]][rand_row,rand_col])
  dd_T[jj,paste("time_",kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[1]][1])
  dd_T[jj,paste("pop_", kk,sep="")]<-as.numeric(list_results[[kk]][[jj]][[5]][rand_row,rand_col])
}
  print(paste("extract_rep_",kk,sep=""))
}

max_time   <- min(min(dd_R[dim(dd_R)[1],1:NREP]),
                  min(dd_A[dim(dd_A)[1],1:NREP]),
                  min(dd_C[dim(dd_C)[1],1:NREP]),
                  min(dd_T[dim(dd_C)[1],1:NREP]))
max_pop_R  <- max(dd_R[1:dim(dd_R)[1],(1+NREP):(2*NREP)])
max_pop_A  <- max(dd_A[1:dim(dd_A)[1],(1+NREP):(2*NREP)])
max_pop_C  <- max(dd_C[1:dim(dd_C)[1],(1+NREP):(2*NREP)])
max_pop_T  <- max(dd_T[1:dim(dd_T)[1],(1+NREP):(2*NREP)])

#############################################################################################
### compute averages of the results of stochastic simulations
#############################################################################################

time_steps   <- seq(delta_t,max_time ,delta_t)
if(COMPUTE_AVERAGES){
  averages_R   <- data.frame(matrix(nrow=length(time_steps),ncol=3,0))
  averages_A   <- data.frame(matrix(nrow=length(time_steps),ncol=3,0))
  averages_C   <- data.frame(matrix(nrow=length(time_steps),ncol=3,0))
  averages_T   <- data.frame(matrix(nrow=length(time_steps),ncol=3,0))
  colnames(averages_R) <-c("time","mean","sd")
  colnames(averages_A) <-colnames(averages_R)
  colnames(averages_C) <-colnames(averages_R)
  colnames(averages_T) <-colnames(averages_R)
  
  #ii <- 100*delta_t
  counter <- 0
  for(ii in time_steps){
    counter <- counter + 1 
    
    mm <- numeric(NREP)
    for(jj in 1:NREP){
      
      bb <- which((dd_R[,jj]<=(ii)))
      aa <- which((dd_R[bb,jj]>(ii-delta_t)))
      # if(length(aa)==0){
      #   kk <- 1
      #   repeat{
      #     kk <- kk+1
      #     aa <- which((dd_R[bb,jj-1]>(ii-kk*delta_t)))
      #     if(length(aa>0))break
      #   }
      # }
      # print(paste("average_resource_time_",jj,sep=""))
      # mm[jj] <- dd_R[aa[length(aa)],jj+NREP]
      
      if(length(aa)==0) mm[jj] <- 0
      if(length(aa)>0)  mm[jj] <- dd_R[aa[length(aa)],jj+NREP]
    }
    
    #print(paste("average_resource_time_",ii,sep=""))
    
    averages_R[counter,"time"] <- ii
    averages_R[counter,"mean"] <- mean(mm)
    averages_R[counter,"sd"]   <- sd(mm)
  }
  counter <- 0
  for(ii in time_steps){
    counter <- counter + 1 
    mm <- numeric(NREP)
    
    for(jj in 1:NREP){
      
      bb <- which((dd_A[,jj]<=(ii)))
      aa <- which((dd_A[bb,jj]>(ii-delta_t)))
      # if(length(aa)==0){
      #   kk <- 1
      #   repeat{
      #     kk <- kk+1
      #     aa <- which((dd_A[bb,jj-1]>(ii-kk*delta_t)))
      #     if(length(aa>0))break
      #   }
      # }
      # mm[jj] <- dd_A[aa[length(aa)],jj+NREP]
      
      if(length(aa)==0)mm[jj] <- 0
      if(length(aa)>0) mm[jj] <- dd_A[aa[length(aa)],jj+NREP]
    }
    #print(paste("average_consumer_time_",ii,sep=""))
    
    averages_A[counter,"time"] <- ii
    averages_A[counter,"mean"] <- mean(mm)
    averages_A[counter,"sd"]   <- sd(mm)
  }
  counter <- 0
  for(ii in time_steps){
    counter <- counter + 1 
    
    mm <- numeric(NREP)
    
    for(jj in 1:NREP){
      
      bb <- which((dd_C[,jj]<=(ii)))
      aa <- which((dd_C[bb,jj]>(ii-delta_t)))
      # if(length(aa)==0){
      #   kk <- 1
      #   repeat{
      #     kk <- kk+1
      #     aa <- which((dd_C[bb,jj-1]>(ii-kk*delta_t)))
      #     if(length(aa>0))break
      #   }
      # }
      # mm[jj] <- dd_C[aa[length(aa)],jj+NREP]
      
      if(length(aa)==0)mm[jj] <- 0
      if(length(aa)>0) mm[jj] <- dd_C[aa[length(aa)],jj+NREP]
    }
    #print(paste("average_compound_time_",ii,sep=""))
    
    averages_C[counter,"time"] <- ii
    averages_C[counter,"mean"] <- mean(mm)
    averages_C[counter,"sd"]   <- sd(mm)
  }
  counter <- 0
  for(ii in time_steps){
    counter <- counter + 1 
    
    mm <- numeric(NREP)
    
    for(jj in 1:NREP){
      
      bb <- which((dd_T[,jj]<=(ii)))
      aa <- which((dd_T[bb,jj]>(ii-delta_t)))
      # if(length(aa)==0){
      #   kk <- 1
      #   repeat{
      #     kk <- kk+1
      #     aa <- which((dd_C[bb,jj-1]>(ii-kk*delta_t)))
      #     if(length(aa>0))break
      #   }
      # }
      # mm[jj] <- dd_C[aa[length(aa)],jj+NREP]
      
      if(length(aa)==0)mm[jj] <- 0
      if(length(aa)>0) mm[jj] <- dd_T[aa[length(aa)],jj+NREP]
    }
    #print(paste("average_compound_time_",ii,sep=""))
    
    averages_T[counter,"time"] <- ii
    averages_T[counter,"mean"] <- mean(mm)
    averages_T[counter,"sd"]   <- sd(mm)
  }
}

#############################################################################################
##### compute deterministic trajectories
#############################################################################################

mean_field <- function(t,state,parms) {
  with(as.list(c(state,parms)),{
    dR <- lambda_R - delta_R*R - (alpha/V)*R*A
    dA <- lambda_A - delta_A*A + 2*nu*C  - (alpha/V)*R*A - chi*C*A + eta*TT
    dC <- (alpha/V)*R*A - nu*C - chi*C*A + eta*TT
    dTT<- chi*C*A - eta*TT
    return(list(c(dR,dA,dC,dTT)))
  })
}

time  <- seq(0,max_time,delta_t)
state <- c(R=N_R,A=N_A,C=N_C,TT=N_T)
parms <- c(lambda_R = lambda_R,
           lambda_A = lambda_A,
           delta_R  = delta_R,
           delta_A  = delta_A,
           alpha    = alpha,
           V        = V,
           nu       = nu,
           chi      = chi,
           eta      = eta
)

det <- rk4(state, time, mean_field, parms)
colnames(det) <- c("time","pop_R","pop_A","pop_C","pop_T")



steady_state <- function(parms){
  
  lambda_A <- parms["lambda_A"]
  lambda_R <- parms["lambda_R"]
  delta_A  <- parms["delta_A"]
  delta_R  <- parms["delta_R"]
  alpha    <- parms["alpha"]/parms["V"]
  nu       <- parms["nu"]
  chi      <- parms["chi"]
  eta      <- parms["eta"]
  
  DELTA_R <- ( alpha*(lambda_A + lambda_R) + delta_A*delta_R)^2 - 4*alpha*lambda_R*delta_A*delta_R
  R_1_R   <- ( (alpha*(lambda_A + lambda_R) + delta_A*delta_R ) + sqrt(DELTA_R) )/(2*alpha*delta_R)
  R_2_R   <- ( (alpha*(lambda_A + lambda_R) + delta_A*delta_R ) - sqrt(DELTA_R) )/(2*alpha*delta_R)
  A_1_R   <- ( lambda_A + lambda_R - delta_R*R_1_R )/delta_A
  A_2_R   <- ( lambda_A + lambda_R - delta_R*R_2_R )/delta_A
  C_1_R   <- alpha*R_1_R*A_1_R/nu 
  C_2_R   <- alpha*R_2_R*A_2_R/nu 
  T_1_R   <- chi*C_1_R*A_1_R/eta 
  T_2_R   <- chi*C_2_R*A_2_R/eta 
  
  DELTA_R2<- ( alpha*(lambda_A - lambda_R) + delta_A*delta_R)^2 + 4*alpha*alpha*lambda_A*lambda_R
  
  DELTA_A <- ( alpha*(lambda_A + lambda_R) - delta_A*delta_R)^2 + 4*alpha*lambda_A*delta_A*delta_R
  A_1_A   <- ( (alpha*(lambda_A + lambda_R) - delta_A*delta_R ) - sqrt(DELTA_A) )/(2*alpha*delta_A)
  A_2_A   <- ( (alpha*(lambda_A + lambda_R) - delta_A*delta_R ) + sqrt(DELTA_A) )/(2*alpha*delta_A)
  R_1_A   <- ( lambda_A + lambda_R - delta_A*A_1_A )/delta_R
  R_2_A   <- ( lambda_A + lambda_R - delta_A*A_2_A )/delta_R
  C_1_A   <- alpha*R_1_A*A_1_A/nu 
  C_2_A   <- alpha*R_2_A*A_2_A/nu 
  T_1_A   <- chi*C_1_A*A_1_A/eta 
  T_2_A   <- chi*C_2_A*A_2_A/eta 
  
  DELTA_A2<- ( alpha*(lambda_A - lambda_R) + delta_A*delta_R)^2 + 4*alpha*alpha*lambda_A*lambda_R
  
  const_1 <- (lambda_R - sqrt(lambda_A*lambda_R) - lambda_A)/2
  const_2 <- delta_A*delta_R/alpha
  const_3 <- (lambda_R + sqrt(lambda_A*lambda_R) - lambda_A)/2
  
  all_results       <- c(DELTA_R, DELTA_A, DELTA_A2,const_1,const_2,const_3,
                         R_1_R,R_2_R,A_1_R,A_2_R,C_1_R,C_2_R,T_1_R,T_2_R,
                         R_1_A,R_2_A,A_1_A,A_2_A,C_1_A,C_2_A,T_1_A,T_2_A)
  names(all_results)<- c("DELTA_R","DELTA_A","DELTA_A2","Const_1","Const_2","Const_3",
                         "R_1_R","R_2_R","A_1_R","A_2_R","C_1_R","C_2_R","T_1_R","T_2_R",
                         "R_1_A","R_2_A","A_1_A","A_2_A","C_1_A","C_2_A","T_1_A","T_2_A")
  return(all_results)
}
steady_state(parms)

#############################################################################################
### plotting stochastic trajectories and add deterministic lines
#############################################################################################

PP  <- 4

b_b <- PP
l_l <- PP
t_t <- PP
r_r <- PP

X11(width=200,height=200)
par(oma=c(6,3,6,3))

layout.matrix <- matrix(seq(1,4,1),
                        nrow = 1, ncol = 4,
                        byrow = F)

layout(mat = layout.matrix,
       heights = rep(1,1),   # Heights of rows
       widths  = rep(1,1))   # Widths of columns

layout.show(4)

max_pop <- max(max_pop_R,max_pop_A,max_pop_C,max_pop_T)

par(mar = c(b_b,l_l,t_t,r_r))
plot(0,0,
     ylab="individuals",
     xlab="Time",
     main="RESOURCES",
     ylim=c(0,max_pop),
     xlim=c(0,max_time))

for(ii in 1:NREP) lines(dd_R[,paste("time_",ii,sep="")],
                        dd_R[,paste("pop_" ,ii,sep="")],type='l',col="red")

if(COMPUTE_AVERAGES)lines(averages_R[,"time"],averages_R[,"mean"], col= "blue",lwd=3)
lines(det[,"time"],det[,"pop_R"],col="black",lwd=3)
if(PLOT_STEADY){
abline(h=steady_state(parms)["R_1_R"],col="blue")
abline(h=steady_state(parms)["R_2_R"],col="green")
}

par(mar = c(b_b,l_l,t_t,r_r))
plot(0,0,
     ylab="individuals",
     xlab="Time",
     main="CONSUMERS",
     ylim=c(0,max_pop),
     xlim=c(0,max_time))

for(ii in 1:NREP) lines(dd_A[,paste("time_",ii,sep="")],
                        dd_A[,paste("pop_" ,ii,sep="")],type='l',col="red")

if(COMPUTE_AVERAGES)lines(averages_A[,"time"],averages_A[,"mean"], col= "blue",lwd=3)
lines(det[,"time"],det[,"pop_A"],col="black",lwd=3)
if(PLOT_STEADY){
abline(h=steady_state(parms)["A_1_R"],col="blue")
abline(h=steady_state(parms)["A_2_R"],col="green")
}

par(mar = c(b_b,l_l,t_t,r_r))
plot(0,0,
     ylab="individuals",
     xlab="Time",
     main="COMPOUNDS",
     ylim=c(0,max_pop),
     xlim=c(0,max_time))

for(ii in 1:NREP) lines(dd_C[,paste("time_",ii,sep="")],
                        dd_C[,paste("pop_" ,ii,sep="")],type='l',col="red")

if(COMPUTE_AVERAGES)lines(averages_C[,"time"],averages_C[,"mean"], col= "blue",lwd=3)
lines(det[,"time"],det[,"pop_C"],col="black",lwd=3)
if(PLOT_STEADY){
abline(h=steady_state(parms)["C_1_R"],col="blue")
abline(h=steady_state(parms)["C_2_R"],col="green")
}

par(mar = c(b_b,l_l,t_t,r_r))
plot(0,0,
     ylab="individuals",
     xlab="Time",
     main="TRIPLETS",
     ylim=c(0,max_pop),
     xlim=c(0,max_time))

for(ii in 1:NREP) lines(dd_T[,paste("time_",ii,sep="")],
                        dd_T[,paste("pop_" ,ii,sep="")],type='l',col="red")

if(COMPUTE_AVERAGES)lines(averages_T[,"time"],averages_T[,"mean"], col= "blue",lwd=3)
lines(det[,"time"],det[,"pop_T"],col="black",lwd=3)
if(PLOT_STEADY){
abline(h=steady_state(parms)["T_1_R"],col="blue")
abline(h=steady_state(parms)["T_2_R"],col="green")
}

mtext(paste("PARAMETERS\n",
            " Lambda_R = ",lambda_R," Lambda_A = ",lambda_A," Delta_R = ",delta_R," Delta_A = ",delta_A,"\n",
            " Alpha/V = ",alpha,"/",V," = ",alpha/V," Nu = ", nu,"\n",
            " Chi = ",chi," Eta = ", eta,
            sep=""),
      side=3,outer=T,at=0.5,cex=1,line=0)
