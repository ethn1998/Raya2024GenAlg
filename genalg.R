
#Initialization
set.seed(20101019)
maxpopsize = 100 #Population size
ngen = 100 #Steady state is reached at 60 generations
target <- "Selamat Hari Raya!"
top <- 6 #Prints top x most fit phenotypes. Obsolete
#target0 <-"Selamat Hari Jadi!" #Small environmental change
#target1 <-"Happy Birthday Ali" #Large environmental change
#target <- "Hello World!"

indivlength <- nchar(target) #Length of string phenotype
glist <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z",
           "A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
           ".",",","?","!",
           " ") #All possible alleles

selstrength <- 10.0 #Selection strength
mut_rate0 <- 1.0/indivlength #Probability of a given element being mutated

# Sanity checks?

Nforms <- length(glist) #Number of all possible alleles
max_entropy <- indivlength*log(Nforms)


# This is only for visualization
get_trait_counts <- function(pop,index){ #Counts number of alleles
  gcounts <- c() #Emptylist
  for(g in glist){
    gcount <- 0
    for(indiv in pop){
      if(g==substr(indiv,index,index)){
        gcount = gcount + 1
      }
    }
    gcounts <- append(gcounts,gcount)
  }
  return(gcounts)
}

get_mode_phenotype <- function(pop){
  mode_phenotype <- ""
  for(trait in 1:indivlength){
    gcounts <- get_trait_counts(pop,trait)
    modetraitindex <- which.max(gcounts)
    mode_phenotype <- paste(mode_phenotype,glist[modetraitindex],sep="")
  }
  return(mode_phenotype)
}

xlogx <- function(x){ #Use for computing entropy
  if (x == 0.0){
    return(0.0)
  } else {
    return(x*log(x))
  }
}

get_entropy <- function(pop){ #Sum of entropy in each element
  popsize <- length(pop)
  H <- 0.0
  for(trait in 1:indivlength){
    gcounts <- get_trait_counts(pop,trait)
    plist <- gcounts/popsize #Normalize into probabilities
    for(p in plist){
      H <- H + xlogx(p)
    }
  }
  return(-H)
}

bin_mismatch <- function(pop,env){ #Bins mismatch of population
  mm <- mapply(get_mismatch,pop,env)
  fout <- c()
  for (v in seq(0,nchar(env))){
    fout <- append(fout,length(which(mm==v)))
  }
  return(fout)
}

# Mutation, this is used for a cheat initialization
mutate <- function(indiv,mut_rate){
  indiv_copy <- indiv
  for(g in 1:indivlength){
    r <- runif(1,0,1)
    if (r<mut_rate){
      g1 <- sample(glist) #Mutated value
      substr(indiv_copy,g,g) <- g1
    }
  }
  return(indiv_copy)
}

### Randomly initialize population
init_pop <- function(popsize){ #Start from semi-fitted population
  pop <- c()
  for (i in seq(1,popsize,1)){ 
    indiv <- mutate(target0, 0.5) 
    pop <- append(pop,indiv)
  }
  print(sprintf("Initial population size: %d",length(pop)))
  
  return(pop)
}

init_pop0 <- function(popsize){ #Randomly initialized population
  pop <- c()
  for (i in seq(1,popsize,1)){ 
    indiv <- ""
    for (j in 1:indivlength){
      g <- sample(glist,1)
      indiv <- paste(indiv,g,sep="")
    }
    pop <- append(pop,indiv)
  }  
  print(sprintf("Initial population size: %d",length(pop)))
    
  return(pop)
}

init_pop1 <- function(popsize){ #Extreme: No adaptive variants at start.
  pop <- c()
  for (i in seq(1,popsize,1)){ 
    indiv <- mutate(target, 1.0) 
    pop <- append(pop,indiv)
  }
  print(sprintf("Initial population size: %d",length(pop)))
  
  return(pop)
}


### Mismatch

get_mismatch <- function(str1,str2){
  s <- 0
  for(i in 1:nchar(str1)){
    v1 <- substr(str1,i,i)
    v2 <- substr(str2,i,i)
    if(v1!=v2){
      s <- s+1
    }
  }
  return(s)
}

get_fitness <- function(s){ #Fitness given mismatch
    return(exp(-selstrength*s))
}

### Selection

select <- function(pop, fits){
  tries <- 0
  parents <- c()
  while(length(parents)<maxpopsize){
    r <- runif(1,0,1)
    indiv <- sample(c(1:maxpopsize),1) #Choose index
    rfit <- fits[indiv]
    if(r<rfit){
      parents <- append(parents,pop[indiv])
    }
    tries <- tries + 1
    if(tries>1000*maxpopsize){ #Too many tries
      break #Terminate loop
    }
  }
  return(parents)
}

### Recombination

mate <- function(mom,dad){ #Recombination between 2 individuals
  kid1 <- mom
  kid2 <- dad
  for(g in 1:indivlength){
    r <- runif(1,0,1)
    if(r<0.5){
      substr(kid1,g,g) <- substr(dad,g,g) #Feed dad's gene instead
      substr(kid2,g,g) <- substr(mom,g,g) #Feed mom's gene instead
    }
  }
  return(c(kid1,kid2))
}

### Reproduction

reproduce <- function(parents){ #Reproduce next generation of offsprings
  popsize <- length(parents)
  if(popsize%%2 == 1){
    parents <- append(parents,parents[popsize])
    popsize <- popsize + 1
  }
  kids <- c()
  for(id in seq(1,popsize,2)){
    kidpair <- mate(parents[id],parents[id+1]) #Recombination
    for (kid in kidpair){
      kids <- append(kids,mutate(kid,mut_rate0)) #Mutation
    }
  }
  return(kids)
}

### Evolution

evolve <- function(pop,maxgen){ #
  traj_mismatch <- c() #For plotting trajectory of mismatch
  traj_entropy <- c()
  for(gen in 1:maxgen){
    #Get statistics before selection
    popmismatch <- mapply(get_mismatch,pop,target) #Get mismatch has two inputs so mapply
    fitvec <- unlist(lapply(popmismatch,get_fitness))
    meanmismatch <- mean(popmismatch)
    entropy <- get_entropy(pop)

    print(sprintf("Gen: %d, Avg. Traits: %s, Target: %s", gen, get_mode_phenotype(pop), target))
    print(sprintf("Stats: Min Mismatch: %d; Mean Mismatch: %.2f, Entropy: %.2f", min(popmismatch), meanmismatch, entropy))
    traj_mismatch <- append(traj_mismatch,meanmismatch)
    traj_entropy <- append(traj_entropy,entropy)

    #Selection
    maxfit <- max(fitvec)
    rfitvec <- fitvec/maxfit #Use rescaled fitness
    parents <- select(pop,rfitvec) #Selected individuals become parents
    
    #Reproduction
    pop <- reproduce(parents) #Kids become population in next generation
  }
  plot(traj_mismatch,type="l",xlab="Gen",ylab="Mismatch",ylim=c(0,indivlength))
  plot(traj_entropy,type="l",xlab="Gen",ylab="Entropy (nats)",ylim=c(0,max_entropy))
  return(pop)
}

pop0 <- init_pop0(maxpopsize) #Initialize
pop1 <- evolve(pop0,ngen) #Evolve for some generations

#Visualize distribution of mismatch before and after evolution

mm0 <- mapply(get_mismatch,pop0,target)
mm1 <- mapply(get_mismatch,pop1,target)
hist(mm0,freq=TRUE,breaks=seq(-0.5,18.5,1),xlab="Mismatch",
     xlim=c(0,20),ylim=c(0,100),col="#D55E00",main="")
hist(mm1,freq=TRUE,breaks=seq(-0.5,18.5,1),xlab="Mismatch",
     xlim=c(0,20),ylim=c(0,100),col="#56B4E9",main="",add=TRUE)
legend("topright",legend=c("Before","After"),fill=c("#D55E00","#56B4E9"))
