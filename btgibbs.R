## function [pi_st, a_st, stats] = btgibbs(w, a, N_Gibbs, N_burn)
###########################################################################
# function [pi_st, a_st, stats] = btgibbs(w, a, N_Gibbs, N_burn)
# Gibbs sampler for the Bradley-Terry model
# Requires the statistics toolbox
##########
# INPUTS
#
# w[i,j] is the number of times i beats j
# a is the shape parameter for the gamma prior (default a = 1)
# If a<0, then it is estimated with a vague prior
# NGibbs: Number of Gibbs iterations
# Nburn: Number of burn-in iterations
##########
# OUTPUTS
# pi_st gives the values of the normalized skills at each iteration
# stats is a structure with some summary statistics on the parameters
###########################################################################
# Reference:
# F. Caron and A. Doucet. Efficient Bayesian inference for generalized
# Bradley-Terry models. To appear in Journal of Computational and Graphical
# Statistics, 2011.
#
# February 2011
# Author: F. Caron (INRIA Bordeaux Sud-Ouest)
# Francois.Caron@inria.fr
# http://www.math.u-bordeaux1.fr/~fcaron/
###########################################################################

btgibbs <- function(w, a, NGibbs, Nburn=1){

    K <- nrow(w)
    if(a<0){
        estimate_a <- TRUE
        a <- 1
    }
    else{
        estimate_a <- FALSE
    }
    b <- K*a - 1

    lambda <- rep(1,K)
    pi_st <- matrix(0,nrow=NGibbs,ncol=K)
    a_st <- rep(0,NGibbs)

    N <- w+t(w) # Number of comparisons between i and j
    ## N <- triupper?
    ak <- a+rowSums(w)

    ## TODO: Find nonzero entries in N
    
    pi_st[1,] <- lambda/sum(lambda)
    a_st[1] <- a
    
    for(i in 1:NGibbs){

    }


    
}
