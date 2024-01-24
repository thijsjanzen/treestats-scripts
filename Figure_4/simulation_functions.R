require(TreeSim)
require(DDD)
require(PBD)
require(diversitree)
require(bissesim)
require(physim)

spec_rate <- 0.5
extinct_rate <- 0.1

simulate_bd_tree <- function(max_t, num_lin, use_extinct) {

  effective_extinction <- use_extinct * extinct_rate

  focal_tree <- TreeSim::sim.bd.taxa.age(n = num_lin,
                                         numbsim = 1,
                                         lambda = spec_rate +
                                                     effective_extinction,
                                         mu = effective_extinction,
                                         age = max_t,
                                         mrca = TRUE)
  return(focal_tree[[1]])
}

simulate_ddd_tree <- function(max_t, num_lin, use_extinct) {

  effective_extinction <- use_extinct * extinct_rate
  effective_diversification <- spec_rate + effective_extinction

  pars <- c(effective_diversification,
            effective_extinction,
            num_lin * 1.2)

  focal_tree <- physim::sim_ddd(lambda = pars[1],
                                mu = pars[2],
                                K = pars[3],
                                max_t = max_t,
                                num_species = num_lin)

  focal_taxa <- length(focal_tree$tip.label)
  while (focal_taxa != num_lin) {
    focal_tree <- physim::sim_ddd(lambda = pars[1],
                                  mu = pars[2],
                                  K = pars[3],
                                  max_t = max_t,
                                  num_species = num_lin)
  }
  return(focal_tree)
}

simulate_pbd_tree <- function(max_t, num_lin, use_extinct) {
  effective_extinction <- use_extinct * extinct_rate
  focal_birth1 <- spec_rate + effective_extinction
  focal_lambda <- 10
  focal_birth2 <- 0.0
  focal_death <-  effective_extinction
  focal_death2 <- effective_extinction

  focal_tree <- physim::pbd_sim_new(pars = c(focal_birth1,
                                      focal_lambda,
                                      focal_birth2,
                                      focal_death,
                                      focal_death2),
                             age = max_t)

  focal_taxa <- length(focal_tree$tip.label)
  while (focal_taxa != num_lin) {
    focal_tree <- physim::pbd_sim_new(pars = c(focal_birth1,
                                               focal_lambda,
                                               focal_birth2,
                                               focal_death,
                                               focal_death2),
                                      age = max_t)
    focal_taxa <- length(focal_tree$tip.label)
  }
  return(focal_tree)
}

simulate_sse_tree <- function(max_t, num_lin, use_extinct) {
  effective_extinction <- use_extinct * extinct_rate

  factor <- 1.5

  lambda0 <- spec_rate + effective_extinction
  lambda1 <- (spec_rate + effective_extinction) * factor
  mu0 <- effective_extinction
  mu1 <- effective_extinction * factor
  q01 <- 0.1
  q10 <- 0.1

  pars <- c(lambda0, lambda1, mu0, mu1, q01, q10)

  focal_tree <- bissesim::bisse_sim(pars = pars,
                                    crown_age = max_t,
                                    num_species = num_lin,
                                    init_state = 0)
  focal_taxa <- length(focal_tree$phy$tip.label)
  while (focal_taxa != num_lin) {
    focal_tree <- bissesim::bisse_sim(pars = pars,
                                      crown_age = max_t,
                                      num_species = num_lin,
                                      init_state = 0)
    focal_taxa <- length(focal_tree$phy$tip.label)
  }

  # remove superfluous information to reduce file size
  return(focal_tree$phy)
}

simulate_tree <- function(model, max_t, num_lin, use_extinct) {
  if (model == "BD") {
    return(simulate_bd_tree(max_t = max_t,
                            num_lin = num_lin,
                            use_extinct = use_extinct))
  }
  if (model == "DDD") {
    return(simulate_ddd_tree(max_t = max_t,
                             num_lin = num_lin,
                             use_extinct = use_extinct))
  }
  if (model == "PBD") {
    return(simulate_pbd_tree(max_t = max_t,
                             num_lin = num_lin,
                             use_extinct = use_extinct))
  }
  if (model == "SSE") {
    return(simulate_sse_tree(max_t = max_t,
                             num_lin = num_lin,
                             use_extinct = use_extinct))
  }
  return(NA)
}


