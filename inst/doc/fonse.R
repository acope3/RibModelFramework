## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----eval = FALSE-------------------------------------------------------------
# library(AnaCoDa)
# genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
# genome <- initializeGenomeObject(file = genome_file)


## ----eval = FALSE-------------------------------------------------------------
# expression_file <- "expression.csv"
# genome <- initializeGenomeObject(file = genome_file,
#                                  observed.expression.file = expression_file,
#                                  match.expression.by.id = TRUE)


## ----eval = FALSE-------------------------------------------------------------
# n_genes      <- length(genome)
# sphi_init    <- 1
# num_mixtures <- 1
# gene_assignment <- rep(1, n_genes)
# 
# parameter <- initializeParameterObject(
#   genome             = genome,
#   sphi               = sphi_init,
#   num.mixtures       = num_mixtures,
#   gene.assignment    = gene_assignment,
#   model              = "FONSE",
#   init.initiation.cost = 4
# )


## ----eval = FALSE-------------------------------------------------------------
# sphi_init    <- c(1, 1)
# num_mixtures <- 2
# gene_assignment <- sample(c(1, 2), size = n_genes, replace = TRUE,
#                            prob = c(0.5, 0.5))
# 
# parameter <- initializeParameterObject(
#   genome             = genome,
#   sphi               = sphi_init,
#   num.mixtures       = num_mixtures,
#   gene.assignment    = gene_assignment,
#   model              = "FONSE",
#   mixture.definition = "allUnique",
#   init.initiation.cost = 4
# )


## ----eval = FALSE-------------------------------------------------------------
# parameter$initMutationCategories(c("mutation_cat1.csv"), numCategories = 1,
#                                   fix = FALSE)
# parameter$initSelectionCategories(c("selection_cat1.csv"), numCategories = 1,
#                                    fix = FALSE)


## ----eval = FALSE-------------------------------------------------------------
# model <- initializeModelObject(parameter, "FONSE")
# 
# mcmc <- initializeMCMCObject(
#   samples        = 5000,
#   thinning       = 10,
#   adaptive.width = 50,
#   est.expression = TRUE,
#   est.csp        = TRUE,
#   est.hyper      = TRUE
# )


## ----eval = FALSE-------------------------------------------------------------
# runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 4)


## ----eval = FALSE-------------------------------------------------------------
# plot(mcmc)


## ----eval = FALSE-------------------------------------------------------------
# trace <- parameter$getTraceObject()
# plot(trace, what = "InitiationCost")


## ----eval = FALSE-------------------------------------------------------------
# # Mutation bias traces for all amino acids
# plotCodonSpecificParameters(parameter, what = "Mutation")
# 
# # Nonsense error coefficient traces
# plotCodonSpecificParameters(parameter, what = "Selection")


## ----eval = FALSE-------------------------------------------------------------
# trace    <- parameter$getTraceObject()
# a1_trace <- trace$getInitiationCostTrace()
# 
# # Discard burn-in (first 50% of samples)
# n_samples <- length(a1_trace)
# a1_post   <- a1_trace[ceiling(n_samples / 2):n_samples]
# 
# cat("a1 posterior mean:", mean(a1_post), "\n")
# cat("a1 95% CI:        ", quantile(a1_post, c(0.025, 0.975)), "\n")


## ----eval = FALSE-------------------------------------------------------------
# csp <- getCSPEstimates(parameter, samples = 1000, mixture = 1)
# 
# # Mutation bias (DeltaM)
# head(csp$Mutation)
# 
# # Nonsense error coefficients (DeltaOmega)
# head(csp$Selection)


## ----eval = FALSE-------------------------------------------------------------
# # Posterior mean of DeltaOmega for codon "AAA" in mixture 1
# # paramType = 1 for selection (DeltaOmega)
# omega_aaa <- parameter$getCodonSpecificPosteriorMean(
#   mixtureElement    = 1,
#   samples           = 1000,
#   codon             = "AAA",
#   paramType         = 1,
#   withoutReference  = TRUE
# )


## ----eval = FALSE-------------------------------------------------------------
# phi_estimates <- getExpressionEstimates(parameter, genome, samples = 1000)
# head(phi_estimates)


## ----eval = FALSE-------------------------------------------------------------
# # Fix mutation bias at current values
# parameter$fixDM()
# 
# # Fix nonsense error coefficients at current values
# parameter$fixDOmega()
# 
# # Fix the initiation cost at its current value
# parameter$fixedInitiationCost()


## ----eval = FALSE-------------------------------------------------------------
# parameter$initMutationCategories(c("mutation_cat1.csv"), numCategories = 1,
#                                   fix = TRUE)


## ----eval = FALSE-------------------------------------------------------------
# # Enable restart file output every 100 MCMC samples
# setRestartSettings(mcmc, filename = "fonse_restart", samples = 100,
#                     overwrite = FALSE)
# 
# # Resume from a restart file
# parameter_resumed <- initializeParameterObject(
#   init.with.restart.file = "fonse_restart.rst"
# )
# model_resumed <- initializeModelObject(parameter_resumed, "FONSE")
# mcmc_resumed  <- initializeMCMCObject(samples = 5000, thinning = 10,
#                                        adaptive.width = 50,
#                                        est.expression = TRUE,
#                                        est.csp        = TRUE,
#                                        est.hyper      = TRUE)
# runMCMC(mcmc = mcmc_resumed, genome = genome, model = model_resumed, ncores = 4)


## ----eval = FALSE-------------------------------------------------------------
# writeMCMCObject(mcmc, file = "fonse_mcmc.Rda")
# mcmc_loaded <- loadMCMCObject(file = "fonse_mcmc.Rda")

