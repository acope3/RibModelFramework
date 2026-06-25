library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Regression tests for multi-file load behavior in loadParameterObject
# and loadMCMCObject.
#
# Addresses #424: combineTwoDimensionalTrace was called positionally
# (3rd arg goes to `start` not `end`) and combineThreeDimensionalTrace
# modified its argument via R's local rebind without returning, so
# loadParameterObject(c(f1, f2, ...)) silently dropped all but the
# first file's CSP / sphi / synth-rate trace data while bumping
# lastIteration to the cumulative count -- producing wrong posteriors
# and the std::vector::max_size() crash from calculateQuantile.
#
# Addresses #388: loadMCMCObject always sliced curLogPostTrace[2:max],
# even for the first file, so the initial-evaluation slot at index 1
# was silently dropped on every load.
#
# Addresses part of #387: provides save/load roundtrip regression
# tests for ROC parameter objects (single-file + multi-file) and
# MCMC objects.
# ======================================================================

genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome <- initializeGenomeObject(file = genome_file)
gene_assignment <- rep(1, length(genome))


# Build a parameter object with traces populated by a tiny ROC MCMC.
buildRunOnce <- function(parameter, mcmc, model) {
    sink(tempfile())  # swallow MCMC log output
    runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()
}

# Build two parameter Rdata files by running MCMC twice on the same
# persistent parameter object.  The C++-side trace is overwritten by
# each runMCMC call, so file2 captures only the second run's trace --
# the realistic chunked-fit scenario.
buildTwoSavedParameters <- function() {
    set.seed(43)
    parameter <- initializeParameterObject(genome = genome, sphi = 1,
                                            num.mixtures = 1,
                                            gene.assignment = gene_assignment)
    mcmc <- initializeMCMCObject(samples = 3, thinning = 2, adaptive.width = 2,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

    buildRunOnce(parameter, mcmc, model)
    f1 <- tempfile(fileext = ".Rdata")
    writeParameterObject(parameter, f1)

    buildRunOnce(parameter, mcmc, model)
    f2 <- tempfile(fileext = ".Rdata")
    writeParameterObject(parameter, f2)

    list(f1 = f1, f2 = f2)
}


test_that("loadParameterObject(c(f1, f2)) concatenates sphi trace (#424)", {
    files <- buildTwoSavedParameters()
    on.exit(unlink(c(files$f1, files$f2)), add = TRUE)

    p1  <- loadParameterObject(files$f1)
    p2  <- loadParameterObject(files$f2)
    p12 <- loadParameterObject(c(files$f1, files$f2))

    sphi1  <- as.numeric(p1$getTraceObject()$getStdDevSynthesisRateTraces()[[1]])
    sphi2  <- as.numeric(p2$getTraceObject()$getStdDevSynthesisRateTraces()[[1]])
    sphi12 <- as.numeric(p12$getTraceObject()$getStdDevSynthesisRateTraces()[[1]])

    # Concatenated length = sum, minus 1 for the overlap-skip slot.
    expect_equal(length(sphi12), length(sphi1) + length(sphi2) - 1)
    # First value comes from file1; last from file2.
    expect_equal(sphi12[1],         sphi1[1])
    expect_equal(tail(sphi12, 1),   tail(sphi2, 1))
})


test_that("loadParameterObject(c(f1, f2)) concatenates CSP traces (#424)", {
    files <- buildTwoSavedParameters()
    on.exit(unlink(c(files$f1, files$f2)), add = TRUE)

    p1  <- loadParameterObject(files$f1)
    p2  <- loadParameterObject(files$f2)
    p12 <- loadParameterObject(c(files$f1, files$f2))

    # Mutation trace (paramType = 0), category 1, codon index 2.
    mut1  <- p1$getTraceObject()$getCodonSpecificParameterTrace(0)[[1]][[2]]
    mut2  <- p2$getTraceObject()$getCodonSpecificParameterTrace(0)[[1]][[2]]
    mut12 <- p12$getTraceObject()$getCodonSpecificParameterTrace(0)[[1]][[2]]

    expect_equal(length(mut12), length(mut1) + length(mut2) - 1)
    expect_equal(mut12[1],         mut1[1])
    expect_equal(tail(mut12, 1),   tail(mut2, 1))

    # Selection trace (paramType = 1), same pattern.
    sel12 <- p12$getTraceObject()$getCodonSpecificParameterTrace(1)[[1]][[2]]
    sel1  <- p1$getTraceObject()$getCodonSpecificParameterTrace(1)[[1]][[2]]
    sel2  <- p2$getTraceObject()$getCodonSpecificParameterTrace(1)[[1]][[2]]
    expect_equal(length(sel12), length(sel1) + length(sel2) - 1)
    expect_equal(sel12[1],         sel1[1])
    expect_equal(tail(sel12, 1),   tail(sel2, 1))
})


test_that("loadParameterObject(c(f1, f2)) sets lastIteration = sum (#424)", {
    files <- buildTwoSavedParameters()
    on.exit(unlink(c(files$f1, files$f2)), add = TRUE)

    p1  <- loadParameterObject(files$f1)
    p2  <- loadParameterObject(files$f2)
    p12 <- loadParameterObject(c(files$f1, files$f2))

    expect_equal(p12$getLastIteration(),
                 p1$getLastIteration() + p2$getLastIteration())
})


test_that("loadMCMCObject preserves the first element on single-file load (#388)", {
    # Construct an MCMC save in the writeMCMCObject convention:
    # logPostTrace length = samples + 1, where index 1 is the initial 0.0
    # evaluation slot.
    logPostTrace  <- c(0, -100, -90, -85)
    logLikeTrace  <- c(0, -200, -180, -170)
    samples       <- 3L
    thinning      <- 1L
    adaptiveWidth <- 10L
    tf <- tempfile(fileext = ".Rda")
    on.exit(unlink(tf), add = TRUE)
    save(logPostTrace, logLikeTrace, samples, thinning, adaptiveWidth,
         file = tf)

    m  <- loadMCMCObject(tf)
    tr <- m$getLogPosteriorTrace()
    expect_equal(length(tr), length(logPostTrace))
    expect_equal(tr[1],         logPostTrace[1])
    expect_equal(tail(tr, 1),   tail(logPostTrace, 1))
})


test_that("loadMCMCObject(c(f1, f2)) concatenates with one overlap-skip (#388)", {
    logPostTrace  <- c(0, -100, -90, -85, -80, -78)
    logLikeTrace  <- c(0, -200, -180, -170, -160, -155)
    samples       <- 5L
    thinning      <- 1L
    adaptiveWidth <- 10L
    tf1 <- tempfile(fileext = ".Rda")
    save(logPostTrace, logLikeTrace, samples, thinning, adaptiveWidth,
         file = tf1)

    logPostTrace  <- c(-78, -75, -72, -70)   # first elem duplicates tf1's last
    logLikeTrace  <- c(-155, -150, -147, -145)
    samples       <- 3L
    tf2 <- tempfile(fileext = ".Rda")
    save(logPostTrace, logLikeTrace, samples, thinning, adaptiveWidth,
         file = tf2)
    on.exit(unlink(c(tf1, tf2)), add = TRUE)

    m  <- loadMCMCObject(c(tf1, tf2))
    tr <- m$getLogPosteriorTrace()
    # 6 + 4 - 1 overlap-skip = 9
    expect_equal(length(tr),       9L)
    expect_equal(tr[1],            0)
    expect_equal(tail(tr, 1),      -70)
    # Spot-check that there's no duplicated overlap value.
    expect_equal(tr[6], -78)
    expect_equal(tr[7], -75)
})
