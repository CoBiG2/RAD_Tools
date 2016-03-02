#!/usr/bin/Rscript

interactive_mode <- TRUE # Change to false if running as a script on a "headless" server.

# Read SNP data:
require(SpaceMix)

counts <- read.csv("~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/Qsuber_indiv90_miss80_maf06_noLD.vcf.counts", sep=" ")
sizes <- read.csv("~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/Qsuber_indiv90_miss80_maf06_noLD.vcf.size", sep=" ")

counts_matrix <- as.matrix(counts)
sizes_matrix <- as.matrix(sizes)


# Read geographic coordinates:.
geo.coords <- read.csv("~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/coords.csv", sep=",", header=TRUE, row.names=1)


# Run SpaceMix analysis:
run.spacemix.analysis(n.fast.reps = 10,
                      fast.MCMC.ngen = 1e5,
                      fast.model.option = "source_and_target",
                      long.model.option = "source_and_target",
                      data.type = "counts",
                      sample.frequencies=NULL,
                      mean.sample.sizes=NULL,
                      counts = counts_matrix,
                      sample.sizes = sizes_matrix,
                      sample.covariance=NULL,
                      target.spatial.prior.scale=NULL,
                      source.spatial.prior.scale=NULL,
                      spatial.prior.X.coordinates = geo.coords[,1],
                      spatial.prior.Y.coordinates = geo.coords[,2],
                      round.earth = FALSE,
                      long.run.initial.parameters=NULL,
                      k = 77,
                      loci = 1356,
                      ngen = 1e6,
                      printfreq = 1e2,
                      samplefreq = 1e3,
                      mixing.diagn.freq = 50,
                      savefreq = 1e5,
                      directory="~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/run1/",
                      prefix = "Qsuber_GBS_run1")

# Post prob. trace plot

load("~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/run1/Qsuber_GBS_run1_LongRun/Qsuber_GBS_run1_space_MCMC_output1.Robj")

if (interactive_mode == FALSE) {svg(filename="MCMC_iterations.svg")}
plot(Prob,
     xlab="MCMC iterations",
     ylab="value",
     main="Posterior probability trace plot",
     type='l')
if (interactive_mode == FALSE) {dev.off()}

# Trace plots of alpha parameters of the spatial covariance function
if (interactive_mode == FALSE) {svg(filename="MCMC_nugget.svg")}
matplot(t(nugget),type='l',
        xlab="MCMC iterations",ylab="Parameter value",
        main="Trace plot of nugget parameters")
if (interactive_mode == FALSE) {dev.off()}

# Joint marginal plot of a0 and a1
#   colored by where in the MCMC these
#   parameters took their values
if (interactive_mode == FALSE) {svg(filename="joint_marginal.svg")}
plot(a0, a1, xlab="a0", ylab="a1",
     main="Joint marginal of a0 and a1", pch=20,
     col=adjustcolor(rainbow(1000, start=4/6, end=6/6), 0.3))
legend(x="bottomright", pch=19, cex=0.8,
       col=rainbow(1000, start=4/6, end=6/6)[c(1, 500, 1000)],
       legend=c("Sampled MCMC iteration 1",
                "Sampled MCMC iteration 500",
                "Sampled MCMC iteration 1000"))
if (interactive_mode == FALSE) {dev.off()}

# Acceptance rate of a0 over the course of the
#   MCMC analysis
if (interactive_mode == FALSE) {svg(filename="acceptance.svg")}
plot(accept_rates$a0_accept_rate,
     xlab="MCMC iterations", ylab="Acceptance rate",
     main="Acceptance rate of a0", type='l',
     ylim=c(0.35, 0.6))
abline(h=0.44, col="gray", lty=2)
if (interactive_mode == FALSE) {dev.off()}


# Acceptance rates of nugget parameters over the
#   course of the MCMC analysis
if (interactive_mode == FALSE) {svg(filename="acceptance_nugget.svg")}
matplot(t(accept_rates$nugget_accept_rate),
        xlab="MCMC iterations", ylab="Acceptance rate",
        main="Acceptance rates of nuggets", type='l',
        ylim=c(0.3, 0.7))
abline(h=0.44, col="gray", lty=2)
if (interactive_mode == FALSE) {dev.off()}


## Model adequacy

load("~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/run1/Qsuber_GBS_run1_LongRun/Qsuber_GBS_run1_MCN.frequencies.list.Robj")

# Now, calculate the sample covariance from the mean centered
#   and normalized sample allele frequencies.
sample.covariance <- cov(t(MCN.frequencies.list$mean.centered.normalized.sample.frequencies),
                         use="pairwise.complete.obs")

# Create a matrix that will perform a mean-centering
#   on the parametric covariance matrix
# Then, mean-center the parametric ovariance matrix.
k <- nrow(MCN.frequencies.list$mean.centered.normalized.sample.frequencies)
MC.matrix <- diag(k) - matrix(1/last.params$inv.mean.sample.sizes /
                                  (sum(1/last.params$inv.mean.sample.sizes)),
                              nrow=k, ncol=k, byrow=TRUE)

MC.parametric.covariance <- (MC.matrix) %*%
    last.params$admixed.covariance %*%
    t(MC.matrix)

# Finally, compare the sample covariance to the parametric
#   covariance.  Ideally, there will be a very tight correspondence
#   between the data and the model.  If there is not, it may
#   be an indication either that the MCMC has not converged on
#   the stationary distribution or that the process that generated
#   the data is only poorly approximated by SpaceMix's model.

# The sample and parametric covariances can be plotted
#   against each other (if model fit is good they should
#   fall on the x=y red line)
index.matrix <- upper.tri(sample.covariance, diag=TRUE)
if (interactive_mode == FALSE) {svg(filename="model_adequacy.svg")}
plot(sample.covariance[index.matrix],
     MC.parametric.covariance[index.matrix],
     col=adjustcolor("black", 0.3), pch=20,
     xlab="Sample covariance",
     ylab="Parametric covariance",
     main="Model adequacy:\n matrix comparison")
abline(0, 1, col="red")
if (interactive_mode == FALSE) {dev.off()}


# Or the patterns of decay of covariance with
#   geographic distance can be compared between
#   the data and the model.
if (interactive_mode == FALSE) {svg(filename="model_adequacy_IBD.svg")}
plot(last.params$D[1:k,1:k][index.matrix],
     sample.covariance[index.matrix],
     pch=19, col="black",
     xlab="geogenetic distance",
     ylab="covariance",
     main="Model adequacy:\n IBD patterns")
points(last.params$D[1:k,1:k][index.matrix],
       MC.parametric.covariance[index.matrix], col="red", pch=20)
legend(x="topright", pch=19, col=c(1, 2),
       legend=c("observed", "model estimate"))
if (interactive_mode == FALSE) {dev.off()}



## Now, for the main course!

# Sample names
sample.names <- as.vector(rownames(counts))
N = length(sample.names)

# Make geo.coords a matrix
geo.coords <- as.matrix(geo.coords)

# And generate a vector of sample colors
sample.colors <- rainbow(n=N,start=1/6,end=6/6)[as.numeric(cut(geo.coords[,1],N))]

# And now generate a sample map list using a 95%
#   credible interval on parameter estimates without
#   `burning' (i.e., discarding) any sampled iterations
#   of the MCMC.
run.spacemix.map.list <- make.spacemix.map.list(MCMC.output.file = "~/Desktop/Qsuber_GBS/Clust5/11-SpaceMix/run1/Qsuber_GBS_run1_LongRun/Qsuber_GBS_run1_space_MCMC_output1.Robj",
                                                    geographic.locations = geo.coords,
                                                    name.vector = sample.names,
                                                    color.vector = sample.colors,
                                                    quantile=0.95,
                                                    burnin=0)

# Now we generate a map of the output showing sample names
#   at the locations of the maximum a posteriori (MAP)
#   geogenetic location parameter estimates
if (interactive_mode == FALSE) {svg(filename="geogenetic_map_names.svg")}
make.spacemix.map(spacemix.map.list = run.spacemix.map.list,
                  text=TRUE,
                  ellipses=FALSE,
                  source.option=FALSE)
if (interactive_mode == FALSE) {dev.off()}

if (interactive_mode == FALSE) {svg(filename="geogenetic_map_elipses.svg")}
make.spacemix.map(spacemix.map.list = run.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=FALSE)
if (interactive_mode == FALSE) {dev.off()}

if (interactive_mode == FALSE) {svg(filename="geogenetic_map_admix.svg")}
make.spacemix.map(spacemix.map.list = run.spacemix.map.list,
                  text=TRUE,
                  ellipses=TRUE,
                  source.option=TRUE)
if (interactive_mode == FALSE) {dev.off()}
