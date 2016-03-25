#!/usr/bin/Rscript
# Copyright 2016 Francisco Pina Martins <f.pinamartins@gmail.com>
# This file is part of Baypass_workflow.R.
# Baypass_workflow.R is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Baypass_workflow.R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Baypass_workflow.R. If not, see <http://www.gnu.org/licenses/>.

require(corrplot)
require(ape)
#source the baypass R functions (check PATH)
source("~/Software/Science/baypass_2.1/utils/baypass_utils.R")

# Define some variables. This is where we define local files & paths.
baypass_executable = "~/Software/Science/baypass_2.1/sources/g_baypass"
popname_file = "~/Desktop/Qsuber_GBS/Clust5/02-filtered_vcfs/popnames_noBul.txt"
envfile = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/ENVFILE"
geno_file = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/Qsuber_GBS_noBul.baypass"
prefix = "Qsuber"
num_pops = 16
num_SNPs = 1067
num_threads = 6


# Everything below this point should be fully automated.

basepath = dirname(geno_file)
coredir = paste(basepath, "/core/", sep="")
mcmc_coredir = paste(basepath, "/mcmc_core/", sep="")
mcmc_stddir = paste(basepath, "/mcmc_std/", sep="")
mcmc_auxdir = paste(basepath, "/mcmc_aux/", sep="")

dir.create(coredir)
dir.create(mcmc_coredir)
dir.create(mcmc_stddir)
dir.create(mcmc_auxdir)

core_omega_file = paste(coredir, prefix, "_core_mat_omega.out", sep="")
core_pi_xtx_file = paste(coredir, prefix, "_core_summary_pi_xtx.out", sep="")
core_summary_beta_params = paste(coredir, prefix, "_core_summary_beta_params.out", sep="")
pod_mat_omega = paste(coredir, prefix, "_core_POD_mat_omega.out", sep="")
pod_summary_beta_params = paste(coredir, prefix, "_core_POD_summary_beta_params.out", sep="")
pod_summary_pi_xtx = paste(coredir, prefix, "_core_POD_summary_pi_xtx.out", sep="")
covis_summary_betai_reg = paste(mcmc_coredir, prefix, "_mcmc_core_summary_betai_reg.out", sep="")
covis2_summary_betai_reg = paste(mcmc_coredir, prefix, "_mcmc_core2_summary_betai_reg.out", sep="")
covmcmc_summary_betai = paste(mcmc_stddir, prefix, "_mcmc_std_summary_betai.out", sep="")
covmcmc_summary_pi_xtx = paste(mcmc_stddir, prefix, "_mcmc_std_summary_pi_xtx.out", sep="")
covaux_summary_betai = paste(mcmc_auxdir, prefix, "_mcmc_aux_summary_betai.out", sep="")
covaux_summary_pi_xtx = paste(mcmc_auxdir, prefix, "_mcmc_aux_summary_pi_xtx.out", sep="")

### Run the first command:
command1 = paste(baypass_executable, " -npop ", num_pops, " -gfile ",
                 geno_file, " -outprefix ", coredir, prefix, "_core",
                 " -nthreads ", num_threads, sep="")
system(command=command1)


#upload estimate of omega
omega=as.matrix(read.table(core_omega_file))
pop.names = c(as.matrix(read.table(popname_file)))

dimnames(omega)=list(pop.names,pop.names)
#Compute and visualize the correlation matrix
cor.mat=cov2cor(omega)
svg(filename=paste(coredir, "omega_corr.svg") )
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))
dev.off()

#Visualize the correlation matrix as hierarchical clustering tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
svg(filename=paste(coredir, "Hier_clust_tree.svg"))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

#Estimates of the XtX differentiation measures
anacore.snp.res=read.table(core_pi_xtx_file,h=T)
svg(filename=paste(coredir, "XtX_diff.svg"))
plot(anacore.snp.res$M_XtX)
dev.off()

#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution
pi.beta.coef=read.table(core_summary_beta_params,h=T)$Mean
#upload the original data to obtain total allele count
current.data<-geno2YN(geno_file)
#Create the POD
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=num_SNPs,
                           sample.size=current.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0,suffix="btapods")


file.rename("G.btapods", paste(coredir, "G.btapods", sep=""))

###
command2 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", coredir,
                 "G.btapods", " -outprefix ", coredir, prefix, "_core_POD",
                 " -nthreads ", num_threads, sep="")
system(command=command2)


#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################
#get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table(pod_mat_omega))
svg(filename=paste(coredir, "Omega_estimate_from_POD.svg"))
plot(pod.omega,omega) ; abline(a=0,b=1)
dev.off()
fmd.dist(pod.omega,omega)

#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table(pod_summary_beta_params,h=T)$Mean
svg(filename=paste(coredir, "POD_estimates.svg"))
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)
dev.off

#######################################################
#XtX calibration
#######################################################
#get the pod XtX
pod.xtx=read.table(pod_summary_pi_xtx,h=T)$M_XtX
#compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.99)
#add the thresh to the actual XtX plot
svg(filename=paste(coredir, "XtX_POD_diff.svg"))
plot(anacore.snp.res$M_XtX)
abline(h=pod.thresh,lty=2)
dev.off()

###
command3 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_coredir, prefix, "_mcmc_core",
                 " -nthreads ", num_threads, " -efile ", envfile, sep="")
system(command=command3)

covis.snp.res=read.table(covis_summary_betai_reg,h=T)
graphics.off()
svg(filename=paste(mcmc_coredir, "BFs_layout.svg"))
layout(matrix(1:3,3,1))
plot(covis.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()


###
command4 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_coredir, prefix, "_mcmc_core2",
                 " -nthreads ", num_threads, " -efile ", envfile,
                 " -omegafile ", core_omega_file, sep="")
system(command=command4)


covis2.snp.res=read.table(covis2_summary_betai_reg,h=T)
graphics.off()
svg(filename=paste(mcmc_coredir, "BFs_layout_pass2.svg"))
layout(matrix(1:3,3,1))
plot(covis2.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis2.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis2.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

###
command5 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_stddir, prefix, "_mcmc_std",
                 " -nthreads ", num_threads, " -efile ", envfile,
                 " -omegafile ", core_omega_file, " -covmcmc", sep="")
system(command=command5)

covmcmc.snp.res=read.table(covmcmc_summary_betai,h=T)
covmcmc.snp.xtx=read.table(covmcmc_summary_pi_xtx,h=T)$M_XtX
graphics.off()
svg(filename=paste(mcmc_stddir, "BFs_layout.svg"))
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res$eBPmc,xlab="SNP",ylab="eBPmc")
plot(covmcmc.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covmcmc.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")
dev.off()

###
command6 = paste(baypass_executable, " -npop ", num_pops, " -gfile ", geno_file,
                 " -outprefix ", mcmc_auxdir, prefix, "_mcmc_aux",
                 " -nthreads ", num_threads, " -efile ", envfile,
                 " -omegafile ", core_omega_file, " -auxmodel", sep="")
system(command=command6)

covaux.snp.res=read.table(covaux_summary_betai,h=T)
covaux.snp.xtx=read.table(covaux_summary_pi_xtx,h=T)$M_XtX
graphics.off()
svg(filename=paste(mcmc_auxdir, "BFs_layout.svg"))
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB.,xlab="SNP",ylab="BFmc (in dB)")
plot(covaux.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")
dev.off()