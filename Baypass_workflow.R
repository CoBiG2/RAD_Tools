################################################################################
#Strat by running this:
#g_baypass -npop 16 -gfile Qsuber_GBS_noBul.baypass -outprefix core/Qsuber_core -nthreads 6 > core/Qsuber_core.log
###############################################################################


require(corrplot) ; require(ape)
#source the baypass R functions (check PATH)
source("~/Software/Science/baypass_2.1/utils/baypass_utils.R")

# Define some variables
popname_file = "~/Desktop/Qsuber_GBS/Clust5/02-filtered_vcfs/popnames_noBul.txt"
geno_file = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/Qsuber_GBS_noBul.baypass"
num_SNPs = 1067
core_omega_file = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/core/Qsuber_core_mat_omega.out"
core_pi_xtx_file = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/core/Qsuber_core_summary_pi_xtx.out"
core_summary_beta_params = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/core/Qsuber_core_summary_beta_params.out"
pod_mat_omega = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/core/Qsuber_core_POD_mat_omega.out"
pod_summary_beta_params = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/core/Qsuber_core_POD_summary_beta_params.out"
pod_summary_pi_xtx = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/core/Qsuber_core_POD_summary_pi_xtx.out"
covis_summary_betai_reg = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/mcmc_core/Qsuber_mcmc_core_summary_betai_reg.out"
covis2_summary_betai_reg = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/mcmc_core/Qsuber_mcmc_core2_summary_betai_reg.out"
covmcmc_summary_betai = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/mcmc_std/Qsuber_mcmc_std_summary_betai.out"
covmcmc_summary_pi_xtx = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/mcmc_std/Qsuber_mcmc_std_summary_pi_xtx.out"
covaux_summary_betai = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/mcmc_aux/Qsuber_mcmc_aux_summary_betai.out"
covaux_summary_pi_xtx = "~/Desktop/Qsuber_GBS/Clust5/04-CenterSNP/08-Baypass/mcmc_aux/Qsuber_mcmc_aux_summary_pi_xtx.out"




#upload estimate of omega
omega=as.matrix(read.table(core_omega_file))
pop.names = c(as.matrix(read.table(popname_file)))

dimnames(omega)=list(pop.names,pop.names)
#Compute and visualize the correlation matrix
cor.mat=cov2cor(omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

#Visualize the correlation matrix as hierarchical clustering tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

#Estimates of the XtX differentiation measures
anacore.snp.res=read.table(core_pi_xtx_file,h=T)
plot(anacore.snp.res$M_XtX)


#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution
pi.beta.coef=read.table(core_summary_beta_params,h=T)$Mean
#upload the original data to obtain total allele count
current.data<-geno2YN(geno_file)
#Create the POD
simu.bta<-simulate.baypass(omega.mat=omega,nsnp=num_SNPs,
                           sample.size=current.data$NN,
                           beta.pi=pi.beta.coef,pi.maf=0,suffix="btapods")

################################################################################
#Run this:
#g_baypass -npop 16 -gfile core/G.btapods -outprefix core/Qsuber_core_POD -nthreads 6 > core/Qsuber_core_POD.log
###############################################################################

#############################################################
# Remember to move G.betapods from ~/ to somewhere sensible!#
#############################################################

#######################################################
#Sanity Check: Compare POD and original data estimates
#######################################################
#get estimate of omega from the POD analysis
pod.omega=as.matrix(read.table(pod_mat_omega))
plot(pod.omega,omega) ; abline(a=0,b=1)
fmd.dist(pod.omega,omega)

#get estimates (post. mean) of both the a_pi and b_pi parameters of
#the Pi Beta distribution from the POD analysis
pod.pi.beta.coef=read.table(pod_summary_beta_params,h=T)$Mean
plot(pod.pi.beta.coef,pi.beta.coef) ; abline(a=0,b=1)

#######################################################
#XtX calibration
#######################################################
#get the pod XtX
pod.xtx=read.table(pod_summary_pi_xtx,h=T)$M_XtX
#compute the 1% threshold
pod.thresh=quantile(pod.xtx,probs=0.99)
#add the thresh to the actual XtX plot
plot(anacore.snp.res$M_XtX)
abline(h=pod.thresh,lty=2)

################################################################################
#Run this:
#g_baypass -npop 16 -gfile Qsuber_GBS_noBul.baypass -outprefix mcmc_core/Qsuber_mcmc_core -nthreads 6 -efile ENVFILE  > mcmc_core/Qsuber_mcmc_core.log
###############################################################################

covis.snp.res=read.table(covis_summary_betai_reg,h=T)
graphics.off()
layout(matrix(1:3,3,1))
plot(covis.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))


################################################################################
#Run this:
#g_baypass -npop 16 -gfile Qsuber_GBS_noBul.baypass -outprefix mcmc_core/Qsuber_mcmc_core2 -nthreads 6 -efile ENVFILE -omegafile core/Qsuber_core_mat_omega.out > mcmc_core/Qsuber_mcmc_core2.log
###############################################################################

covis2.snp.res=read.table(covis2_summary_betai_reg,h=T)
graphics.off()
layout(matrix(1:3,3,1))
plot(covis2.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis2.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis2.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))

################################################################################
#Run this:
#g_baypass -npop 16 -gfile Qsuber_GBS_noBul.baypass -outprefix mcmc_std/Qsuber_mcmc_std -nthreads 6 -efile ENVFILE -omegafile core/Qsuber_core_mat_omega.out -covmcmc > mcmc_std/Qsuber_mcmc_std.log
###############################################################################

covmcmc.snp.res=read.table(covmcmc_summary_betai,h=T)
covmcmc.snp.xtx=read.table(covmcmc_summary_pi_xtx,h=T)$M_XtX
graphics.off()
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res$eBPmc,xlab="SNP",ylab="eBPmc")
plot(covmcmc.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covmcmc.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")


################################################################################
#Run this:
#g_baypass -npop 16 -gfile Qsuber_GBS_noBul.baypass -outprefix mcmc_aux/Qsuber_mcmc_aux -nthreads 6 -efile ENVFILE -omegafile core/Qsuber_core_mat_omega.out -auxmodel > mcmc_aux/Qsuber_mcmc_aux.log
###############################################################################

covaux.snp.res=read.table(covaux_summary_betai,h=T)
covaux.snp.xtx=read.table(covaux_summary_pi_xtx,h=T)$M_XtX
graphics.off()
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB.,xlab="SNP",ylab="BFmc (in dB)")
plot(covaux.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")
