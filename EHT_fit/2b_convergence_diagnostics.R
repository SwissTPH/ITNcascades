rm(list=ls())
library(dplyr)
library(ggplot2)
library(readxl)
library(lubridate)
library(stringr)
library(cowplot)
library(tidyr)
library(bayesplot)
library(rstan)


mainDir="."
scriptDir=file.path(".","EHT_fit")

stanDir=file.path(mainDir, "results/stan_outputs")
plotDir=file.path(mainDir, "plots")
dataDir=file.path(scriptDir, "processed_data")
outputsDir=file.path(mainDir, "results/csv_outputs")

niter=6000
nwarmup=3000
nchains=3
bayesplot::color_scheme_set("teal")

######################################################
# KIBONDO 2022
######################################################

results_kibondo_w<- readRDS(file.path(stanDir,"/stan_kibondo_washed_72.rds"))
results_kibondo_unw<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_72.rds"))

kibondo_names=c("xi[0]", "xi[Pyrethroid]",  "xi[IG2]",
               "kappa[0]", "kappa[Pyrethroid]", "kappa[IG2]",
               "pi[0]","pi[Pyrethroid]",  "pi[IG2]", 
               "Phi[0]", "Phi[Pyrethroid]",  "Phi[IG2]",
               "alpha_0","mu_0","pc_0", "lp__")
kibondo_names_plot=c("pi[Pyrethroid]","kappa[Pyrethroid]","xi[Pyrethroid]",
                     "pi[IG2]","kappa[IG2]","xi[IG2]")

names(results_kibondo_unw)=kibondo_names
names(results_kibondo_w)=kibondo_names

posterior_kibondo_unw<- rstan::extract(results_kibondo_unw,
                                      pars= kibondo_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_kibondo_unw=bayesplot::mcmc_trace(posterior_kibondo_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 2, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_kibondo_w<- rstan::extract(results_kibondo_w,
                                    pars= kibondo_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_kibondo_w=bayesplot::mcmc_trace(posterior_kibondo_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 2, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_kibondo_unw,trace_kibondo_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_kibondo.png"), width = 12, height = 12, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_kibondo_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_kibondo_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# BIT055
######################################################


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))

BIT055_names=c("xi[0]", "xi[Pyrethroid]",  "xi[OlysetPlus]",
               "kappa[0]", "kappa[Pyrethroid]", "kappa[OlysetPlus]",
               "pi[0]","pi[Pyrethroid]",  "pi[OlysetPlus]", 
               "Phi[0]", "Phi[Pyrethroid]",  "Phi[OlysetPlus]",
               "alpha_0","mu_0","pc_0", "lp__")
BIT055_names_plot=c("pi[Pyrethroid]","kappa[Pyrethroid]","xi[Pyrethroid]",
                    "pi[OlysetPlus]","kappa[OlysetPlus]","xi[OlysetPlus]")

names(results_bit055_unw)=BIT055_names
names(results_bit055_w)=BIT055_names

posterior_bit055_unw<- rstan::extract(results_bit055_unw,
                                      pars= BIT055_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit055_unw=bayesplot::mcmc_trace(posterior_bit055_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 2, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit055_w<- rstan::extract(results_bit055_w,
                                    pars= BIT055_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_bit055_w=bayesplot::mcmc_trace(posterior_bit055_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 2, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit055_unw,trace_bit055_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit055.png"), width = 12, height = 12, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_bit055_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit055_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))

######################################################
# Assenga, Tanzania
######################################################


results_bit103_unw<- readRDS(file.path(stanDir,"/stan_assenga_unwashed_72_all.rds"))
results_bit103_w<- readRDS(file.path(stanDir,"/stan_assenga_washed_72_all.rds"))

BIT103_names=c("xi[0]", "xi[Pyrethroid]", "xi[OlysetPlus]", "xi[IG2]",
               "kappa[0]","kappa[Pyrethroid]", "kappa[OlysetPlus]", "kappa[IG2]",
               "pi[0]", "pi[Pyrethroid]","pi[OlysetPlus]", "pi[IG2]", 
               "Phi[0]","Phi[Pyrethroid]", "Phi[OlysetPlus]", "Phi[IG2]",
               "alpha_0","mu_0", "pc_0", "lp__")
BIT103_names_plot=c("pi[Pyrethroid]","kappa[Pyrethroid]","xi[Pyrethroid]",
                    "pi[OlysetPlus]","kappa[OlysetPlus]","xi[OlysetPlus]",
                    "pi[IG2]","kappa[IG2]","xi[IG2]")

names(results_bit103_unw)=BIT103_names
names(results_bit103_w)=BIT103_names

posterior_bit103_unw<- rstan::extract(results_bit103_unw,
                                  pars= BIT103_names_plot,
                                  inc_warmup = TRUE,permuted = FALSE)



trace_bit103_unw=bayesplot::mcmc_trace(posterior_bit103_unw,
                      n_warmup = nwarmup,
                      facet_args = list(nrow = 3, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit103_w<- rstan::extract(results_bit103_w,
                                      pars= BIT103_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit103_w=bayesplot::mcmc_trace(posterior_bit103_w,
                      n_warmup = 7000,
                      facet_args = list(nrow = 3, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit103_unw,trace_bit103_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit103.png"), width = 12, height = 18, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_bit103_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit103_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# Assenga, C6ote d'Ivoire
######################################################


results_bit103_unw_72_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_unwashed_72.rds"))
results_bit103_w_72_CI<- readRDS(file.path(stanDir,"/stan_assengaCotedIvoire_washed_72.rds"))

names(results_bit103_unw_72_CI)=BIT103_names
names(results_bit103_w_72_CI)=BIT103_names

posterior_bit103_unw_CI<- rstan::extract(results_bit103_unw_72_CI,
                                      pars= BIT103_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit103_unw_CI=bayesplot::mcmc_trace(posterior_bit103_unw_CI,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 3, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit103_w_CI<- rstan::extract(results_bit103_w_72_CI,
                                    pars= BIT103_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_bit103_w_CI=bayesplot::mcmc_trace(posterior_bit103_w_CI,
                                     n_warmup = 7000,
                                     facet_args = list(nrow = 3, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit103_unw_CI,trace_bit103_w_CI, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit103_CI.png"), width = 12, height = 18, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_bit103_unw_72_CI, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit103_w_72_CI, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# BIT059
######################################################


results_bit059_w<- readRDS(file.path(stanDir,"/stan_odufuwa_washed_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_odufuwa_unwashed_24.rds"))

BIT059_names=BIT055_names
BIT059_names_plot=BIT055_names_plot

names(results_bit059_unw)=BIT059_names
names(results_bit059_w)=BIT059_names

posterior_bit059_unw<- rstan::extract(results_bit059_unw,
                                      pars= BIT059_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit059_unw=bayesplot::mcmc_trace(posterior_bit059_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 2, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit059_w<- rstan::extract(results_bit059_w,
                                    pars= BIT059_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_bit059_w=bayesplot::mcmc_trace(posterior_bit059_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 2, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit059_unw,trace_bit059_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_odufuwa.png"), width = 12, height = 12, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_bit059_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit059_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# BIT080
######################################################

results_bit080_w<- readRDS(file.path(stanDir,"/stan_bit080_washed_72.rds"))
results_bit080_unw<- readRDS(file.path(stanDir,"/stan_bit080_unwashed_72.rds"))

bit080_names=kibondo_names
bit080_names_plot=kibondo_names_plot

names(results_bit080_unw)=bit080_names
names(results_bit080_w)=bit080_names

posterior_bit080_unw<- rstan::extract(results_bit080_unw,
                                       pars= bit080_names_plot,
                                       inc_warmup = TRUE,permuted = FALSE)



trace_bit080_unw=bayesplot::mcmc_trace(posterior_bit080_unw,
                                        n_warmup = nwarmup,
                                        facet_args = list(nrow = 2, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit080_w<- rstan::extract(results_bit080_w,
                                     pars= bit080_names_plot,
                                     inc_warmup = TRUE,permuted = FALSE)



trace_bit080_w=bayesplot::mcmc_trace(posterior_bit080_w,
                                      n_warmup = nwarmup,
                                      facet_args = list(nrow = 2, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit080_unw,trace_bit080_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit080.png"), width = 12, height = 12, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_bit080_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit080_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))



#############
# Martin, gambiae 
##############

results_martin_w72<- readRDS(file.path(stanDir,"/stan_martin_36m_72_gambiae.rds"))
results_martin_unw72<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_gambiae.rds"))

martin_names=BIT103_names
martin_names_plot=BIT103_names_plot

names(results_martin_unw72)=martin_names
names(results_martin_w72)=martin_names

posterior_martin_unw<- rstan::extract(results_martin_unw72,
                                      pars= martin_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_martin_unw=bayesplot::mcmc_trace(posterior_martin_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 3, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_martin_w<- rstan::extract(results_martin_w72,
                                    pars= martin_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_martin_w=bayesplot::mcmc_trace(posterior_martin_w,
                                     n_warmup = 7000,
                                     facet_args = list(nrow = 3, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_martin_unw,trace_martin_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_martin_gambiae.png"), width = 12, height = 18, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_martin_unw72, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_martin_w72, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


#############
# Martin , funestus
##############

results_martinf_w72<- readRDS(file.path(stanDir,"/stan_martin_36m_72_funestus.rds"))
results_martinf_unw72<- readRDS(file.path(stanDir,"/stan_martin_unwashed_72_funestus.rds"))

names(results_martinf_unw72)=martin_names
names(results_martinf_w72)=martin_names

posterior_martinf_unw<- rstan::extract(results_martinf_unw72,
                                      pars= martin_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_martinf_unw=bayesplot::mcmc_trace(posterior_martinf_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 3, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_martinf_w<- rstan::extract(results_martinf_w72,
                                    pars= martin_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_martinf_w=bayesplot::mcmc_trace(posterior_martinf_w,
                                     n_warmup = 7000,
                                     facet_args = list(nrow = 3, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_martinf_unw,trace_martinf_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_martin_funestus.png"), width = 12, height = 18, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_martinf_unw72, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_martinf_w72, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# Nguessan 2016
######################################################

results_nguessan_w<- readRDS(file.path(stanDir,"/stan_nguessan_washed.rds"))
results_nguessan_unw<- readRDS(file.path(stanDir,"/stan_nguessan_unwashed.rds"))

nguessan_names=kibondo_names
nguessan_names_plot=kibondo_names_plot

names(results_nguessan_unw)=nguessan_names
names(results_nguessan_w)=nguessan_names

posterior_nguessan_unw<- rstan::extract(results_bit080_unw,
                                      pars= nguessan_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_nguessan_unw=bayesplot::mcmc_trace(posterior_nguessan_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 2, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_nguessan_w<- rstan::extract(results_nguessan_w,
                                    pars= nguessan_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_nguessan_w=bayesplot::mcmc_trace(posterior_nguessan_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 2, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_nguessan_unw,trace_nguessan_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_nguessan.png"), width = 12, height = 12, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_nguessan_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_nguessan_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))




######################################################
# Sovegnon 2024
######################################################

results_sovegnonIG2_unw<- readRDS(file.path(stanDir,"/stan_sovegnonIG2_unwashed.rds"))
results_sovegnonIG2_w<- readRDS(file.path(stanDir,"/stan_sovegnonIG2_aged.rds"))
results_sovegnonIG1_unw<- readRDS(file.path(stanDir,"/stan_sovegnonIG1_unwashed.rds"))
results_sovegnonIG1_w<- readRDS(file.path(stanDir,"/stan_sovegnonIG1_aged.rds"))


sovegnonIG1_names=c("xi[0]", "xi[Pyrethroid]",
                "kappa[0]", "kappa[Pyrethroid]",
                "pi[0]","pi[Pyrethroid]", 
                "Phi[0]", "Phi[Pyrethroid]",
                "alpha_0","mu_0","pc_0", "lp__")
sovegnonIG2_names=c("xi[0]", "xi[IG2]",
                    "kappa[0]",  "kappa[IG2]",
                    "pi[0]", "pi[IG2]", 
                    "Phi[0]",   "Phi[IG2]",
                    "alpha_0","mu_0","pc_0", "lp__")

sovegnonIG1_names_plot=c("pi[Pyrethroid]","kappa[Pyrethroid]","xi[Pyrethroid]")

sovegnonIG2_names_plot=c("pi[IG2]","kappa[IG2]","xi[IG2]")

names(results_sovegnonIG1_unw)=sovegnonIG1_names
names(results_sovegnonIG1_w)=sovegnonIG1_names
names(results_sovegnonIG2_unw)=sovegnonIG2_names
names(results_sovegnonIG2_w)=sovegnonIG2_names

posterior_sovegnonIG1_unw<- rstan::extract(results_sovegnonIG1_unw,
                                        pars= sovegnonIG1_names_plot,
                                        inc_warmup = TRUE,permuted = FALSE)



trace_sovegnonIG1_unw=bayesplot::mcmc_trace(posterior_sovegnonIG1_unw,
                                         n_warmup = nwarmup,
                                         facet_args = list(nrow = 1, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_sovegnonIG1_w<- rstan::extract(results_sovegnonIG1_w,
                                      pars= sovegnonIG1_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_sovegnonIG1_w=bayesplot::mcmc_trace(posterior_sovegnonIG1_w,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 1, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

posterior_sovegnonIG2_unw<- rstan::extract(results_sovegnonIG2_unw,
                                           pars= sovegnonIG2_names_plot,
                                           inc_warmup = TRUE,permuted = FALSE)



trace_sovegnonIG2_unw=bayesplot::mcmc_trace(posterior_sovegnonIG2_unw,
                                            n_warmup = nwarmup,
                                            facet_args = list(nrow = 1, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_sovegnonIG2_w<- rstan::extract(results_sovegnonIG2_w,
                                         pars= sovegnonIG2_names_plot,
                                         inc_warmup = TRUE,permuted = FALSE)



trace_sovegnonIG2_w=bayesplot::mcmc_trace(posterior_sovegnonIG2_w,
                                          n_warmup = nwarmup,
                                          facet_args = list(nrow = 1, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_sovegnonIG1_unw,trace_sovegnonIG2_unw,trace_sovegnonIG1_w,trace_sovegnonIG2_w, ncol=1, 
          labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_sovegnon.png"), width = 12, height = 12, dpi=100)

# count divergences
sum(sapply(get_sampler_params(results_sovegnonIG1_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_sovegnonIG1_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_sovegnonIG2_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_sovegnonIG2_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))

