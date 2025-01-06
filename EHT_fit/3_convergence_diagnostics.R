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
stanDir=file.path(mainDir, "stan_outputs")
plotDir=file.path(mainDir, "plots")
dataDir=file.path(mainDir, "data")
outputsDir=file.path(mainDir, "csv_outputs")
scriptDir=file.path(mainDir,"EHT_fit")

niter=6000
nwarmup=3000
nchains=3
bayesplot::color_scheme_set("teal")

######################################################
# KIBONDO 2022
######################################################

results_kibondo_w<- readRDS(file.path(stanDir,"/stan_kibondo_washed_controlUnw_72.rds"))
results_kibondo_unw<- readRDS(file.path(stanDir,"/stan_kibondo_unwashed_72.rds"))

kibondo_names=c("xi[0]",  "xi[IG2]",
               "kappa[0]", "kappa[IG2]",
               "pi[0]",  "pi[IG2]", 
               "Phi[0]",  "Phi[IG2]",
               "alpha_0","mu_0", "lp__")
kibondo_names_plot=c("pi[IG2]","kappa[IG2]","xi[IG2]")

names(results_kibondo_unw)=kibondo_names
names(results_kibondo_w)=kibondo_names

posterior_kibondo_unw<- rstan::extract(results_kibondo_unw,
                                      pars= kibondo_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_kibondo_unw=bayesplot::mcmc_trace(posterior_kibondo_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 1, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_kibondo_w<- rstan::extract(results_kibondo_w,
                                    pars= kibondo_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_kibondo_w=bayesplot::mcmc_trace(posterior_kibondo_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 1, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_kibondo_unw,trace_kibondo_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_kibondo.png"), width = 12, height = 12)

# count divergences
sum(sapply(get_sampler_params(results_kibondo_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_kibondo_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# BIT055
######################################################


results_bit055_w<- readRDS(file.path(stanDir,"/stan_BIT055_washed_controlUnw_24.rds"))
results_bit055_unw<- readRDS(file.path(stanDir,"/stan_BIT055_unwashed_24.rds"))

BIT055_names=c("xi[0]", "xi[OlysetPlus]",
               "kappa[0]", "kappa[OlysetPlus]",
               "pi[0]", "pi[OlysetPlus]", 
               "Phi[0]", "Phi[OlysetPlus]",
               "alpha_0","mu_0", "lp__")
BIT055_names_plot=c("pi[OlysetPlus]","kappa[OlysetPlus]","xi[OlysetPlus]")

names(results_bit055_unw)=BIT055_names
names(results_bit055_w)=BIT055_names

posterior_bit055_unw<- rstan::extract(results_bit055_unw,
                                      pars= BIT055_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit055_unw=bayesplot::mcmc_trace(posterior_bit055_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 1, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit055_w<- rstan::extract(results_bit055_w,
                                    pars= BIT055_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_bit055_w=bayesplot::mcmc_trace(posterior_bit055_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 1, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit055_unw,trace_bit055_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit055.png"), width = 12, height = 12)

# count divergences
sum(sapply(get_sampler_params(results_bit055_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit055_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))

######################################################
# BIT103
######################################################


results_bit103_unw<- readRDS(file.path(stanDir,"/stan_BIT103_unwashed_72_ifakara.rds"))
results_bit103_w<- readRDS(file.path(stanDir,"/stan_BIT103_washed_controlUnw_72_ifakara.rds"))

BIT103_names=c("xi[0]", "xi[OlysetPlus]", "xi[IG2]",
               "kappa[0]", "kappa[OlysetPlus]", "kappa[IG2]",
               "pi[0]", "pi[OlysetPlus]", "pi[IG2]", 
               "Phi[0]", "Phi[OlysetPlus]", "Phi[IG2]",
               "alpha_0","mu_0", "lp__")
BIT103_names_plot=c("pi[OlysetPlus]","kappa[OlysetPlus]","xi[OlysetPlus]",
                    "pi[IG2]","kappa[IG2]","xi[IG2]")

names(results_bit103_unw)=BIT103_names
names(results_bit103_w)=BIT103_names

posterior_bit103_unw<- rstan::extract(results_bit103_unw,
                                  pars= BIT103_names_plot,
                                  inc_warmup = TRUE,permuted = FALSE)



trace_bit103_unw=bayesplot::mcmc_trace(posterior_bit103_unw,
                      n_warmup = nwarmup,
                      facet_args = list(nrow = 2, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit103_w<- rstan::extract(results_bit103_w,
                                      pars= BIT103_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit103_w=bayesplot::mcmc_trace(posterior_bit103_w,
                      n_warmup = 7000,
                      facet_args = list(nrow = 2, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit103_unw,trace_bit103_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit103.png"), width = 12, height = 12)

# count divergences
sum(sapply(get_sampler_params(results_bit103_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit103_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))



######################################################
# BIT059
######################################################


results_bit059_w<- readRDS(file.path(stanDir,"/stan_BIT059_washed_controlUnw_24.rds"))
results_bit059_unw<- readRDS(file.path(stanDir,"/stan_BIT059_unwashed_24.rds"))

BIT059_names=c("xi[0]", "xi[OlysetPlus]",
               "kappa[0]", "kappa[OlysetPlus]",
               "pi[0]", "pi[OlysetPlus]", 
               "Phi[0]", "Phi[OlysetPlus]",
               "alpha_0","mu_0", "lp__")
BIT059_names_plot=c("pi[OlysetPlus]","kappa[OlysetPlus]","xi[OlysetPlus]")

names(results_bit059_unw)=BIT059_names
names(results_bit059_w)=BIT059_names

posterior_bit059_unw<- rstan::extract(results_bit059_unw,
                                      pars= BIT059_names_plot,
                                      inc_warmup = TRUE,permuted = FALSE)



trace_bit059_unw=bayesplot::mcmc_trace(posterior_bit059_unw,
                                       n_warmup = nwarmup,
                                       facet_args = list(nrow = 1, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit059_w<- rstan::extract(results_bit059_w,
                                    pars= BIT059_names_plot,
                                    inc_warmup = TRUE,permuted = FALSE)



trace_bit059_w=bayesplot::mcmc_trace(posterior_bit059_w,
                                     n_warmup = nwarmup,
                                     facet_args = list(nrow = 1, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit059_unw,trace_bit059_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit059.png"), width = 12, height = 12)

# count divergences
sum(sapply(get_sampler_params(results_bit059_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit059_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))


######################################################
# KIBONDO 2022
######################################################

results_bit080_w<- readRDS(file.path(stanDir,"/stan_bit080_washed_controlUnw_72.rds"))
results_bit080_unw<- readRDS(file.path(stanDir,"/stan_bit080_unwashed_72.rds"))

bit080_names=c("xi[0]",  "xi[IG2]",
                "kappa[0]", "kappa[IG2]",
                "pi[0]",  "pi[IG2]", 
                "Phi[0]",  "Phi[IG2]",
                "alpha_0","mu_0", "lp__")
bit080_names_plot=c("pi[IG2]","kappa[IG2]","xi[IG2]")

names(results_bit080_unw)=bit080_names
names(results_bit080_w)=bit080_names

posterior_bit080_unw<- rstan::extract(results_bit080_unw,
                                       pars= bit080_names_plot,
                                       inc_warmup = TRUE,permuted = FALSE)



trace_bit080_unw=bayesplot::mcmc_trace(posterior_bit080_unw,
                                        n_warmup = nwarmup,
                                        facet_args = list(nrow = 1, labeller = ggplot2::label_parsed))+
  bayesplot::facet_text(size=15)

posterior_bit080_w<- rstan::extract(results_bit080_w,
                                     pars= bit080_names_plot,
                                     inc_warmup = TRUE,permuted = FALSE)



trace_bit080_w=bayesplot::mcmc_trace(posterior_bit080_w,
                                      n_warmup = nwarmup,
                                      facet_args = list(nrow = 1, labeller = ggplot2::label_parsed), np_style = trace_style_np())+
  bayesplot::facet_text(size=15)

plot_grid(trace_bit080_unw,trace_bit080_w, ncol=1, labels = c("Unwashed", "Washed"))
ggsave(file.path(plotDir, "traceplot_bit080.png"), width = 12, height = 12)

# count divergences
sum(sapply(get_sampler_params(results_bit080_w, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
sum(sapply(get_sampler_params(results_bit080_unw, inc_warmup = FALSE), function(x) sum(x[, "divergent__"])))
