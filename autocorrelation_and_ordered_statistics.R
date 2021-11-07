#' ---
#' title: "Rank statistics with correlated draws"
#' author: "Aki Vehtari"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     theme: readable
#'     toc: true
#'     toc_depth: 2
#'     toc_float: true
#'     code_download: true
#' ---

#' ## Introduction
#' 
#' Rank statistics with correlated draws is not trivial. This notebook
#' illustrates the problematic behavior, but also why it's not a
#' problem in SBC. The results are useful, for example, in the context
#' of SBC, PPC, LOO-PIT and MCMC convergence diagnostics. The related
#' paper is *Graphical Test for Discrete Uniformity and its
#' Applications in Goodness of Fit Evaluation and Multiple Sample
#' Comparison* by Säilynoja, Bürkner and Vehtari (2021).
#'
#' If we have draws $y \sim g(y)$ and draws $x \sim p(x)$, and want to
#' assess whether $g=p$, we can use discrete uniformity test of rank
#' statistics.
#' 

#+ setup, include=FALSE
knitr::opts_chunk$set(message=FALSE, error=FALSE, warning=FALSE, comment=NA,
                      cache=TRUE)
library(posterior)
library(tibble)

#' ## Bias in extreme ranks when using the same Markov chain many times
#'
#' The first example illustrates how extreme ranks are overrepresented
#' if independent draws $y \sim g(y)$ are compared to dependent draws
#' $x \sim p(y)$. We simulate MCMC with an AR process with phi=0.95
#' that will give relative efficiency of 2.5% compared to the
#' independent draws. The dependent draws tend to be more clustered
#' and have smaller support than corresponding number of independent
#' draws. The issue of distribution of maximum value of correlated
#' observations is well known in extreme value analysis.
#'
#' This example is related to PPC and LOO-PIT. $y$ would be the
#' observed data and $x$ would be dependent draws from the posterior
#' predictive or LOO predictive distribution.
#' 

#' Rank statistic function
rank_stat = function(element, vec) { sum(vec < element) + 1 }
#' We choose N to be 2^a-1 to avoid binning artifacts
N = 16384-1
#' Independent draws from the distribution g(y) = normal(y|0,1)
y = rnorm(N, 0, 1)
#' Pre-allocate storage
ranks = matrix(nrow=1000, ncol=N)
ranksthin = matrix(nrow=1000, ncol=N)
#' 1000 repetitions
for (j in 1:1000) {

  # arima.sim to simulate dependent draws x from AR process that has
  # marginal distribution p(x) = normal(y|0,1).
  # phi=0.95 gives average ESS of 417 for the sample size 16383 (2.5%).
  # Here we generate new dependent draws $x$ for each $y_j$
  phi = 0.95
  x = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))

  # rank stats for all y with respect to dependent draws x
  ranks[j,] = sapply(y, function(ele) { rank_stat(ele, x)})

  # rank stats for all y with respect to thinned draws x (with thinning
  # 64 the draws have relative efficiency of 90%, ie, almost
  # independent)
  ranksthin[j,] = sapply(y, function(ele) { rank_stat(ele, x[seq(1,N,by=64)])})

}

#' When x are correlated, we see spikes in the extreme ranks
hist(ranks, breaks=seq(0.5,N+1+0.5,length=N+1+1))

#' When x are correlated, we see spikes in the extreme ranks,
#' zoom to ranks 1:100
hist(ranks[ranks<=100], breaks=seq(0.5,100+1+0.5,length=100+1+1))

#' When x are thinnd to be almost independent, we see no spikes in the
#' extreme ranks
hist(ranksthin, breaks=seq(0.5,256+1+0.5,length=256+1+1))

#' When x are thinnd to be almost independent, we see no spikes in the
#' extreme ranks, zooom to ranks 1:100
hist(ranksthin[ranksthin<=100], breaks=seq(0.5,100+1+0.5,length=100+1+1))

#' ## Bias in extreme ranks when using one Markov chain only once
#'
#' The second example illustrates the behavior in SBC, where we
#' compare one $y_j \sim g(y)$ to dependent draws $x \sim p(y)$, but
#' for each $y_j$ we generate new dependent draws $x$. This means that
#' the rank statistics of $y_j$ are independent from each other, but
#' the rank statistics are still biased.
#' 

#' 100 repetitions
sx = matrix(nrow=1000, ncol=N)
sy = matrix(nrow=1000, ncol=N)
sranks = matrix(nrow=1000, ncol=N)
sranksthin = matrix(nrow=1000, ncol=N)
sranksthin16 = matrix(nrow=1000, ncol=N)
for (k in 1:100) {
  print(k)
  # Rank statistic function
  rank_stat = function(element, vec) { sum(vec < element) + 1 }
  # We choose N to be 2^a-1 to avoid binning artifacts
  N = 16384-1
  # Independent draws from the distribution g(y) = normal(y|0,1).
  # In SBC, these draws are from the prior.
  y = rnorm(N, 0, 1)
  # Pre-allocate storage
  tbsbc = tibble(draw=0,
                 rank=0,
                 rankthin=0,
                 ess=0,
                 esstail=0,
                 essthin=0,
                 .rows=N)
  # Repeat for all $y_j$
  for (j in 1:N) {
    
    # arima.sim to simulate dependent draws x from AR process that has
    # marginal distribution p(x) = normal(y|0,1).
    # phi=0.95 gives average ESS of 417 for the sample size 16383 (2.5%).
    # Here we generate new dependent draws $x$ for each $y_j$
    phi = 0.95
    x = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))
    tbsbc$draw[j]=x[N]
    
    # Rank stat for one draw $y_j \sim g(y)$ (ie prior in SBC) with
    # respect to dependent draws x
    tbsbc$rank[j]=rank_stat(y[j], x)
    
    # Rank stat for one draw $y_j \sim g(y)$ (ie prior in SBC) with
    # respect to thinned draws x.
    # Thinning by 64 gives average ESS of 228 for the sample size 256 (89%)
    tbsbc$rankthin[j]=rank_stat(y[j], x[seq(1,N,by=64)])
    tbsbc$rankthin16[j]=rank_stat(y[j], x[seq(1,N,by=16)])
    
    # ESSs
    tbsbc$ess[j]=ess_basic(x)
    tbsbc$esstail[j]=ess_tail(x)
    tbsbc$essthin[j]=ess_basic(x[seq(1,N,by=64)])
    tbsbc$essthin16[j]=ess_basic(x[seq(1,N,by=16)])
  }
  sx[k,] = tbsbc$draw
  sy[k,] = y
  sranks[k,] = tbsbc$rank
  sranksthin[k,] = tbsbc$rankthin
  sranksthin16[k,] = tbsbc$rankthin16
}

#' Rank statistics are independent but biased and we see spikes
shist=matrix(nrow=100,ncol=N+1);
for (i in 1:100) {
  shist[i,]=hist(sranks[i,], breaks=seq(0.5,N+1+0.5,length=N+1+1),plot=F)$counts
}
qplot(1:16384,colMeans(shist))

#' Thinning reduces the bias so that we don't see spikes
shistthin=matrix(nrow=100,ncol=N+1);
for (i in 1:100) {
  shistthin[i,]=hist(sranksthin[i,], breaks=seq(0.5,N+1+0.5,length=N+1+1),plot=F)$counts
}
qplot(1:16384,colMeans(shistthin))

#' ## Bias in extreme ranks when using the same Markov chain many times
#'
#' The third example illustrates if dependent draws $y \sim g(y)$ are
#' compared to similarly dependent draws $x \sim p(y)$.  We simulate
#' MCMC with an AR process with phi=0.95 that will give relative
#' efficiency of 2.5% compared to the independent draws.
#'
#' This example is related to MCMC rank plots. $y$ and $x$ are from
#' two independent Markov chains with the same target distribution.
#' 

#' Rank statistic function
rank_stat = function(element, vec) { sum(vec < element) + 1 }
#' We choose N to be 2^a-1 to avoid binning artifacts
N = 16384-1
#' Pre-allocate storage
mranks = matrix(nrow=1000, ncol=N)
mranksthin = matrix(nrow=1000, ncol=N)
#' 1000 repetitions
for (j in 1:1000) {

  # arima.sim to simulate dependent draws y from AR process that has
  # marginal distribution g(y) = normal(y|0,1)
  # phi=0.95 gives average ESS of 417 for the sample size 16383 (2.5%).
  phi = 0.95
  y = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))
  # arima.sim to simulate dependent draws x from AR process that has
  # marginal distribution g(y) = normal(y|0,1)
  # phi=0.95 gives average ESS of 417 for the sample size 16383 (2.5%).
  phi = 0.95
  x = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))
  
  # rank stats for dependent draws y with respect to dependent draws x
  mranks[j,] = sapply(y, function(ele) { rank_stat(ele, x)})
  
  # rank stats for thinned y with respect to thinned draws x (with thinning
  # 64 the draws have relative efficiency of 90%, ie, almost
  # independent)
  mranksthin[j,] = sapply(y, function(ele) { rank_stat(ele, x[seq(1,N,by=64)])})
  
}

#' When x are correlated, we see spikes in the extreme ranks
hist(mranks, breaks=seq(0.5,N+1+0.5,length=N+1+1))

#' When x are correlated, we see spikes in the extreme ranks,
#' zoom to ranks 1:100
hist(mranks[mranks<=100], breaks=seq(0.5,100+1+0.5,length=100+1+1))

#' When x are thinnd to be almost independent, we see no spikes in the
#' extreme ranks
hist(mranksthin, breaks=seq(0.5,256+1+0.5,length=256+1+1))

#' When x are thinnd to be almost independent, we see no spikes in the
#' extreme ranks, zooom to ranks 1:100
hist(mranksthin[mranksthin<=100], breaks=seq(0.5,100+1+0.5,length=100+1+1))


#' ## Ilustration of the bias in the ordered statistics
#'
#' 
N = 1000

#' store 100 first ordered values of iid normal
ordn=matrix(nrow=1000,ncol=100);
for (i in 1:1000) {
  ordn[i,]=sort(rnorm(N,0,1))[1:100]
};
#' 10 first expected ordered statistics for iid normal
colMeans(ordn)

#' store 100 first ordered values of AR(phi=0.95)
ordar=matrix(nrow=1000,ncol=100);
phi=0.95
for (i in 1:1000) {
  ordar[i,]=sort(arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))[seq(1,N,by=1)])[1:100]
}
#' 10 first expected ordered statistics
colMeans(ordar)[1:10]

#' store 100 first ordered values of AR(phi=-0.95)
ordarm=matrix(nrow=1000,ncol=100);
phi= - 0.95
for (i in 1:1000) {
  ordarm[i,]=sort(arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))[seq(1,N,by=1)])[1:100]
}
#' 10 first expected ordered statistics
colMeans(ordarm)[1:10]

#' store 100 first ordered values of AR(phi=0.95) thinned by mean ess
ordar.thin=matrix(nrow=1000,ncol=100);
phi=0.95
for (i in 1:1000) {
  ordar.thin[i,]=sort(arima.sim(n = 18 * N, list(ar = c(phi)), sd = sqrt((1-phi^2)))[seq(1,18 * N,by=18)])[1:100]
}
#' 10 first expected ordered statistics
colMeans(ordar.thin)[1:10]

#' store 100 first ordered values of AR(phi=-0.95)
ordarm.thin=matrix(nrow=1000,ncol=100);
phi= - 0.95
for (i in 1:1000) {
  ordarm.thin[i,]=sort(arima.sim(n = 7 * N, list(ar = c(phi)), sd = sqrt((1-phi^2)))[seq(1,7 * N,by=7)])[1:100]
}
#' 10 first expected ordered statistics
colMeans(ordarm.thin)[1:10]

# plot 100 first ordered statistics for normal and AR(phi=0.95)
ggplot() + geom_abline(alpha=.5, size=1)+
  geom_point(mapping=aes_(x = colMeans(ordn), y=colMeans(ordar), colour="A"), size=4)+
  labs( # title="First 100 order statistics",
        x="Normal(0, 1)",
        y=TeX("AR($\\varphi = 0.95$)"),
        ) +
  theme(text = element_text(size=22)) +
  scale_y_continuous(limits=c(-3.3,-1)) +
  theme(legend.position = "none") +
  scale_color_bright()
ggsave("figures/order_statistics_sigma_95.pdf", width=6, height=5)

# plot 100 first ordered statistics for normal and AR(phi=-0.95)
ggplot() + geom_abline(alpha = .5, size = 1)+
  geom_point(mapping=aes_(x = colMeans(ordn), y=colMeans(ordarm), colour="A"), size=4)+
  labs( # title="First 100 order statistics",
        x="Normal(0, 1)",
        y=TeX("$AR(\\varphi = -0.95)$")) +
  theme(text = element_text(size=22)) +
  scale_y_continuous(limits=c(-3.3,-1)) +
  theme(legend.position = "none") +
  scale_color_bright()
ggsave("figures/order_statistics_sigma_-95.pdf", width=6, height=5)
# plot 100 first ordered statistics for normal and AR(phi=0.95)
ggplot() + geom_abline(alpha=.5, size = 1) +
  geom_point(mapping=aes_(colMeans(ordn),colMeans(ordar.thin), color="a"), size=4)+
  labs(# title= "First 100 order statistics",
    x="Normal(0, 1)",
    y=TeX("Thinned $AR(\\varphi = 0.95)$")) +
  theme(text = element_text(size=22)) +
  scale_y_continuous(limits=c(-3.3,-1)) +
  theme(legend.position = "none") +
  scale_color_bright()
ggsave("figures/order_statistics_sigma_95_thinned.pdf", width=6, height=5)

# plot 100 first ordered statistics for normal and AR(phi=-0.95)
ggplot() + geom_abline(alpha=.5, size = 1) +
  geom_point(mapping=aes_(x = colMeans(ordn), y = colMeans(ordarm.thin), color="a"), size=4)+
  labs(# title= "First 100 order statistics",
    x="Normal(0, 1)",
    y=TeX("Thinned $AR(\\varphi = -0.95)$")) +
  theme(text = element_text(size=22)) +
  scale_y_continuous(limits=c(-3.3,-1)) +
  theme(legend.position = "none") +
  scale_color_bright()
ggsave("figures/order_statistics_sigma_-95_thinned.pdf", width=6, height=5)

