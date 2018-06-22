
#library(multipanelfigure)
library(tidyverse)
library(ggplot2)
library(topicmodels)
library(cowplot)


###############################################################################
# Functions


#' create simulated data for demonstrating LDATS
#' 2 topics, nonuniform distribution of species in two community-types
#'     
#' @param tsteps = number of [monthly] time steps
#' 
#' @return 
#'    beta = matrix of species composition of the groups
#'    gamma = matrix of topic composition over time
#'            3 simulations of gamma: uniform, slow transition, and fast transitions
create_sim_data_ldats = function(tsteps=400) {
    
    topics = 2
    nspecies = 12
    
    # beta: species composition of topics
    # I calculated this distribution by taking the average of each Portal sampling sp distribution (periods 1:436)
    distribution = c(27,13,7, 5, 3, 2, 1, 1, 1, 0, 0, 0)
    # simple permutation of the first distribution
    distribution2 = c(3,1, 0, 1, 0, 13,2, 0, 1,27, 5, 7)
    
    beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
    beta[1,] = distribution/sum(distribution)
    beta[2,] = distribution2/sum(distribution2)
    
    # gamma for a constant topic prevalence through time: topic1 at 90% and topic2 at 10%
    gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
    gamma_constant[,1] = rep(.9,tsteps)
    gamma_constant[,2] = rep(.1,tsteps)
    
    # gamma for a slow gradual transition from topic1 to topic2 
    slow_series = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
    # proportions are constant before change
    slow_series[1:50,1] = rep(1)
    slow_series[1:50,2] = rep(0)
    # duration of change
    slow_series[(50+1):350,1] = seq(300)*(-1/300)+1
    slow_series[(50+1):350,2] = seq(300)*(1/300)+0
    # proportions are constant for rest of time series
    slow_series[(350+1):400,1] = rep(0)
    slow_series[(350+1):400,2] = rep(1)
    
    
    # gamma for two fast transitions
    fast_series = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
    # proportions are constant before change
    fast_series[1:100,1] = rep(1)
    fast_series[1:100,2] = rep(0)
    # change 1: lasts 12 timesteps
    fast_series[(100+1):112,1] = seq(12)*(-.6/12)+1
    fast_series[(100+1):112,2] = seq(12)*(.6/12)+0
    # stable at .6 and .4 for a while
    fast_series[113:250,1] = rep(.4)
    fast_series[113:250,2] = rep(.6)
    # change 2: lasts 12 timesteps
    fast_series[(250+1):262,1] = seq(12)*(.3/12)+.4
    fast_series[(250+1):262,2] = seq(12)*(-.3/12)+.6
    # constant for the rest
    fast_series[263:400,1] = rep(.7)
    fast_series[263:400,2] = rep(.3)
    
    return(list(beta,gamma_constant,slow_series,fast_series))
}

make_pop_plot <- function(pop_data, date_labels = " ", 
                          date_breaks = seq.Date(as.Date("1980-01-01"), as.Date("2010-01-01"), by = "10 years"))
{
    ggplot(data1,aes(x=time,y=n,colour=species)) +
        geom_line() +
        ylab('') +
        scale_x_date(name='', date_labels = date_labels, 
                     breaks = date_breaks) +
        scale_y_continuous(expand = expand_scale(mult = c(0.02, 0.1))) + 
        theme_bw() + 
        theme(legend.position="none", 
              plot.margin = unit(c(0.5, 2.5, 0.5, 0), "mm"), 
              panel.grid.minor = element_line(size = 0.5))
}


make_spcomp_plot <- function(ldamodel, labels = NULL)
{
    species_comp <- data.frame(sp_names, t(exp(ldamodel@beta)))
    names(species_comp) <- c("species", "t1", "t2")
    species_comp %>%
        tidyr::gather(topic, percent, t1:t2) %>%
        ggplot(aes(x = species, y = percent, fill = topic)) +
        geom_bar(stat='identity', 
                 position = position_dodge(width = .33), alpha = 0.8) +
        scale_x_discrete(name='', labels = labels) +
        scale_y_continuous(name='', expand = expand_scale(mult = c(0.02, 0.1))) +
        scale_fill_manual(values=cbPalette[c(1,3)]) +
        theme_bw() + 
        theme(plot.margin = unit(c(0.5, 2.5, 0.5, 0), "mm"), 
              legend.position = 'none')
}

#' @param ldamodel lda model output
#' @param sim_dates xaxis values
#' @param x_labels T/F whether you want to plot x ticks
#' 
#' 
gg_plot_gamma = function(ldamodel,sim_dates,x_labels=F) {
    z = posterior(ldamodel)
    xticks=sim_dates
    ldaplot = data.frame(date=c(),relabund = c(), community = c())
    for (t in 1:2) {
        ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=z$topics[,t],community = as.factor(rep(t,length(z$topics[,1])))))
    }
    if (x_labels==F) {x_text = c('','','','');xname=''} else {x_text=c('1980','1990','2000','2010');xname=''}
    g = ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
        geom_line(size = 1) +
        scale_y_continuous(limits=c(0,1),name='', expand = expand_scale(mult = c(0.02, 0.1))) +
        scale_x_date(name=xname,breaks=c(as.Date('1980-01-01'),as.Date('1990-01-01'),as.Date('2000-01-01'),as.Date('2010-01-01')),labels=x_text) +
        theme_bw() + 
        theme(panel.grid.minor = element_line(size = 0.5), 
              plot.margin = unit(c(0.5, 2.5, 0.5, 0), "mm"), 
              legend.position='none') +
        scale_colour_manual(breaks=as.character(seq(2)),
                            values=cbPalette[c(1,3)],
                            guide=FALSE)
    return(g)
}

##################################################################################

#cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")
cbPalette <- RColorBrewer::brewer.pal(3, "BrBG")

N = 200   # total number of individuals
set.seed(1)

n <- 400
N <- cumsum(sample(c(-10:10), n, TRUE)) + 150

output = create_sim_data_ldats()

# distribution of species in the sample communities follows a species abundance distribution derived from Portal data (average of sampling periods 1:436)
# the distribution of species in the second sample community is a simple permutation of the first
# the two communities both contain 9 out of 12 species
beta = output[[1]]
sp_names <- c('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12')
colnames(beta) <- sp_names

# plot three types of simulated dynamics for the 2 sample communities
gamma_constant = output[[2]]
gamma_slow = output[[3]]
gamma_fast = output[[4]]

sim_dates = seq.Date(from=as.Date('1977-01-01'),by=30,length.out = 400) 

slow =  data.frame(date = rep(sim_dates,dim(gamma_slow)[2]),
                   relabund = as.vector(gamma_slow),
                   community = as.factor(c(rep(1,dim(gamma_slow)[1]),rep(2,dim(gamma_slow)[1]))))
fast =  data.frame(date = rep(sim_dates,dim(gamma_fast)[2]),
                   relabund = as.vector(gamma_fast),
                   community = as.factor(c(rep(1,dim(gamma_fast)[1]),rep(2,dim(gamma_fast)[1]))))
const = data.frame(date = rep(sim_dates,dim(gamma_constant)[2]),
                   relabund = as.vector(gamma_constant),
                   community = as.factor(c(rep(1,dim(gamma_constant)[1]),rep(2,dim(gamma_constant)[1]))))


# create data sets from beta and gamma; data must be in integer form (simulating species counts)
dataset1 = round(as.data.frame(gamma_slow %*% beta) *N,digits=0)
dataset2 = round(as.data.frame(gamma_fast %*% beta) *N,digits=0)
dataset3 = round(as.data.frame(gamma_constant %*% beta) *N,digits=0)

dataset1$time = sim_dates
dataset2$time = sim_dates
dataset3$time = sim_dates

# plot species pop over time for 3 datasets
data1 <- tidyr::gather(dataset1, species, n, S1:S12, factor_key=TRUE)
data2 <- tidyr::gather(dataset2, species, n, S1:S12, factor_key=TRUE)
data3 <- tidyr::gather(dataset3, species, n, S1:S12, factor_key=TRUE)

pop1 <- make_pop_plot(data1)
pop2 <- make_pop_plot(data2)
pop3 <- make_pop_plot(data3, date_labels = "%Y")

#############################################################################################
# Run LDA
SEED  = 1

ldamodel1 = LDA(dataset1[,-13],k=2, control = list(seed = SEED,estimate.alpha=F,alpha=.1),method='VEM')
ldamodel2 = LDA(dataset2[,-13],k=2, control = list(seed = SEED,estimate.alpha=F,alpha=10),method='VEM')
ldamodel3 = LDA(dataset3[,-13],k=2, control = list(seed = SEED,estimate.alpha=F,alpha=.1),method='VEM')

beta1 <- exp(ldamodel1@beta)
beta2 <- exp(ldamodel2@beta)
beta3 <- exp(ldamodel3@beta)

#  plot gammas
g1 = gg_plot_gamma(ldamodel1,sim_dates)
g2 = gg_plot_gamma(ldamodel2,sim_dates)
g3 = gg_plot_gamma(ldamodel3,sim_dates, TRUE)

spcomp1 <- make_spcomp_plot(ldamodel1)
spcomp2 <- make_spcomp_plot(ldamodel2)
spcomp3 <- make_spcomp_plot(ldamodel3, labels = LETTERS[1:12])
print(spcomp3)

# Adjust formatting on panels

pop1 <- pop1 + labs(y = "Gradual Change", title = "Community Time Series") + 
    theme(plot.title = element_text(hjust = 0.5, size = 11))
pop2 <- pop2 + labs(y = "Rapid Change")
pop3 <- pop3 + labs(y = "Steady State")
gg_spcomp1 <- gg_spcomp1 + labs(title = "Species Assemblages") + 
    theme(plot.title = element_text(hjust = 0.5, size = 11))
g1 <- g1 + labs(title = "Biodiversity Dynamics") + 
    theme(plot.title = element_text(hjust = 0.5, size = 11))

plot_final <- plot_grid(pop1, gg_spcomp1, g1, 
                        pop2, gg_spcomp2, g2, 
                        pop3, gg_spcomp3, g3, nrow = 3, 
                        align = "hv")
ggsave("figure_LDATS.pdf", plot_final, width = 7.5, height = 4.5)
