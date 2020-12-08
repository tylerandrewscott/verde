
require(data.table)
require(pbapply)
require(statnet)
require(lsa)
library(fields)
library(ggraph)
require(Matrix)
require(purrr)
require(tidyverse)
lsa::cosine()
#specify working parameters
n_agents = 40
n_issues = 100

alpha = runif(1,0,1)
beta = runif(1,0,1)
gamma = 1 - alpha - beta
beta1 = 0.90
beta2 = 0.05
periods = 1
n_info_pieces = 4
info_uncertainty = 0.1
payoff_uncertainty = 1
n_agent_per_issue = n_agents
range_of_motivation = c(0.1,0.9)
range_of_communication = c(0.1,0.9)
range_of_capacity = c(0.1,0.9)
min_payout = 0;max_payout = 10
createCGR = TRUE





## robustness:
#knockout nodes
#shuffle payoffs

#CGR DRIVERS:
#initiating leadership: smart invitations?
#uncertainty: sampling distribution for payout information

### THESE TWO ARE SORT OF IMPLICITY IN THE PUBLIC GOODS GAME -- MORE NEEDED?
#consequential incentives: needed? maybe not -- players do win and lose here -- how can we adjust BATNA?
#interdependence: how can we tie goods together?

#initial agent attributes
#agent_attributes = list(motivation = runif(n_agents,range_of_motivation[1],range_of_motivation[2]),
#                        communication = runif(n_agents,range_of_communication[1],range_of_communication[2]),
#                        capacity = runif(n_agents,range_of_capacity[1],range_of_capacity[2]))
agent_attribute_start = runif(n_agents,range_of_motivation[1],range_of_motivation[2])
agent_attributes = list(motivation = agent_attribute_start,
                        communication = agent_attribute_start,
                        capacity = agent_attribute_start)


#initialize agent network
agent_net = network(matrix(0,ncol=n_agents,nrow=n_agents),directed = T,vertex.attr = agent_attributes)
#initialize agent-issue network with placeholder 1 and 0 values to fix density


issue_agent_list = lapply(seq(n_issues),function(x) sort(sample(x = seq(n_agents),size = n_agent_per_issue,replace = F)))


agent_issue_matrix = sapply(issue_agent_list,function(x) {x %in% seq(n_agents)}+0)

# start_agent_issue_net = cbind(matrix(1,ncol = n_agent_focal_issues,nrow = n_agents),matrix(0,ncol = n_issues - n_agent_per_issue,nrow = n_agents))
# agent_issue_network = do.call(cbind,replicate(n_agents,sample(c(rep(1,n_agent_focal_issues),rep(0,n_issues - n_agent_per_issue)),replace = F),simplify = F))

# agent_issue_network = simulate(as.network(agent_issue_network,bipartite = T,directed = F,matrix.type = 'incidence')~ edges,
#                                constraints = ~edges + bd(maxin =10,maxout = 10,minout = 5))


agent_issue_network = agent_issue_matrix
#compute agent homophily based on shared issues                     
agent_issue_homophily <- cosine(t(agent_issue_network))
#compute issue homophily based on shared agents
issue_agent_homophily <- cosine(agent_issue_network)

seed_network = network(rgraph(n_agents,tprob = 0.04),directed = F)
sqrt_cosim = sqrt(agent_issue_homophily)
diag(sqrt_cosim) <- 0

if(n_agent_per_issue==n_agents){
  agent_network = (simulate(seed_network ~ edges + isolates + gwdegree(0.5,fixed = T) + gwesp(0.5,fixed = T),coef = c(-2,-Inf,-0.75,0.75),constraints = ~edges))
}
if(n_agent_per_issue!=n_agents){
  agent_network = (simulate(seed_network ~ edges + isolates + gwdegree(0.5,fixed = T) + gwesp(0.5,fixed = T) + 
                            edgecov(sqrt_cosim),coef = c(-2,-Inf,-0.75,0.75,1),constraints = ~edges))
}

# issue_network = rgraph(n = nrow(issue_agent_homophily),tprob = 0.025, mode = 'graph')
# if(n_agent_per_issue==n_agents){
# issue_network = (simulate(network(issue_network,directed = F) ~ edges + isolates,coef = c(-2,-Inf),
#                           constraints = ~ edges))
# }
# if(n_agent_per_issue!=n_agents){
#   issue_network = (simulate(network(issue_network,directed = F) ~ edges + edgecov(sqrt(issue_agent_homophily)) + isolates,coef = c(-5,2,-Inf),
#                             constraints = ~ edges + bd(minin = 2)))
# }
issue_payoffs = runif(n = n_issues,min_payout,max_payout)

require(abind)


agent_edge_index = data.table(which(as.sociomatrix(agent_network)==1,arr.ind = T))
agent_alters = lapply(seq(n_agents),function(x) agent_edge_index[row == x,]$col)

agent_issue_edge_index = data.table(which(agent_issue_network==1,arr.ind = T))
agent_issues = lapply(seq(n_agents),function(x) agent_issue_edge_index[row == x,]$col)



agent_payoff_estimates = payoff_estimates_matrix  = buyin_matrix = expected_contrib = expected_payout  = group_expected_payout = contribution_matrix = contributor_matrix = payout_matrix = reciprocity_matrix  = array(NA,dim = c(n_agents,n_issues,periods))
true_payoff_matrix = matrix(NA,nrow = n_issues,ncol = periods)
motivation_matrix = engagement_matrix = capacity_matrix = matrix(NA,nrow = n_agents,ncol = periods)

motivation_matrix[,1] <- capacity_matrix[,1] <- engagement_matrix[,1] <- agent_attributes$motivation 


for(i in seq(n_issues)){
  true_payoff_matrix[i,] <- issue_payoffs[i] + arima.sim(n = periods, model = list(ar = c(1,-0.1)),#list(ar = c(1, -0.4858), ma = c(-0.2279, 0.2488)),
        sd = 0.25,n.start = 100)}


#### FIX: TIE GOODS CONTRIBUTIONS AND PAYOUTS TO BIPARTITE NETWORK
### CGR INNOVATION IS TO PAIR THINGS UP

#CGR_issues = sample(seq(n_issues),20)
#CGR_participants = sample(seq(n_agents),20)

#communication --- sharing and learning distortion
#shared motivation --- willingness to contribute
#joint action --- thinking like a team (group payout)


agent_issue_binary = as.sociomatrix(agent_issue_network)

for(t in seq(periods)){
  print(t)
  agent_info_matrices  = replicate(n_agents,sapply(true_payoff_matrix[,t],rnorm, n = n_info_pieces,sd =info_uncertainty),simplify = F)
  agent_info_shared = do.call(rbind,lapply(agent_info_matrices,function(x) colMeans(cbind(x))))

  
  agent_info_received = mapply(function(x,y) x * y, 
                               x = lapply(agent_alters,function(x) colMeans(matrix(agent_info_shared[x,],ncol = n_issues,byrow = T))), 
                               y = as.list(sqrt(engagement_matrix[,t])),SIMPLIFY = F)

  agent_post_sharing_payoff_estimates = mapply(function(x,y) colMeans(rbind(x,y),na.rm = T),x = agent_info_matrices,y = agent_info_received,SIMPLIFY = F)
  payoff_estimates = do.call(rbind,agent_post_sharing_payoff_estimates)
  payoff_estimates_matrix[,,t] <- payoff_estimates

  issue_rank_t = rank(-rowSums(apply(payoff_estimates_matrix[,,t],1,rank)),ties.method = 'first')
  
  for(ir in seq(n_issues)){
    i = which(issue_rank_t==ir)
    #update estimates
    #agents revise up or down based on their payout from last time
  #  if(t == 1){payoff_estimates_matrix[,i,t]<-ifelse(is.na(payoff_estimates[,i]),0,payoff_estimates[,i])}
  #  if(t != 1){payoff_estimates_matrix[,i,t] <- payoff_estimates_matrix[,i,t-1]} 
    # + {payout_matrix[,i,t-1] - contribution_matrix[,i,t-1]}} 
    #contribution = p(random>motivation) * estimated payoff * motivation
    
    expected_contrib[,i,t] <- payoff_estimates_matrix[,i,t]  * motivation_matrix[,t] * motivation_matrix[,t]
    expected_payout[,i,t] <- sum(expected_contrib[,i,t]) *  payoff_estimates_matrix[,i,t] / n_agents
 
    contribution_matrix[,i,t]<- (runif(1,min = 0,max = 1) < motivation_matrix[,t]) * payoff_estimates_matrix[,i,t]  * motivation_matrix[,t]
    contributor_matrix[,i,t] <- (contribution_matrix[,i,t]>0)+0
    payout_matrix[,i,t] <- {sum(contribution_matrix[,i,t] * true_payoff_matrix[i,t] ) / {sum(agent_issue_binary[,i])}} * agent_issue_binary[,i]
    #reciprocity_matrix[,i,t] <- payout_matrix[,i,t] > {contribution_matrix[,i,t] * payoff_estimates_matrix[,i,t] * motivation_matrix[,i,t]}
    #reciprocity_matrix[,i,t] <- payout_matrix[,i,t] >= expected_payout[,i,t] - motivation_matrix[,t]
    reciprocity_matrix[,i,t] <- contribution_matrix[,i,t] * payoff_estimates_matrix[,i,t] * motivation_matrix[,t] < payout_matrix[,i,t]
    buyin_matrix[,i,t] <- sapply(seq(n_agents),function(x) sum(contributor_matrix[-x,1:i,t])/({n_agents-1}*i))
    if(t!=periods){
    motivation_matrix[,t+1] = motivation_matrix[,t] * beta1 + rowSums(cbind(reciprocity_matrix[,i,t])) * beta2
    engagement_matrix[,t+1] = engagement_matrix[,t] * beta1 + rowSums(cbind(contribution_matrix[,i,t]>0))/(t*n_agents) * beta2
    }
  }
}





pg_results = data.table(melt(payout_matrix))

goods_by_time = pg_results[,sum(value,na.rm = T),by=.(Var2,Var3)]
setnames(goods_by_time,c('Var2','Var3','V1'),c('Issue','Time','goods'))

payoffs = data.table(melt(true_payoff_matrix))
setnames(payoffs,c('Var1','Var2','value'),c('Issue','Time','multiplier'))


goods_by_time = merge(goods_by_time,payoffs,all = T)

ggplot(goods_by_time,
       aes(x = Time,y = goods,group = Issue)) + 
  geom_path(aes(col = multiplier)) + 
  geom_point(aes(col = multiplier)) + #geom_point() + 
  scale_color_viridis_c() + 
  ylab('Public goods produced') + xlab('time period') + 
  ggtitle('100 public goods over 20 periods','no CGR inserted') 



ggplot(melt(issue_payoff_matrix),aes(x = Var2,col = as.factor(Var1),y = value)) + geom_point()


mot_melt = melt(motivation_matrix)
mot_melt$dynamic = 'motivation'
eng_melt = melt(engagement_matrix)
eng_melt$dynamic = 'engagement'
dynamic_melt = data.table(rbind(mot_melt,eng_melt))

ggplot(dynamic_melt,aes(x = Var2,
                    y = value,col = dynamic)) + 
  geom_jitter()

issue_payoffs
issue_payoff_matrix[10,]
lapply(seq(periods),function(x) table(contribution_matrix[,,x] < payout_matrix[,,x]))


test = data.table(melt(motivation_matrix))

ggplot(data = test[Var3%in% c(1,100),]) + geom_boxplot(aes(x = as.factor(Var2),y = value)) + 
  facet_wrap(~Var3)



melt(contribution_matrix)

table(motivation_matrix[,,27])

payout_matrix - contribution_matrix
test = data.table(melt(contribution_matrix>0))

table(agent_issue_binary[,1])

ggplot(test[Var2==1,]) + geom_tile(aes(x = as.factor(Var1),y = as.factor(Var3),fill = value))
summary(motivation_matrix[,,100])

rowSums(agent_issue_binary)





goods_by_time = pg_results[,sum(value,na.rm = T),by=.(Var2)]
setnames(goods_by_time,c('Var2','V1'),c('Issue','goods'))

payoffs = data.table(melt(true_payoff_matrix))
setnames(payoffs,c('Var1','Var2','value'),c('Issue','Time','multiplier'))

goods_by_time = merge(goods_by_time,payoffs,all = T)


ggplot(goods_by_time,
       aes(x = rank(multiplier),y = goods)) + 
  geom_point() +
  scale_color_viridis_c() + 
  ylab('Public goods produced') + xlab('issue') + 
  ggtitle('100 public goods over 20 periods','no CGR inserted') 


