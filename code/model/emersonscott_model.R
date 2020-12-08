
packs = c('data.table','tidyverse','lsa','statnet','igraph','multinets','tidygraph','intergraph','tidyr')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

###### 3.1 Creating the Policy System ######

##### functions for 3.1 ######
createAgents = function(simulation.control) {c(rep('stakeholders',simulation.control$stakeholders),
                                    rep('regulators',simulation.control$regulators),
                                    rep('convenors',simulation.control$convenors))}
#createIssues = function(control) {runif(n = simulation.control$n_issues,min = simulation.control$min.payout,max = simulation.control$max.payout)}
createMotivation = function(simulation.control) {runif(length(agents), simulation.control$bottom.incentive,simulation.control$top.incentive)}


createAgentInformationMatrix = function(payoffs =true.values,
                                        u = simulation.control$uncertainty,
                                        n_pieces = simulation.control$n_pieces){
    infoFull = do.call(cbind,sapply(seq_along(payoffs),function(x) 
    rnorm(n = simulation.control$n_pieces*length(agents),mean = payoffs[x],sd = u),simplify = F))
    infoFull[infoFull<0] <- 0
    agentInfo = rbind(infoFull[n_pieces,])
  return(agentInfo)}


# createCGR = function(select = NULL,number.of.issues.to.join = 2,issue.network = issue.network){
#   issue.graph = intergraph::asIgraph(issue.network)
#   if(select=='betweenness'){issue.weight = igraph::betweenness(issue.graph)}
#   picked.issue = sample(seq(simulation.control$n_issues),size = 1,prob = issue.weight)
#   issues = ego(issue.graph,nodes = picked.issue)
#   agents = which(rowSums(agent.issue.incidence.matrix[,issues[[1]]])>=number.of.issues.to.join)
#   list('issues' = sort(issues[[1]]),'agents' = sort(agents))
#   }

EmersonScottModel = function(simulation.control = NULL) {
#### actions for 3.1 #######
agents = createAgents(simulation.control = simulation.control)
true.values = {runif(n = simulation.control$n_issues,min = simulation.control$min.payout,max = simulation.control$max.payout)}
agent.motivation = createMotivation(simulation.control = simulation.control)

agent.issue.incidence.matrix = matrix(0+(runif(simulation.control$n_issues*length(agents))>0.75),ncol = simulation.control$n_issues)
issue.homophily <- lsa::cosine(t(agent.issue.incidence.matrix))
agent.homophily <- lsa::cosine(agent.issue.incidence.matrix)

seed_network = network(rgraph(length(agents),tprob = 0.04),directed = F)
sqrt.agent.homophily = sqrt(agent.homophily)
agent.network = (simulate(seed_network ~ edges + isolates + gwdegree(0.5,fixed = T) + gwesp(0.5,fixed = T) + 
                            edgecov(sqrt.agent.homophily),coef = c(-2,-Inf,-0.75,0.75,1),constraints = ~edges))

seed_network2 = network(rgraph(length(true.values),tprob = 0.02),directed = F)
sqrt.issue.homophily = sqrt(issue.homophily)
issue.network = (simulate(seed_network2 ~ edges + isolates + gwdegree(0.5,fixed = T) + gwesp(0.5,fixed = T) + edgecov(sqrt.issue.homophily),
                          coef = c(-2,-Inf,-0.75,0.75,1),constraints = ~edges))

agent.edges=melt(as.sociomatrix(agent.network))
issue.edges=melt(as.sociomatrix(issue.network))
issue.edges$Var1 <- issue.edges$Var1 + length(agents)
issue.edges$Var2 <- issue.edges$Var2 + length(agents)

issue.agent.edges = data.table(melt(agent.issue.incidence.matrix))[value>0,]
issue.agent.edges$Var2 <- issue.agent.edges$Var2 + length(agents)


complete.edges = rbindlist(list(agent.edges,issue.edges,issue.agent.edges))[value == 1&Var1!=Var2,]

#create CGR

issue.graph = intergraph::asIgraph(issue.network)
if(simulation.control$CGselector=='betweenness'){issue.weight = igraph::betweenness(issue.graph)}
picked.issue = sample(seq(simulation.control$n_issues),size = 1,prob = issue.weight)
cgr = list('issues' = sort(ego(issue.graph,nodes = picked.issue)[[1]]),'agents' = sort(which(rowSums(agent.issue.incidence.matrix[,issues[[1]]])>=simulation.control$number.of.issues.to.join)))

#cgr = createCGR(select = 'betweenness',number.of.issues.to.join = 2)

###### 3.2 Running the Policy System ######
##### functions for 3.2 ######

shared.info.matrix = private.info.matrix = issue.matrix = incentive.matrix = contrib.matrix = orig.contrib.matrix = payout.matrix = reciprocity.matrix  = array(NA,dim = list(length(agents),length(true.values),simulation.control$t))
contributors.matrix = array(NA,dim = list(1,length(true.values),simulation.control$t))
incentive.matrix[,,1] <- agent.motivation
issue.matrix[,,1] <- agent.issue.incidence.matrix

dynamics.tracker = list('principled.engagement'=NULL,'capacity.for.joint.action'=NULL,'shared.motivation'=NULL)

for(t in 1:simulation.control$t){
  info.pieces = replicate(length(agents),createAgentInformationMatrix(),simplify = F)
  ### 3.2.1 Agents share information #### 
  private.info.matrix[,,t]  = do.call(rbind,info.pieces)
  shared.info.matrix[,,t] = private.info.matrix[,,t]  / sqrt(incentive.matrix[,,t])
  if(t != 1){
  private.info.matrix[,,t][issue.matrix[,,t]==1]<-NA
  shared.info.matrix[,,t][issue.matrix[,,t]==1]<-NA
  }
  all.private.info = data.table(melt(private.info.matrix))
  setnames(all.private.info,c('Var1','Var2','Var3'),c('Agent','Issue','Time'))
  
  all.shared.info = data.table(melt(shared.info.matrix))
  setnames(all.shared.info,c('Var1','Var2','Var3'),c('Sharing.Agent','Issue','Time'))
  all.shared.info = all.shared.info[,mean(value,na.rm = T),by=.(Sharing.Agent,Issue)]
  agent.edges = data.table(agent.edges)
  agent.alters = lapply(seq_along(agents),function(a) agent.edges[value==1&Var1==a,]$Var2)

  
#agent network info
personal.network.info.t = do.call(rbind,lapply(seq(length(agents)),function(x) all.shared.info[Sharing.Agent %in% agent.alters[[x]],mean(V1,na.rm = T),by=.(Issue)]$V1))

agent.private.info.t = dcast(all.private.info[,mean(value,na.rm = T),by=.(Agent,Issue)],Agent ~ Issue,value.var = 'V1')
agent.private.info.t[,Agent:=NULL]

### 3.2.2 Agents make investment decisions ####
### 3.2.3 Agents respond and adjust ####
  if(min(c(incentive.matrix[,,t]))<0){
  print(t);print(summary(c(incentive.matrix[,,t])))}
  #calculate principled engagement at time t, used as basis for information sharing
  dynamics.tracker$principled.engagement[t] <- mean(incentive.matrix[cgr$agents,cgr$issues,t])
 
  #cgr info
    cgr.info = agent.private.info.t + incentive.matrix[,,t] * (1-dynamics.tracker$principled.engagement[t])
    cgr.payout.estimates = colMeans(cgr.info[cgr$agents,])
    info.in.cgr = (agent.private.info.t + personal.network.info.t + cgr.payout.estimates) / 3
    info.out.cgr = (agent.private.info.t + personal.network.info.t)/2
    operating.info = do.call(rbind,lapply(seq(length(agents)),function(a) {if(a %in% cgr$agents){info.in.cgr[a,]}else{info.out.cgr[a,]}}))
    agent.best.guesses = operating.info
    agent.issue.knowledge.pieces = lapply(info.pieces,function(x) colSums(!is.na(x),na.rm = T))
    agent.binary.issue.knowledge = lapply(info.pieces,function(x) (colSums(!is.na(x),na.rm = T)>0)+0)

    draws = do.call(rbind,lapply(seq_along(agents),function(a) 
    rbinom(n = simulation.control$n_issues,size = 1,prob = incentive.matrix[a,,t])))
   
  ctbs = incentive.matrix[,,t] * issue.matrix[,,t] * agent.best.guesses  * draws *  {(agent.best.guesses>1)+0}
  orig.contrib.matrix[,,t] <-  as.matrix(ctbs)
  orig.contrib.matrix[,,t][is.na(orig.contrib.matrix[,,t])]<-0
  
  cgr.reallocate.contrib = t(sapply(cgr$agents,function(a) sum(orig.contrib.matrix[a,cgr$issues,t]) * { cgr.payout.estimates[cgr$issues]/sum(cgr.payout.estimates[cgr$issues])}))
  original.contrib = orig.contrib.matrix[cgr$agents,cgr$issues,t]
  
  #calculate capacity for joint action at time t, used as basis for reallocation
  dynamics.tracker$capacity.for.joint.action[t] <- mean(original.contrib>0)
  cgr.contribs = original.contrib * (1-dynamics.tracker$capacity.for.joint.action[t]) + cgr.reallocate.contrib * (dynamics.tracker$capacity.for.joint.action[t])
  #cgr reallocated contribs
  contrib.matrix[,,t] <- orig.contrib.matrix[,,t]
  contrib.matrix[cgr$agents,cgr$issues,t] <- cgr.contribs
  #cgr.contribs = contrib.matrix[cgr$agents,cgr$issues,t]
  payout.matrix[,,t] <- t(replicate(length(agents),colSums({contrib.matrix[,,t] * t(replicate(length(agents),true.values))},na.rm = T)/
                                      colSums(issue.matrix[,,t],na.rm = T))) * issue.matrix[,,t]
  reciprocity.matrix[,,t] = orig.contrib.matrix[,,t] * agent.best.guesses * incentive.matrix[,,t] < payout.matrix[,,t]
  reciprocity.matrix[,,t][is.na( reciprocity.matrix[,,t])]<-F
  contributors.matrix[,,t] <- colSums(orig.contrib.matrix[,,t]>0,na.rm=T)
  if(simulation.control$behavior == 'consistent'){beta = 0.90;alpha1 = 0.05;alpha2 = 0.05}
  if(simulation.control$behavior  == 'contingent'){beta = 0.05;alpha1 = 0.90;alpha2 = 0.05}
  if(simulation.control$behavior == 'conforming'){beta = 0.05;alpha1 = 0.05;alpha2 = 0.90}
  
    #cumulative
    #incentiveMatrix[,t+1]<-incentiveMatrix[,t]*beta+alpha1*(rowSums(cbind(reciprocityMatrix[,1:t]),na.rm = T)/t)+alpha2*{cumsum(contributorsMatrix[1:t])[t]/(t*n_agents)}
    #this time
    recip = Reduce('+',lapply(max(1,t-3):max(1,t),function(x) (reciprocity.matrix[,,x]+0)))/t
    all.contrib = apply(orig.contrib.matrix[,,max(1,t-3):max(1,t)]>0,2,sum) 
    all.contrib = t(replicate(length(agents),all.contrib))
    self.contrib = Reduce('+',lapply(max(1,t-3):max(1,t),function(x) (orig.contrib.matrix[,,x]>0)+0))
    
    all.in.issue = apply(issue.matrix[,,max(1,t-3):max(1,t)]>0,2,sum)
    all.in.issue = t(replicate(length(agents),all.in.issue))
    self.issue = Reduce('+',lapply(max(1,t-3):max(1,t),function(x) (issue.matrix[,,x]>0)+0))
    contrib = (all.contrib - self.contrib) / (all.in.issue - self.issue)
    cgr_payout = apply(contrib.matrix[,,t] * incentive.matrix[,,t] * t(replicate(length(agents),0 + {seq(simulation.control$n_issues) %in% cgr$issues})) ,1,sum) 
    #calculate shared motivation, used as basis for probabilistic issue uptake
    dynamics.tracker$shared.motivation[t] <- mean(issue.matrix[cgr$agents,cgr$issues,t])
    if(t<simulation.control$t){
    incentive.matrix[,,t+1]<-   incentive.matrix[,,t] * beta + alpha1 * recip + alpha2 * contrib
    issue.matrix[,,t+1] <- issue.matrix[,,t]
    issue.matrix[,cgr$issues,t+1] <- {(issue.matrix[,cgr$issues,t] + 
      {(matrix(rbinom(length(cgr$issues) * length(agents),1,
                      dynamics.tracker$shared.motivation[t]),
           ncol = length(cgr$issues)) | 
      {(issue.matrix[,cgr$issues,t]==1) + 0})+0} * 
      ifelse(1:length(agents) %in% cgr$agents,T,F))>0}+0
    }
  #possible perturbations
  if(t == simulation.control$perturb.time&!is.null(simulation.control$perturb.type)){
    if(simulation.control$perturb.type=='payoff.change'){
      new.true.values = {runif(n = simulation.control$n_issues,min = simulation.control$min.payout,max = simulation.control$max.payout)}
      new.info.pieces = replicate(length(agents),createAgentInformationMatrix(mask = simulation.control$mask),simplify = F)
      true.values[cgr$issues] <- new.true.values[cgr$issues]
    }
    if(simulation.control$perturb.type=='lose.contributor'){
      cgr$dropped.agent = sample(x = cgr$agents,size = 1,prob = rowSums(contrib.matrix[cgr$agents,,t]))
      cgr$agents = cgr$agents[cgr$agents!=cgr$dropped.agent]
    }
    if(simulation.control$perturb.type=='add.agents'){
      new_members = which(seq(length(agents)) %in% cgr$agents)[rbinom(size = 1,n = length(agents)-length(cgr$agents),prob = 0.3)==1]
      cgr$agents = union(cgr$agents,new_members)
    }
  }
}
return(list('payoffs' = true.values,'cgr' = cgr,
            'dynamics' = dynamics.tracker,sim = simulation.control,
            'payout' = payout.matrix,'incentive'=incentive.matrix))
}
