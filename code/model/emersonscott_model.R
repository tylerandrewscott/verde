
packs = c('data.table','tidyverse','lsa','statnet','igraph','multinets','EnvStats','tidygraph','intergraph','tidyr','reshape2')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                          incentive.range = runif(n = 1,min = 0.1,max = 0.9),
                          capacity.range = runif(n = 1,min = 0.1,max = 0.9),
                          uncertainty = runif(n = 1,min = 0.25,max = 1.75), 
                          n_pieces = 1,beta = 0.9,alpha = 0.1, 
                          min.payout = 0,
                          max.payout = 10,n_issues = 100,
                          number.of.issues.to.join = 2,
                          t = 50,perturb.time = 15,perturb.type = NULL,
                          CGselector='betweenness')

###### 3.1 Creating the Policy System ######
setDTthreads(threads = 1)
##### functions for 3.1 ######
createAgents = function(simulation.control) {c(rep('stakeholders',simulation.control$stakeholders),
                                               rep('regulators',simulation.control$regulators),
                                               rep('convenors',simulation.control$convenors))}
createIssues = function(simulation.control) {runif(n = simulation.control$n_issues,min = simulation.control$min.payout,max = simulation.control$max.payout)}
createAttribute = function(simulation.control,agents,skill) {temp = rnorm(n = length(agents),mean = simulation.control[[skill]],sd = 0.5);temp[temp<=0.01]<-0.01;temp[temp>=0.99]<-0.99;return(temp)}


createAgentInformationMatrix = function(payoffs =true.values,ag = agents,
                                        u = simulation.control$uncertainty,
                                        n_pieces = simulation.control$n_pieces){
  infoFull = do.call(cbind,sapply(seq_along(payoffs),function(x) 
    rnorm(n = simulation.control$n_pieces*length(ag),mean = payoffs[x],sd = u),simplify = F))
  infoFull[infoFull<0] <- 0
  agentInfo = rbind(infoFull[n_pieces,])
  return(agentInfo)}


createCGR = function(select = simulation.control$CGselector,
                     number.to.join = simulation.control$number.of.issues.to.join,
                     issue.net = issue.network,
                     agent.issue.graph= agent.issue.incidence.matrix){
  issue.graph = asIgraph(issue.net)
  if(select=='betweenness'){issue.weight = igraph::betweenness(issue.graph)}
  picked.issue = sample(seq(simulation.control$n_issues),size = 1,prob = issue.weight)
  issues = ego(issue.graph,nodes = picked.issue)
  props = rowSums(agent.issue.graph[,issues[[1]]])/length(issues[[1]])
  props = ifelse(props==0,0.025,ifelse(props==1,0.975,props))
  agents = which(sapply(props,function(x) rbinom(1,1,prob = x))==1)
  list('issues' = sort(issues[[1]]),'agents' = sort(agents))
}

EmersonScottModel = function(simulation.control = NULL,debug = F) {
  #### actions for 3.1 #######
  agents = createAgents(simulation.control = simulation.control)
  true.values = createIssues(simulation.control = simulation.control )
  agent.motivation = createAttribute(simulation.control = simulation.control,agents = agents,skill = 'incentive.range')
  agent.capacity = createAttribute(simulation.control = simulation.control,agents = agents,skill = 'capacity.range')
  
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
  
  cgr = createCGR(select = simulation.control$CGselector,issue.net = issue.network,
                  number.to.join = simulation.control$number.of.issues.to.join,
                  agent.issue.graph= agent.issue.incidence.matrix)
  #cgr = createCGR(select = 'betweenness',number.of.issues.to.join = 2)
  
  ###### 3.2 Running the Policy System ######
  ##### functions for 3.2 ######
  
  #joint.action.matrix = principled.engagement.matrix = 
  agent.best.guesses = shared.info.matrix = private.info.matrix = cgr.wouldpayout.matrix = issue.matrix = incentive.matrix = capacity.matrix = expected.contrib.matrix = contrib.matrix = orig.contrib.matrix = payout.matrix = orig.payout.matrix = reciprocity.matrix  = array(NA,dim = list(length(agents),length(true.values),simulation.control$t))
  contributors.matrix = array(NA,dim = list(1,length(true.values),simulation.control$t))
  incentive.matrix[,,1] = agent.motivation
  capacity.matrix[,,1] = agent.capacity
  issue.matrix[,,1] <- agent.issue.incidence.matrix
#  joint.action.matrix[,,1] = principled.engagement.matrix[,,1] = 
  dynamics.tracker = list('principled.engagement'=NULL,'capacity.for.joint.action'=NULL,'shared.motivation'=NULL,'issue.alignment'=NULL,'cgr.contributions' = NULL)
  info.pieces = replicate(length(agents),createAgentInformationMatrix(payoffs =true.values,ag = agents,
                                                                      u = simulation.control$uncertainty,
                                                                      n_pieces = simulation.control$n_pieces),simplify = F)
  
  #cgr.info.list = list()
  
  
 for(t in 1:simulation.control$t){
   agent_issue_dummy = matrix(NA,nrow = length(agents),ncol = simulation.control$n_issues,byrow = T)
   agent_issue_dummy[,cgr$issues] <- 1
   agent_issue_dummy[,-cgr$issues] <- 0
   agent_issue_dummy[-cgr$agents,]<-0
   true.matrix = matrix(true.values,nrow = length(agents),ncol = simulation.control$n_issues,byrow = T)
    #print(t)
    ### 3.2.1 Agents share information #### 

    private.info.matrix[,,t]<- do.call(rbind,info.pieces) 
    
    ### if new informaiton is passed out at every time t, this ensures that new information only goes to people in the group
     # if(t > 1){
     #   private.info.matrix[,,t][orig.contrib.matrix[,,t-1]==0]<-NA
     #  }

    all.private.info = data.table(reshape2::melt(private.info.matrix))
    setnames(all.private.info,c('Var1','Var2'),c('Agent','Issue'))
    all.private.info = all.private.info[!is.na(value),]
    all.private.info[,r:=rank(-Var3),by=.(Agent,Issue)]
    all.private.info = all.private.info[r %in% 1:3,]
    
    
    pim = dcast( all.private.info[,mean(value,na.rm = T),by=.(Agent,Issue)],Agent~Issue,value.var = 'V1')
    pim = data.table(pim)
    pim[,Agent:=NULL]

    
    shared.info.matrix[,,t] =  as.matrix(pim / {1/sqrt(incentive.matrix[,,t])})
 
    all.shared.info = data.table(reshape2::melt(shared.info.matrix[,,t]))
    setnames(all.shared.info,c('Var1','Var2'),c('Sharing.Agent','Issue'))
    all.shared.info = all.shared.info[,mean(value,na.rm = T),by=.(Sharing.Agent,Issue)]
    agent.edges = data.table(agent.edges)
    agent.alters = lapply(seq_along(agents),function(a) agent.edges[value==1&Var1==a,]$Var2)
    
    #agent network info
    personal.network.info.t = do.call(rbind,lapply(seq(length(agents)),function(x) all.shared.info[Sharing.Agent %in% agent.alters[[x]],mean(V1,na.rm = T),by=.(Issue)]$V1))
    
    #agent.private.info.t = pim#dcast(all.private.info[,mean(value,na.rm = T),by=.(Agent,Issue)],Agent ~ Issue,value.var = 'V1')

    cgr.payout.estimates = colMeans(shared.info.matrix[cgr$agents,,t])
    info.in.cgr = as.matrix((pim + personal.network.info.t + do.call(rbind,replicate(length(agents),cgr.payout.estimates,simplify = F)))/3)
    info.out.cgr = as.matrix((pim + personal.network.info.t)/2)
    operating.info = do.call(rbind,lapply(seq(length(agents)),function(a) {if(a %in% cgr$agents){info.in.cgr[a,]}else{info.out.cgr[a,]}}))
    agent.best.guesses[,,t] = operating.info
   # agent.issue.knowledge.pieces = lapply(info.pieces,function(x) colSums(!is.na(x),na.rm = T))
  #  agent.binary.issue.knowledge = lapply(info.pieces,function(x) (colSums(!is.na(x),na.rm = T)>0)+0)
    
    draws = do.call(rbind,lapply(seq_along(agents),function(a) 
      rbinom(n = simulation.control$n_issues,size = 1,prob = incentive.matrix[a,,t])))
  
    #ctbs = incentive.matrix[,,t] * issue.matrix[,,t] * agent.best.guesses[,,t]  
    ctbs = incentive.matrix[,,t] * issue.matrix[,,t] * agent.best.guesses[,,t]  
    expected.contrib.matrix[,,t] <- ctbs
  
    orig.contrib.matrix[,,t] <-  as.matrix(ctbs * draws)
    orig.contrib.matrix[,,t][is.na(orig.contrib.matrix[,,t])]<-0
    
  
    cgr.reallocate.contrib = t(sapply(cgr$agents,function(a) sum(orig.contrib.matrix[a,cgr$issues,t]) * { cgr.payout.estimates[cgr$issues]/sum(cgr.payout.estimates[cgr$issues])}))
    original.contrib = orig.contrib.matrix[cgr$agents,cgr$issues,t]
    
    
    
    #calculate capacity for joint action at time t, used as basis for reallocation
    cgr.contribs = (original.contrib * (1 - colMeans(incentive.matrix[cgr$agents,cgr$issues,t]))) + cgr.reallocate.contrib *colMeans(incentive.matrix[cgr$agents,cgr$issues,t])
    
    
    
    #total.cgr.contrib.expectation =  do.call(rbind,lapply(seq_along(agents), function(x) 
    #  colSums(issue.matrix[,,t] * agent_issue_dummy * orig.contrib.matrix[,,t] *  incentive.matrix[,,t]  * t(replicate(length(agents),agent.best.guesses[x,,t]))))) 
    
    #individual.cgr.payout.expectation = total.cgr.contrib.expectation/length(cgr$agents)
    

    
    #cgr reallocated contribs
    contrib.matrix[,,t] <- orig.contrib.matrix[,,t]
    contrib.matrix[cgr$agents,cgr$issues,t] <- cgr.contribs
  
    payoff.by.issue.t = colSums(contrib.matrix[,,t] * true.values)
    split.by.issue.t = payoff.by.issue.t/colSums(issue.matrix[,,t])
  
    payout.matrix[,,t] = matrix(rep(split.by.issue.t,each = length(agents)),byrow=F,nrow = length(agents)) * issue.matrix[,,t]
    reciprocity.matrix[,,t] = contrib.matrix[,,t] * agent.best.guesses[,,t] * incentive.matrix[,,t] < payout.matrix[,,t]
    
   # reciprocity.matrix[cgr$agents,cgr$issues,t] <- individual.cgr.payout.expectation[cgr$agents,cgr$issues] < payout.matrix[cgr$agents,cgr$issues,t]
    
    contributors.matrix[,,t] <- colSums(orig.contrib.matrix[,,t]>0,na.rm=T)
    
    #cgr.reciprocity = rowSums(contrib.matrix[cgr$agents,cgr$issues,t] * agent.best.guesses[cgr$agents,cgr$issues,t] * incentive.matrix[cgr$agents,cgr$issues,t]) < rowSums(payout.matrix[cgr$agents,cgr$issues,t])
    #reciprocity.matrix[cgr$agents,cgr$issues,t] <- cgr.reciprocity
    
    cgr.better = rowSums(original.contrib * agent.best.guesses[cgr$agents,cgr$issues,t] * incentive.matrix[cgr$agents,cgr$issues,t]) < 
      rowSums(contrib.matrix[cgr$agents,cgr$issues,t] * agent.best.guesses[cgr$agents,cgr$issues,t] * incentive.matrix[cgr$agents,cgr$issues,t])
    #joint.action.matrix[cgr$agents,cgr$issues,t] <- cgr.better
    
   
    #principled engagement: delta between group information and network information

    cgr.error = sum((agent.best.guesses[cgr$agents,cgr$issues,t] - true.matrix[cgr$agents,cgr$issues])^2) / length(true.matrix[cgr$agents,cgr$issues])

    #network.error = sum((agent.best.guesses[,cgr$issues,t] - true.matrix)^2) / length(true.matrix)
    #network.error = sum((info.out.cgr[cgr$agents,cgr$issues] - true.matrix)^2) / length(true.matrix)
    dynamics.tracker$principled.engagement[t] <- 1/cgr.error #- network.error
  
    ### issue alignment is cosine similarity of 1/0 issue matrix
    issue_cosim = cosine(t(issue.matrix[cgr$agents,cgr$issues,t]))
    dynamics.tracker$issue.alignment[t] <-  mean(issue_cosim[upper.tri(issue_cosim)])
    
    #shared motivation: variance in total contributions

    #dynamics.tracker$shared.motivation[t] <- sum(rowSums(original.contrib)>0)/length(cgr$agents)
  
    #mean(c(contrib.matrix[cgr$agents,cgr$issues,t]))
    #contrib.error = sum((contrib.matrix[cgr$agents,cgr$issues,t] - expected.contrib.matrix[cgr$agents,cgr$issues,t])^2) / length(contrib.matrix[cgr$agents,cgr$issues,t])
    
    #network.error = sum((agent.best.guesses[,cgr$issues,t] - true.matrix)^2) / length(true.matrix)
    #network.error = sum((info.out.cgr[cgr$agents,cgr$issues] - true.matrix)^2) / length(true.matrix)
    dynamics.tracker$shared.motivation[t] <- geoMean(rowSums(incentive.matrix[cgr$agents,cgr$issues,t]))
    dynamics.tracker$cgr.contributions[t] <- sum(contrib.matrix[cgr$agents,cgr$issues,t])
    #capacity for joint action is difference between payout generated with and without reallocation  -- i.e., value creation
    no.reallocate.payout = sum(original.contrib * agent.best.guesses[cgr$agents,cgr$issues,t] * incentive.matrix[cgr$agents,cgr$issues,t]) 
    cgr.reallocate.payout = sum(contrib.matrix[cgr$agents,cgr$issues,t] * agent.best.guesses[cgr$agents,cgr$issues,t] * incentive.matrix[cgr$agents,cgr$issues,t])
    cgr.reallocate.payout - no.reallocate.payout
    #dynamics.tracker$capacity.for.joint.action[t] <- sum((contrib.matrix[cgr$agents,cgr$issues,t] - original.contrib)^2)/length(original.contrib)
    dynamics.tracker$capacity.for.joint.action[t] <- sum(c(original.contrib))/length(cgr$agents)
   # pe <- min(sum(incentive.matrix[cgr$agents,cgr$issues,t]*issue.matrix[cgr$agents,cgr$issues,t])/sum(orig.contrib.matrix[cgr$agents,cgr$issues,t]>0),1)
  
     recip = Reduce('+',lapply(max(1,t-3):max(1,t),function(x) (reciprocity.matrix[,,x]+0)))/t
   
    if(t < simulation.control$t){
      issue.matrix[,,t+1] <- issue.matrix[,,t]
      incentive.matrix[,,t+1]<- incentive.matrix[,,t] * simulation.control$beta + simulation.control$alpha * recip
    #  incentive.matrix[,,t+1]<- incentive.matrix[,,t] * simulation.control$beta + simulation.control$alpha * recip
      #principled.engagement.matrix[,,t+1] <- principled.engagement.matrix[,,t] * simulation.control$beta + simulation.control$alpha * pe
      #joint.action.matrix[cgr$agents,cgr$issues,t+1] <- joint.action.matrix[cgr$agents,cgr$issues,t] * simulation.control$beta + simulation.control$alpha * mean(cgr.better)
    }

    if(t == simulation.control$perturb.time&!is.null(simulation.control$perturb.type)){
      if(simulation.control$perturb.type=='payoff.change'){
        new.true.values = createIssues(simulation.control = simulation.control)
        new.info.pieces = replicate(length(agents),createAgentInformationMatrix(payoffs =true.values,ag = agents,
                                                                                u = simulation.control$uncertainty,
                                                                                n_pieces = simulation.control$n_pieces),simplify = F)
        true.values[cgr$issues] <- new.true.values[cgr$issues]
      }
      if(simulation.control$perturb.type=='lose.contributor'){
        cgr$dropped.agent = sample(x = cgr$agents,size = 1,prob = rowSums(contrib.matrix[cgr$agents,,t]))
        cgr$agents = cgr$agents[cgr$agents!=cgr$dropped.agent]
      }
      if(simulation.control$perturb.type=='add.agents'){
        cgr$new_members = which(!seq(length(agents)) %in% cgr$agents)[rbinom(size = 1,n = length(which(!seq(length(agents)) %in% cgr$agents)),prob = 0.25)==1]
        cgr$agents = union(cgr$agents,cgr$new_members)
      }
    }
  }
  if(!debug){
  return(list('payoffs' = true.values,'cgr' = cgr,
              'dynamics' = dynamics.tracker,
              'issue_sd' = summary(apply(private.info.matrix[,,1],2,sd)[cgr$issues]),
              sim = simulation.control,
              incentive.mean = simulation.control$incentive.range,
              capacity.mean = simulation.control$capacity.range,
              'cgr.payout.t' = apply(payout.matrix[,cgr$issues,],3,sum),
              'starting.incentive'=summary(incentive.matrix[,1,1])))}
  if(debug){
    return(list('incentive.matrix' = incentive.matrix,'cgr' = cgr,'mults' = true.values,
  #'joint.action' = joint.action.matrix,'principled.engagement' = principled.engagement.matrix,
  'contribs' = contrib.matrix, 'payouts' = payout.matrix, 'sim' = simulation.control))}

}





