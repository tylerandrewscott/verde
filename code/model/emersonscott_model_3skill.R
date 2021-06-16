
packs = c('data.table','tidyverse','lsa','statnet','igraph','multinets','tidygraph','intergraph','tidyr')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                          motivation.set = runif(1,0.1,0.9),
                          skill.set = runif(1,0.1,0.9),
                          capacity.set = runif(1,0.1,0.9),
                          uncertainty = runif(1,min = 0.25,2.25), 
                          n_pieces = 1,
                          min.payout = 0,
                          max.payout = 10,n_issues = 100,
                          number.of.issues.to.join = 2,
                          beta = 0.9,alpha = 0.1, 
                          t = 50,perturb.time = 15,perturb.type = NULL,
                          CGselector='betweenness',behavior = "consistent")

###### 3.1 Creating the Policy System ######
setDTthreads(threads = 1)
##### functions for 3.1 ######
createAgents = function(simulation.control) {c(rep('stakeholders',simulation.control$stakeholders),
                                               rep('regulators',simulation.control$regulators),
                                               rep('convenors',simulation.control$convenors))}
createIssues = function(simulation.control) {runif(n = simulation.control$n_issues,min = simulation.control$min.payout,max = simulation.control$max.payout)}

createAttribute = function(simulation.control,agents,skill) {temp = rnorm(n = length(agents),mean = simulation.control[[skill]]);temp[temp<=0.01]<-0.01;temp[temp>=0.99]<-0.99}
#createMotivation = function(simulation.control,agents) {runif(length(agents), simulation.control$motivation.set[1],simulation.control$motivation.set[2])}
#createSkill= function(simulation.control,agents) {runif(length(agents), simulation.control$skill.set[1],simulation.control$skill.set[2])}
#createCapacity = function(simulation.control,agents) {runif(length(agents), simulation.control$capacity.set[1],simulation.control$capacity.set[2])}



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
  
  agent.capacity=createAttribute(simulation.control = simulation.control,agents = agents,skill = 'capacity.set')
    agent.communication = createAttribute(simulation.control = simulation.control,agents = agents,skill = 'skill.set')
    agent.motivation = createAttribute(simulation.control = simulation.control,agents = agents,skill = 'motivation.set')
  #agent.communication = createMotivation(simulation.control = simulation.control,agents = agents)
  #agent.capacity = createMotivation(simulation.control = simulation.control,agents = agents)
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
  
  capacity.matrix = skill.matrix = shared.info.matrix = private.info.matrix = issue.matrix = motivation.matrix = contrib.matrix = orig.contrib.matrix = payout.matrix = orig.payout.matrix = reciprocity.matrix  = array(NA,dim = list(length(agents),length(true.values),simulation.control$t))
  contributors.matrix = array(NA,dim = list(1,length(true.values),simulation.control$t))
  motivation.matrix[,,1] = agent.motivation
  capacity.matrix[,,1] = agent.capacity
  skill.matrix[,,1] = agent.communication
  issue.matrix[,,1] <- agent.issue.incidence.matrix
  
  dynamics.tracker = list('principled.engagement'=NULL,'capacity.for.joint.action'=NULL,'shared.motivation'=NULL)
  info.pieces = replicate(length(agents),createAgentInformationMatrix(payoffs =true.values,ag = agents,
                                                                      u = simulation.control$uncertainty,
                                                                      n_pieces = simulation.control$n_pieces),simplify = F)
  
 for(t in 1:simulation.control$t){
    #print(t)
    ### 3.2.1 Agents share information #### 
    
    

    private.info.matrix[,,t]<- do.call(rbind,info.pieces) 
     
    if(t > 1){
       private.info.matrix[,,t][orig.contrib.matrix[,,t-1]==0]<-NA
      }

    all.private.info = data.table(melt(private.info.matrix))
    setnames(all.private.info,c('Var1','Var2'),c('Agent','Issue'))
    all.private.info = all.private.info[!is.na(value),]
    all.private.info[,r:=rank(-Var3),by=.(Agent,Issue)]
    all.private.info = all.private.info[r %in% 1:3,]
    pim = dcast( all.private.info[,mean(value,na.rm = T),by=.(Agent,Issue)],Agent~Issue,value.var = 'V1')
    pim[,Agent:=NULL]
    shared.info.matrix[,,t] =  as.matrix(pim / sqrt(skill.matrix[,,t]))
    
    cgr.info = pim / sqrt(skill.matrix[,,t])

 all.shared.info = data.table(melt(shared.info.matrix[,,t]))
    setnames(all.shared.info,c('Var1','Var2'),c('Sharing.Agent','Issue'))
    all.shared.info = all.shared.info[,mean(value,na.rm = T),by=.(Sharing.Agent,Issue)]
    agent.edges = data.table(agent.edges)
    agent.alters = lapply(seq_along(agents),function(a) agent.edges[value==1&Var1==a,]$Var2)
    
    #agent network info
    personal.network.info.t = do.call(rbind,lapply(seq(length(agents)),function(x) all.shared.info[Sharing.Agent %in% agent.alters[[x]],mean(V1,na.rm = T),by=.(Issue)]$V1))
    
    #agent.private.info.t = pim#dcast(all.private.info[,mean(value,na.rm = T),by=.(Agent,Issue)],Agent ~ Issue,value.var = 'V1')

    cgr.payout.estimates = colMeans(cgr.info[cgr$agents,])
    info.in.cgr = as.matrix((pim + personal.network.info.t + cgr.payout.estimates)/3)
    info.out.cgr = as.matrix((pim + personal.network.info.t)/2)
  
    
    
    dynamics.tracker$principled.engagement[t] <- sqrt(sum((info.in.cgr[cgr$agents,cgr$issues] - 
                              do.call(rbind,replicate(expr = true.values[cgr$issues],n = length(cgr$agents),simplify =F)))^2) / 
                                {length(cgr$agents)*length(cgr$issues)})
    

    operating.info = do.call(rbind,lapply(seq(length(agents)),function(a) {if(a %in% cgr$agents){info.in.cgr[a,]}else{info.out.cgr[a,]}}))
    agent.best.guesses = operating.info
   # agent.issue.knowledge.pieces = lapply(info.pieces,function(x) colSums(!is.na(x),na.rm = T))
  #  agent.binary.issue.knowledge = lapply(info.pieces,function(x) (colSums(!is.na(x),na.rm = T)>0)+0)
    
    draws = do.call(rbind,lapply(seq_along(agents),function(a) 
      rbinom(n = simulation.control$n_issues,size = 1,prob = motivation.matrix[a,,t])))
  
    ctbs = motivation.matrix[,,t] * issue.matrix[,,t] * agent.best.guesses  * draws
    orig.contrib.matrix[,,t] <-  as.matrix(ctbs)
    orig.contrib.matrix[,,t][is.na(orig.contrib.matrix[,,t])]<-0
    
    cgr.reallocate.contrib = t(sapply(cgr$agents,function(a) sum(orig.contrib.matrix[a,cgr$issues,t]) * { cgr.payout.estimates[cgr$issues]/sum(cgr.payout.estimates[cgr$issues])}))
    original.contrib = orig.contrib.matrix[cgr$agents,cgr$issues,t]
    
    #calculate capacity for joint action at time t, used as basis for reallocation
    cgr.contribs = (original.contrib * (1 - colMeans(capacity.matrix[cgr$agents,cgr$issues,t]))) + cgr.reallocate.contrib *colMeans(capacity.matrix[cgr$agents,cgr$issues,t])
    
     #cgr reallocated contribs
    contrib.matrix[,,t] <- orig.contrib.matrix[,,t]
    contrib.matrix[cgr$agents,cgr$issues,t] <- cgr.contribs
  
    payoff.by.issue.t = colSums(contrib.matrix[,,t] * true.values)
    split.by.issue.t = payoff.by.issue.t/colSums(issue.matrix[,,t])
    
    payout.matrix[,,t] = matrix(rep(split.by.issue.t,each = length(agents)),byrow=F,nrow = length(agents)) * issue.matrix[,,t]
    reciprocity.matrix[,,t] = contrib.matrix[,,t] * agent.best.guesses * motivation.matrix[,,t] < payout.matrix[,,t]
    contributors.matrix[,,t] <- colSums(orig.contrib.matrix[,,t]>0,na.rm=T)
    
    cgr.reciprocity = rowSums(contrib.matrix[cgr$agents,cgr$issues,t] * agent.best.guesses[cgr$agents,cgr$issues] * motivation.matrix[cgr$agents,cgr$issues,t]) < rowSums(payout.matrix[cgr$agents,cgr$issues,t])
    reciprocity.matrix[cgr$agents,cgr$issues,t] <- cgr.reciprocity
    
    cgr.better = rowSums(original.contrib * agent.best.guesses[cgr$agents,cgr$issues] * capacity.matrix[cgr$agents,cgr$issues,t]) < 
      rowSums(contrib.matrix[cgr$agents,cgr$issues,t] * agent.best.guesses[cgr$agents,cgr$issues] * capacity.matrix[cgr$agents,cgr$issues,t])
    #joint.action.matrix[cgr$agents,cgr$issues,t] <- cgr.better
    
    dynamics.tracker$shared.motivation[t] <- mean(reciprocity.matrix[cgr$agents,cgr$issues,t]) #mean(reciprocity.matrix[cgr$agents,cgr$issues,t])
    dynamics.tracker$capacity.for.joint.action[t] <-mean(payout.matrix[cgr$agents,cgr$issues,t] - contrib.matrix[cgr$agents,cgr$issues,t])
  #  pe <- min(sum(motivation.matrix[cgr$agents,cgr$issues,t]*issue.matrix[cgr$agents,cgr$issues,t])/sum(orig.contrib.matrix[cgr$agents,cgr$issues,t]>0),1)
   # dynamics.tracker$principled.engagement[t] <-  pe
    
     recip = Reduce('+',lapply(max(1,t-3):max(1,t),function(x) (reciprocity.matrix[,,x]+0)))/t
    if(t < simulation.control$t){
      issue.matrix[,,t+1] <- issue.matrix[,,t]
      motivation.matrix[,,t+1]<- motivation.matrix[,,t] * simulation.control$beta + simulation.control$alpha * recip
      skill.matrix[,,t+1] <- skill.matrix[,,t] * simulation.control$beta + simulation.control$alpha * dynamics.tracker$principled.engagement[t]
      capacity.matrix[cgr$agents,cgr$issues,t+1] <- capacity.matrix[cgr$agents,cgr$issues,t] * simulation.control$beta + simulation.control$alpha * mean(cgr.better)
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
              'cgr.payout.t' = apply(payout.matrix[,cgr$issues,],3,sum),
              'starting.incentive'=summary(motivation.matrix[,1,1])))}
  if(debug){
    return(list('motivation.matrix' = motivation.matrix,'cgr' = cgr,'mults' = true.values,
  'joint.action' = joint.action.matrix,'principled.engagement' = principled.engagement.matrix,
  'contribs' = contrib.matrix, 'payouts' = payout.matrix, 'sim' = simulation.control))}

}



dynamics.tracker$capacity.for.joint.action
dynamics.tracker$shared.motivation



