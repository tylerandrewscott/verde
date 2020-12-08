
packs = c('data.table','tidyverse','lsa','statnet','igraph','multinets','tidygraph','intergraph','tidyr')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

storage = list()

simulation.control = list(n_issues = 100,n_pieces = 30,stakeholders = 50,
                          #ranges from 0 (no weight on group estimate) to 1 (full weight on group estimate + 0 distortion on group estimate)
                          principled.engagement = 1,
                          #ranges from 0 (no reallocation) to 1 (full reallocation)
                          capacity.for.joint.action = 1,
                          #shared motivation ranges from 0 (no sharing) to 1 (full sharing)
                          shared.motivation = 1,
                          #types of perturbations
                          perturb.type = NULL,
                          perturb.time = 15,
                          uncertainty = 1,
                          regulators = 0, convenors = 0,min_keep = 30,max_keep = 30,
                          max.payout = 10,min.payout = 0,agent.issue.tie.thresh = 0.9,
                          bottom.incentive = 0.3,top.incentive = 0.7,t = 50,
                          behavior = 'consistent',
                          insertCG = T,CGselector = 'central')

###### 3.1 Creating the Policy System ######

##### functions for 3.1 ######
createAgents = function(control) {c(rep('stakeholders',simulation.control$stakeholders),
                                    rep('regulators',simulation.control$regulators),
                                    rep('convenors',simulation.control$convenors))}
createIssues = function(control) {runif(n = simulation.control$n_issues,min = simulation.control$min.payout,max = simulation.control$max.payout)}
createMotivation = function(control) {runif(length(agents), simulation.control$bottom.incentive,simulation.control$top.incentive)}

createAgentInformationMatrix = function(payoffs =true.values,u = simulation.control$uncertainty,
                                        min_keep = simulation.control$min_keep,
                                        max_keep = simulation.control$max_keep,
                                        n_pieces = simulation.control$n_pieces){
  if(min_keep==max_keep){num_info_pieces_retained = max_keep}
  if(min_keep<max_keep){num_info_pieces_retained = sample(min_keep:max_keep,1)}
  #simulate info 
  infoFull = do.call(cbind,sapply(seq_along(payoffs),function(x) rnorm(n = simulation.control$n_pieces,mean = payoffs[x],sd = u),simplify = F))
  infoFull[ infoFull < 0 ] <- 0 
  agentInfo = rbind(infoFull[1,])
  return(agentInfo)}


createCGR = function(select = NULL,number.of.issues.to.join = 2){
  issue.graph = intergraph::asIgraph(issue.network)
  if(select=='betweenness'){issue.weight = igraph::betweenness(issue.graph)}
  picked.issue = sample(seq(simulation.control$n_issues),size = 1,prob = issue.weight)
  issues = ego(issue.graph,nodes = picked.issue)
  agents = which(rowSums(agent.issue.incidence.matrix[,issues[[1]]])>=number.of.issues.to.join)
  list('issues' = sort(issues[[1]]),'agents' = sort(agents))
}

require(ggnetwork)


require(doParallel)
seed = 24
  #### actions for 3.1 #######
  agents = createAgents()
  true.values = createIssues()
  agent.motivation = createMotivation()
  
  agent.issue.incidence.matrix = matrix(0+(runif(simulation.control$n_issues*length(agents))>0.75),ncol = simulation.control$n_issues)
  issue.homophily <- lsa::cosine(t(agent.issue.incidence.matrix))
  agent.homophily <- lsa::cosine(agent.issue.incidence.matrix)
  
  seed_network = network(rgraph(length(agents),tprob = 0.04),directed = F)
  sqrt.agent.homophily = sqrt(agent.homophily)
  agent.network = (simulate(seed_network ~ edges + isolates + gwdegree(0.5,fixed = T) + gwesp(0.5,fixed = T) + 
                              edgecov(sqrt.agent.homophily),coef = c(-2,-Inf,-0.75,0.75,1),constraints = ~edges))
  
  anet = data.table(ggnetwork(agent.network))
  hphily = mapply(function(i,j) agent.homophily[i,j],i = a.elist[,1],j = a.elist[,2],SIMPLIFY = T)

  anet$hphily <- NA
  anet$iso = ifelse(anet$x==anet$xend&anet$y==anet$yend,T,F)
  anet$hphily[!anet$iso]<-hphily
  
  ggplot(anet,aes(x = x,y=y,xend=xend,yend=yend)) + geom_edges(aes(size = hphily)) + 
    geom_nodes(pch=21,fill = 'white',size = 1.5) + ggtitle('Simulated agent network\nwith homophily based on shared interests')+
    scale_size(range=c(0.2,1.2),name = 'cosine similarity\nbetween agents') + theme_blank() +
    theme(legend.position = 'bottom',legend.direction = 'horizontal')
    
  
  
  anet[!(x==xend&y==yend),]$isonode
  hphily
  
  
  
  
  a.elist = as.edgelist(agent.network)
  
  
  
  dim(ggnetwork(agent.network))
  a.elist[,1],a.elist[,2]
  apply(a.elist,1,function(x) agent.homophily[a.elist[x,1],a.elist[x,2]])
  
  agent.homophily[as.edgelist(agent.network)[,1],as.edgelist(agent.network)[,2]]
  
  an = ggnetwork(agent.network)
  
  
  
  
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

  require(mlergm)
  node.df = data.table(NodeID = 1:(length(agents)+simulation.control$n_issues),Node.Type =c(rep('Actor',length(agents)),rep('Issue',simulation.control$n_issues)))
  
  mln = mlnet(complete.edges,node_memb = node.df$Node.Type)
  require(ggthemes)
  ggplot(ggnetwork(mln),aes(x = x,y = y,xend = xend,yend = yend)) + 
    geom_edges(alpha = 0.5,lwd = 0.2,col = 'grey50') + 
    geom_nodes(aes(shape = node_memb),fill = 'white') + 
    theme_blank() +
    scale_shape_manual(values= c(24,21)) + 
    theme(legend.position = c(0.75,0.1),legend.title = element_blank(),
          text = element_text(size = 12),legend.direction = 'horizontal') + 
    ggtitle('Sample multinode network with agents and issues')
  
  

plot(mln)
  dim(complete.edges)
  install.packages('mlergm')# 
   ig = graph_from_data_frame(data.frame(complete.edges[,c(1:2)]),directed = F,vertices = data.table(vid = 1:(length(agents)+length(true.values)),type = c(rep(FALSE,length(agents)),rep(TRUE,length(true.values)))))
  # 
  ig = igraph::simplify(ig)
  # # Set the layout coordinates
   l <- layout_multilevel(ig, layout = layout_with_kk)
  # # Set different colors and shapes for each level vertices
   ig <- set_color_multilevel(ig,color.true = 'black',color.false = 'red')
   ig <- set_shape_multilevel(ig)

  plot(extract_lowlevel(ig))
   
  # # Plot
   plot(ig, layout = l, vertex.size = 5, vertex.label = NA)
  
  require(networkD3)
  
  
  
  complete.edges$Var1 = complete.edges$Var1 - 1
  complete.edges$Var2 = complete.edges$Var2 - 1
  

node.df
  forceNetwork(Links = complete.edges, Nodes = node.df,
               Source = "Var1", Target = "Var2",
               Value = "value", NodeID = "NodeID",
               Group = "Node.Type", opacity = 0.8)
  
  
  
  
  
  
  
  
  
  
    