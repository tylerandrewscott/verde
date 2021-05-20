
packs = c('data.table','tidyverse','lsa','statnet','igraph','multinets','tidygraph','intergraph','tidyr','mlergm')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

incentive.set = list(c(0.1,0.5),c(0.3,0.7),c(0.5,0.9),c(0.1,0.9))
source('code/model/emersonscott_model.R')
simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                          incentive.set = sample(incentive.set,1)[[1]],
                          uncertainty = runif(1,min = 0.25,2.25), 
                          n_pieces = 1,
                          min.payout = 0,
                          max.payout = 10,n_issues = 100,
                          beta = 0.9,alpha = 0.1, 
                          t = 50,perturb.time = 15,perturb.type = NULL,
                          CGselector='betweenness',behavior = "consistent")

res = tryCatch({EmersonScottModel(simulation.control = simulation.control,debug = T)},error = function(e) NULL)


py = data.table(melt(res$payouts))
setnames(py,c('Var1','Var2','Var3'),c('Agent','Issue','Time'))
py = py[,sum(value),by=.(Issue,Time)]
py$cgr_issue = py$Issue %in% res$cgr$issues
py$mult = res$mults[py$Issue]
ggplot(py,aes(x = Time,y = V1,group = Issue,col = mult)) + 
  geom_path(lwd = 0.2) + scale_color_viridis_c(name = 'PG multiplier') + 
  theme_bw() + theme(legend.position = c(0.8,0.5)) + 
  ggtitle('Public goods created over time') + 
  ylab("Public goods generated")


g1 = ggplot(py,aes(x = Time,y = V1,group = Issue,col = as.factor(cgr_issue),alpha = cgr_issue,size = cgr_issue)) + 
  geom_path() + scale_color_manual(values=c('grey80','black'),name = 'issue in CGR')+
  theme_bw() + theme(legend.position = c(0.8,0.5)) + 
  scale_alpha_manual(values= c(0.6,0.9),name = 'issue in CGR')+
  scale_size_manual(values =c(0.2,0.4),name = 'issue in CGR')+
  ggtitle('Individual public goods') + 
  scale_x_continuous(limits=c(1,50),expand = c(0,0))+
  ylab("Public goods generated")

g2 = ggplot(py[,sum(V1),by=.(Time,cgr_issue)][cgr_issue==T,],
            aes(x = Time,y = V1,group = cgr_issue,ymin = 0,ymax = V1,fill = cgr_issue)) + 
  geom_path(lwd = 0.2) + 
  geom_ribbon() + 
  scale_color_manual(values=c('grey80','black'),name = 'issue in CGR')+
  scale_fill_manual(values = c('grey20'),name = '',labels =c('CGR issues'))+
  scale_x_continuous(limits=c(1,50),expand = c(0,0))+
  theme_bw() + theme(legend.position = c(0.8,0.5)) + 
  ggtitle('Total public goods (CGR issues only)') + 
  ylab("Public goods generated") + theme(legend.title = element_blank())


grid.arrange(g1,g2,ncol = 2,top = 'Public goods generated over time')


###### 3.1 Creating the Policy System ######

require(ggnetwork)


require(doParallel)
seed = 24
  #### actions for 3.1 #######
  agents = createAgents(simulation.control = simulation.control)
  true.values = createIssues(simulation.control = simulation.control)
  agent.motivation = createMotivation(simulation.control = simulation.control,agents = agents)
  
  agent.issue.incidence.matrix = matrix(0+(runif(simulation.control$n_issues*length(agents))>0.75),ncol = simulation.control$n_issues)
  issue.homophily <- lsa::cosine(t(agent.issue.incidence.matrix))
  agent.homophily <- lsa::cosine(agent.issue.incidence.matrix)
  
  seed_network = network(rgraph(length(agents),tprob = 0.04),directed = F)
  sqrt.agent.homophily = sqrt(agent.homophily)
  agent.network = (simulate(seed_network ~ edges + isolates + gwdegree(0.5,fixed = T) + gwesp(0.5,fixed = T) + 
                              edgecov(sqrt.agent.homophily),coef = c(-2,-Inf,-0.75,0.75,1),constraints = ~edges))
  
  anet = data.table(ggnetwork(agent.network))
  hphily = mapply(function(i,j) agent.homophily[i,j],i = anet[,1],j = anet[,2],SIMPLIFY = T)

  anet$hphily <- NA
  anet$iso = ifelse(anet$x==anet$xend&anet$y==anet$yend,T,F)
  anet$hphily[!anet$iso]<-hphily
  
  ggplot(anet,aes(x = x,y=y,xend=xend,yend=yend)) + geom_edges(aes(size = hphily)) + 
    geom_nodes(pch=21,fill = 'white',size = 1.5) + ggtitle('Simulated agent network\nwith homophily based on shared interests')+
    scale_size(range=c(0.2,1.2),name = 'cosine similarity\nbetween agents') + theme_blank() +
    theme(legend.position = 'bottom',legend.direction = 'horizontal')
    

  #a.elist = as.edgelist(agent.network)
  #an = ggnetwork(agent.network)
  
  
  
  
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


  node.df = data.table(NodeID = 1:(length(agents)+simulation.control$n_issues),Node.Type =c(rep('Actor',length(agents)),rep('Issue',simulation.control$n_issues)))
  
  mln = mlnet(complete.edges,node_memb = node.df$Node.Type)
  require(ggthemes)
  
  mln = network(complete.edges,matrix.type = 'edgelist')
  mln %v% 'node_memb' <- c(rep('agent',length(agents)),rep('issue',simulation.control$n_issues))
  cgr = createCGR()
  
  mln %v% 'CGR' <- (1:network.size(mln)) %in% c(cgr$agents,cgr$issues + length(agents))
  
  ggmln = ggnetwork(mln)
  
  
  match(paste(as.character(ggmln$xend),as.character(ggmln$yend)),
        paste(as.character(ggmln$x),as.character(ggmln$y)))
  ggmln[!(duplicate)]
  
  ggplot(ggnetwork(mln),aes(x = x,y = y,xend = xend,yend = yend,col = paste(CGR,node_memb),
                            fill = paste(CGR,node_memb),shape = paste(CGR,node_memb))) + 
    geom_edges(alpha = 0.5,lwd = 0.2,col = 'grey50') + 
    geom_nodes() + 
    theme_blank() +
    scale_shape_manual(values= c(24,21,24,21),labels=c('agent','issue','agent in CGR','issue in CGR')) + 
    scale_fill_manual(values=c('white','white','red','red'),labels=c('agent','issue','agent in CGR','issue in CGR'))+
    scale_colour_manual(values=c('black','black','red','red'),labels=c('agent','issue','agent in CGR','issue in CGR'))+
    theme(legend.position = 'bottom',legend.title = element_blank(),
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
  
  
  
  
  
  
  
  
  
  
    