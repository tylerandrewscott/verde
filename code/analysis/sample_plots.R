
packs = c('data.table','tidyverse','gridExtra','lsa','statnet','igraph','multinets','tidygraph','intergraph','tidyr')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

incentive.set = list(c(0.1,0.5),c(0.3,0.7),c(0.5,0.9),c(0.1,0.9))
source('code/model/emersonscott_model.R')
simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                          incentive.range = runif(n = 1,min = 0.1,max = 0.9),
                          capacity.range = runif(n = 1,min = 0.1,max = 0.9),
                          # skill.set = runif(1,0.1,0.9),
                          #capacity.set = runif(1,0.1,0.9),
                          uncertainty = runif(1,min = 0.25,2.25), 
                          n_pieces = 1,
                          min.payout = 0,
                          max.payout = 10,n_issues = 100,
                          number.of.issues.to.join = 2,
                          beta = 0.8,alpha = 0.2, 
                          t = 50,perturb.time = 15,perturb.type = NULL,
                          CGselector='betweenness',behavior = "consistent")
res = tryCatch({EmersonScottModel(simulation.control = simulation.control,debug = F)},error = function(e) NULL)



res$dynamics$principled.engagement
res$dynamics$shared.motivation
res$dynamics$capacity.for.joint.action



res$dynamics$cgr.contributions

par(mfrow = c(2,2))
plot(res$dynamics$capacity.for.joint.action)
plot(res$dynamics$shared.motivation,type = 'lines')
plot(res$dynamics$principled.engagement,type = 'lines')

dyn = res$dynamics
dyn_all = as.data.frame(dyn)
dyn_all$t =  seq(nrow(dyn_all))
dyn_all$cgr.payout = res$cgr.payout.t

dyn_all$total.payout <- sum(res$cgr.payout.t)
dyn_all$iter = 1

require(plotly)
require(viridis)


axis <- list(
  range = c(0, 1)
)

scene <- list(xaxis = axis,yaxis = axis, zaxis = axis)


fig2 <- plot_ly(dyn_all[dyn_all$iter%in%1,], 
                x = ~principled.engagement, z = ~capacity.for.joint.action, color = ~t,
                y = ~shared.motivation, type = 'scatter3d',mode = 'lines',showlegend = F)
fig2 %>% layout(#scene = scene,
  showlegend=T,legend = list(x = 0.5,y =0.1,z = 0.9),margin = list(l = 0,r = -10,t = 0,b = 0))


plot(res$dynamics$principled.engagement)
plot(res$dynamics$shared.motivation)
plot(res$dynamics$capacity.for.joint.action)

plot(res$dynamics$shared.motivation)

py = data.table(melt(res$payouts))
setnames(py,c('Var1','Var2','Var3'),c('Agent','Issue','Time'))
py = py[,sum(value),by=.(Issue,Time)]
py$cgr_issue = py$Issue %in% res$cgr$issues
py$mult = res$mults[py$Issue]

# figure 2
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

#figure 5
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
  
  anet$head <-  anet$vertex.names[match(anet$xend,anet$x)]
  
  anet$hphily = mapply(function(i,j) agent.homophily[i,j],i = anet$vertex.names,j = anet$head,SIMPLIFY = T)
  
  
  agent.homophily[anet$vertex.names,anet$head]
  
  
  hphily = mapply(function(i,j) agent.homophily[i,j],i = anet[,1],j = anet[,2],SIMPLIFY = T)

  anet$hphily <- NA
  anet$iso = ifelse(anet$x==anet$xend&anet$y==anet$yend,T,F)
  anet$hphily[!anet$iso]<-hphily
  
  #figure 1
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

  #figure 3
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
  
  

  
  mults = res$payoffs
  issues = res$cgr$issues
  agents = res$cgr$agents
  info.sd = res$issue_sd[['Mean']]
  dt = data.table(total.payoff = sum(res$cgr.payout.t),
                  final.payoff = res$cgr.payout.t[[50]],
                  issue.sd = unlist(info.sd),
                  med.issue = median(mults[issues]),
                  min.issue = min(mults[issues]),
                  max.issue = max(mults[issues]),
                  mean.issue = mean(mults[issues]),
                  mean.incentive =      mean(res$starting.incentive),
                  med.incentive =      median(res$starting.incentive),
                  min.incentive =      min(res$starting.incentive),
                  max.incentive =      max(res$starting.incentive),
                  agents = length(agents),issues =  length(issues),id = 1)
  
  
  
  

  
  
    