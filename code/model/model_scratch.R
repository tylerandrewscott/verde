packs = c('data.table','tidyverse')
need = packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

source('code/model/emersonscott_model.R')
create_agents = function(convenors=0, regulators=0,stakeholders=10){
return(c(rep('convener',convenors),
           rep('regulator',regulators),
           rep('stakeholder',stakeholders)))}

n_issues = 50
n_pieces = 30


cell_interations = 1000
require(data.table)

b1 = data.table(behavior  = c('consistent','contingent','conforming'))
b2 = data.table(variance = 'large',avg = c('low','medium','high'),
                min = c(0.1,0.25,0.4),
           max = c(0.6,0.75,.90))
b3 = data.table(variance='small',avg = c('low','medium','high'),
                min = c(0.3,0.45,0.6),
                max = c(0.4,0.55,0.7))
inputs = full_join(b1,rbind(b2,b3),by = character())


require(tidyverse)
require(pbapply)
results = pblapply(rep(1:nrow(inputs),each = 1000),function(y){
RunSimulation(
behavior = inputs$behavior[y],
bottomIncentive = inputs$min[y],
topIncentive = inputs$max[y],
n_issues = 50,
n_pieces = 30,
stakeholders = 0,regulators = 10,study_phase = FALSE)},cl = 1)






resdt = inputs[rep(1:nrow(inputs),each = 1000),]
resdt$avg <- fct_inorder(resdt$avg)
resdt$behavior <- fct_inorder(resdt$behavior)
resdt$PG_Generated = unlist(lapply(results,function(x) x$PG_Generated))
resdt$Cor_Coef= unlist(lapply(results,function(x) x$rankOrderCor))
resdt$Avg_Motivation= unlist(lapply(results,function(x) x$finalIncentiveAvg))
resdt$Start_Avg_Motivation= unlist(lapply(results,function(x) x$startingIncentiveAvg))
resdt$Var_Motivation= unlist(lapply(results,function(x) x$finalIncentiveVar))


ggplot(resdt,aes(y = Cor_Coef,x = avg)) + 
  scale_y_continuous(limits = c(-0.5,1))+
  facet_grid(variance~ behavior) + geom_jitter()


ggplot(resdt,aes(y = Avg_Motivation,x = avg)) + 
  scale_y_continuous(limits = c(0,1))+
  facet_grid(variance~ behavior) + geom_jitter()


ggplot(resdt,aes(y = Var_Motivation,x = avg)) + 
  scale_y_continuous(limits = c(-.05,.1))+
  facet_grid(variance~ behavior) + geom_jitter()

ggplot(resdt,aes(y = PG_Generated,x = avg)) + 
  facet_grid(variance~ behavior) + geom_jitter(pch = 21,alpha = 0.5) + 
  geom_boxplot()



ggplot(resdt,aes(y = PG_Generated,x = avg)) + facet_grid(variance~ behavior) + geom_jitter()

resdt
rep(1:nrow(inputs),each = 100)
