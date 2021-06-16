
run_base = T
run_perturbs = F
require(pbapply)
source('code/model/emersonscott_model.R')
require(parallel)
require(doParallel)

cores = detectCores() /2
cl = makeCluster(cores)
reps = 1e3
registerDoParallel(cl)

if(run_base){
res = foreach(i = 1:reps) %dopar% {
  source('code/model/emersonscott_model.R')
  simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                            motivation.set =sort(runif(2,0.1,0.9)),
                            skill.set = sort(runif(2,0.1,0.9)),
                            capacity.set = sort(runif(2,0.1,0.9)),
                            uncertainty = runif(1,min = 0.25,2.25), 
                            n_pieces = 1,
                            min.payout = 0,
                            max.payout = 10,n_issues = 100,
                            number.of.issues.to.join = 2,
                            beta = 0.9,alpha = 0.1, 
                            t = 50,perturb.time = 15,perturb.type = NULL,
                            CGselector='betweenness',behavior = "consistent")
  EmersonScottModel(simulation.control = simulation.control)}
saveRDS(res,paste0('../bucket_mount/verde_scratch/baserun.',reps/1e3,'k.RDS'))
rm(res);gc()
}
stopCluster()


if(run_perturbs){
simulation.control$perturb.type<-'lose.contributor'
res = pblapply(1:reps,function(i){
  EmersonScottModel(simulation.control = simulation.control)
},cl = cores)
saveRDS(res[sapply(res,class)=='list'],paste0('../bucket_mount/verde_scratch/lose.contributor.',reps/1e3,'k.RDS'))
rm(res);gc()

simulation.control$perturb.type<-'add.agents'
res = pblapply(1:reps,function(i){
  EmersonScottModel(simulation.control = simulation.control)
},cl = cores)

saveRDS(res[sapply(res,class)=='list'],paste0('../bucket_mount/verde_scratch/add.agents.',reps/1e3,'k.RDS'))
rm(res);gc()

simulation.control$perturb.type<-'payoff.change'
res = pblapply(1:1000,function(i){
  EmersonScottModel(simulation.control = simulation.control)
},cl = cores)

saveRDS(res[sapply(res,class)=='list'],paste0('../bucket_mount/verde_scratch/payoff.change.',reps/1e3,'k.RDS'))
rm(res);gc()
}


