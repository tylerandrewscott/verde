

require(pbapply)
source('code/model/emersonscott_model.R')
require(parallel)
require(doParallel)

cores = detectCores() / 1.25
cl = makeCluster(cores)
reps = 1e3
registerDoParallel(cl)

incentive.set = list(c(0.1,0.5),c(0.3,0.7),c(0.5,0.9),c(0.1,0.9))

clusterEvalQ(cl,expr =  source('code/model/emersonscott_model.R'))


res = foreach(i = 1:reps) %dopar% {
  simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                            incentive.set = sample(incentive.set,1)[[1]],
                            uncertainty = runif(1,min = 0.25,2.25), 
                            n_pieces = 1,
                            min.payout = 0,
                            max.payout = 10,n_issues = 100,
                            number.of.issues.to.join = 2,
                            beta = 0.9,alpha = 0.1, 
                            t = 50,perturb.time = 10,perturb.type = 'payoff.change',
                            CGselector='betweenness',behavior = "consistent")
  tryCatch(EmersonScottModel(simulation.control = simulation.control),error = function(e) NULL)}
saveRDS(res[sapply(res,class)=='list'],paste0('../bucket_mount/verde_scratch/payoff.change.',reps/1e3,'k.RDS'))
rm(res);gc()


# 
res = foreach(i = 1:reps) %dopar% {
  simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                            incentive.set = sample(incentive.set,1)[[1]],
                            uncertainty = runif(1,min = 0.25,2.25),
                            n_pieces = 1,
                            min.payout = 0,
                            max.payout = 10,n_issues = 100,
                            number.of.issues.to.join = 2,
                            beta = 0.9,alpha = 0.1,
                            t = 50,perturb.time = 10,perturb.type = 'lose.contributor',
                            CGselector='betweenness',behavior = "consistent")
  tryCatch(EmersonScottModel(simulation.control = simulation.control),error = function(e) NULL)}

saveRDS(res[sapply(res,class)=='list'],paste0('../bucket_mount/verde_scratch/lose.contributor.',reps/1e3,'k.RDS'))
rm(res);gc()
# 
# 
res = foreach(i = 1:reps) %dopar% {
  simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                            incentive.set = sample(incentive.set,1)[[1]],
                            uncertainty = runif(1,min = 0.25,2.25),
                            n_pieces = 1,
                            min.payout = 0,
                            max.payout = 10,n_issues = 100,
                            number.of.issues.to.join = 2,
                            beta = 0.9,alpha = 0.1,
                            t = 50,perturb.time = 10,perturb.type = 'add.agents',
                            CGselector='betweenness',behavior = "consistent")
  tryCatch(EmersonScottModel(simulation.control = simulation.control),error = function(e) NULL)}
saveRDS(res[sapply(res,class)=='list'],paste0('../bucket_mount/verde_scratch/add.agents.',reps/1e3,'k.RDS'))
rm(res);gc()
# 


