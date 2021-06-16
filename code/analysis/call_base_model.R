

require(pbapply)
source('code/model/emersonscott_model.R')
require(parallel)
require(doParallel)

cores = floor(detectCores() / 1.2)
cl = makeCluster(cores)
reps = 1e3
registerDoParallel(cl)

clusterEvalQ(cl,expr =  source('code/model/emersonscott_model.R'))

incentive.set = list(c(0.1,0.5),c(0.3,0.7),c(0.5,0.9),c(0.1,0.9))

res = foreach(i = 1:reps) %do% {
  simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                            incentive.range  = runif(n = 1,min = 0.1,max = 0.9),
                           # skill.set = runif(1,0.1,0.9),
                            #capacity.set = runif(1,0.1,0.9),
                            uncertainty = runif(1,min = 0.25,2.25), 
                            n_pieces = 1,
                            min.payout = 0,
                            max.payout = 10,n_issues = 100,
                            number.of.issues.to.join = 2,
                            beta = 0.9,alpha = 0.1, 
                            t = 50,perturb.time = 15,perturb.type = NULL,
                            CGselector='betweenness',behavior = "consistent")
  tryCatch({EmersonScottModel(simulation.control = simulation.control)},error = function(e) NULL)
}

saveRDS(res,paste0('scratch/test_baserun.',reps,'reps.RDS'))
rm(res);gc()
stopCluster(cl)


