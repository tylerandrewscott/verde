



require(pbapply)

res = pblapply(1:12,function(i){
  simulation.control = list(stakeholders = 50,regulators = 0, convenors = 0 ,
                            bottom.incentive = runif(1,0.1,0.5),
                            top.incentive = runif(1,0.5,0.9),
                            uncertainty = runif(1,min = 0.25,1.75), 
                            n_pieces = 1,
                            min.payout = 0,
                            max.payout = 10,n_issues = 100,
                            number.of.issues.to.join = 2,
                            t = 50,perturb.time = 15,perturb.type = NULL,
                            CGselector='betweenness',behavior = "consistent")
  EmersonScottModel(simulation.control = simulation.control)
},cl = 4)



summary(c(test$incentive))
which(test$incentive>1,arr.ind = T)
dim(test$incentive)

which(test$incentive==1)
summary(c(test$incentive)>1)
