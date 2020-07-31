create_agent_information_matrix = function(values,min_keep = min_keep,max_keep = max_keep){
  if(min_keep==max_keep){num_info_pieces_retained = max_keep}
  if(min_keep<max_keep){num_info_pieces_retained = sample(min_keep:max_keep,1)}
  info_mask = matrix(sample(c(rep(NA,length(values) - num_info_pieces_retained),rep(1,num_info_pieces_retained))),ncol = ncol(values),byrow = F)
  agent_info = values * info_mask
  return(agent_info)}

agentGive = function(incentive,estimate){
  rbinom(n=1,size = 1,prob = incentive) * estimate}

collectiveContributions = function(agent_incentives,agent_estimates){
  sapply(seq_along(agent_incentives),function(x){agentGive(incentive= agent_incentives[x],
                                                           estimate = agent_estimates[x])})}

RunSimulation = function(behavior = 'consistent',n_issues = 50,n_pieces = 30,bottomIncentive = 0.3,topIncentive = 0.4,keep_pieces_min  = 30,keep_pieces_max = 30,
                         study_budget = 5,study_phase = FALSE,...){
  #behavior = 'consistent';n_issues = 50;n_pieces = 30;bottomIncentive = 0.3;topIncentive = 0.4;keep_pieces_min  = 30;keep_pieces_max = 30; study_budget = 5
  trueValues = runif(n = n_issues,min = 0,max = 10)
  n_agents = length(create_agents())
  agentIncentive = runif(n_agents,min = bottomIncentive,max = topIncentive)
  infoFull = do.call(cbind,sapply(seq_along(trueValues),function(x) rnorm(n = n_pieces,mean = trueValues[x],sd = 0.1),simplify = F))
  infoPieces = replicate(n_agents,create_agent_information_matrix(values = infoFull,min_keep = keep_pieces_min,max_keep =keep_pieces_max),simplify = F)
  agentInfoDistortion = agentIncentive
  #CI_i = I_i * 1/sqrt(S_i)
  agentInfoShared = mapply(function(info,distortion) info * distortion,info = agentInfo,distortion = agentInfoDistortion,SIMPLIFY = F)
  groupEstimates = rowMeans(simplify2array(agentInfoShared),na.rm = T)
  agentInterpretationOfOthers = lapply(seq_along(agentInfoShared),function(i) {simplify2array(agentInfoShared[-i]) * sqrt(agentInfoDistortion[i])})
  
  agentLearnedInfo = lapply(seq_along(agentInterpretationOfOthers),function(i){
    colMeans(rbind(agentInfo[[i]],
                   do.call(rbind,lapply(seq(dim(agentInterpretationOfOthers[[i]])[3]),function(z) agentInterpretationOfOthers[[i]][,,z])))
             ,na.rm = T)
  })
  
  agentIssueKnowledgePieces = lapply(agentInfo,function(x) colSums(!is.na(x),na.rm = T))
  agentBinaryIssueKnowledge = lapply(agentInfo,function(x) (colSums(!is.na(x),na.rm = T)>0)+0)
  
  issuesKept = {colSums(do.call(rbind,agentBinaryIssueKnowledge)) / n_agents} > 0.5
  issuesUncertainty = apply(do.call(rbind,agentInfoShared),2,var,na.rm = T)
  
  if(study_phase){
  moreStudy = rank(issuesUncertainty)<=5 & issuesKept
  studyInfo = lapply(seq_along(moreStudy),function(x) if(moreStudy[x]){rnorm(mean = trueValues[x],n=study_budget)}else{NA})
  studyInfoMatrix = do.call(cbind,studyInfo)
  agentPostStudyLearnedInfo = lapply(seq_along(agentInterpretationOfOthers),function(i){
    colMeans(rbind(agentInfo[[i]],studyInfoMatrix,
                   do.call(rbind,lapply(seq(dim(agentInterpretationOfOthers[[i]])[3]),function(z) agentInterpretationOfOthers[[i]][,,z])))
             ,na.rm = T)})
  postSharingEstimates = (do.call(rbind,agentPostStudyLearnedInfo))
  }
  if(!study_phase){postSharingEstimates = (do.call(rbind,agentLearnedInfo ))}
  
  individualRanks = apply(postSharingEstimates,1,function(x) rank(-x))
  groupRank = rowMeans(individualRanks)
  
  estimatesRankedIssues = postSharingEstimates[,order(groupRank)]
  trueValuesRanked = trueValues[order(groupRank)]
  startingIncentive = agentIncentive
  
  incentiveMatrix = contribMatrix = payoutMatrix = reciprocityMatrix  = matrix(NA,ncol = length(trueValuesRanked),nrow = n_agents)
  contributorsMatrix = matrix(NA,ncol = length(trueValuesRanked))
  incentiveMatrix[,1] <- startingIncentive 
  for(t in 1:ncol(estimatesRankedIssues)){
    if(all(estimatesRankedIssues[,t]<1)){break}
    if(t>1 & all(contributorsMatrix[t-(1:min(3,t))]==0)){break}
    contribMatrix[,t] <- collectiveContributions(incentiveMatrix[,t],estimatesRankedIssues[,t])
    payoutMatrix[,t] = sum(contribMatrix[,t]) * trueValuesRanked[t] / n_agents
    reciprocityMatrix[,t] = contribMatrix[,t] * incentiveMatrix[,t] < payoutMatrix[,t]
    contributorsMatrix[t] <- colSums(contribMatrix>0)[t]
    if(behavior == 'consistent'){beta = 0.90;alpha1 = 0.05;alpha2 = 0.05}
    if(behavior == 'contingent'){beta = 0.05;alpha1 = 0.90;alpha2 = 0.05}
    if(behavior == 'conforming'){beta = 0.05;alpha1 = 0.05;alpha2 = 0.90}
    if(t<ncol(estimatesRankedIssues)){
      #cumulative
      #incentiveMatrix[,t+1]<-incentiveMatrix[,t]*beta+alpha1*(rowSums(cbind(reciprocityMatrix[,1:t]),na.rm = T)/t)+alpha2*{cumsum(contributorsMatrix[1:t])[t]/(t*n_agents)}
      #this time
      incentiveMatrix[,t+1]<-incentiveMatrix[,t]*beta+alpha1*(rowSums(cbind(reciprocityMatrix[,t]),na.rm = T)/t)+alpha2*{sum(contributorsMatrix[t])/(t*n_agents)}
      
       }
  }
  list(rankOrderCor = cor(groupRank,rank(-trueValues),method = 'kendall'),startingIncentiveAvg = mean(startingIncentive),finalIncentiveAvg = mean(incentiveMatrix[,t]),
       finalIncentiveVar = var(incentiveMatrix[,t]), PG_Generated = sum(payoutMatrix,na.rm = T))
}



