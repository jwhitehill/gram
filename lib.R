# ----------------------------------------------------------------------------
#   
#   GRAM - Gaussian Rating Aggregation Model
#
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# Cecilia Aguerrebere <caguerrebere@gmail.com>
# ----------------------------------------------------------------------------
  
MStep <- function(d, workers, e_step_out, m_step_out, useMu, sigmaType ){
  
  for ( worker in workers ) {
    idxs = which(d$worker == worker)
    mu = 0
    for ( idx in idxs ) {
      video = d$video[idx]
      label = d$helpOthers[idx]
      mu = mu + label - e_step_out$ms[e_step_out$video==video]
    }
    mu = mu/length(idxs)
    sigma2 = 0
    for ( idx in idxs ) {
      video = d$video[idx]
      label = d$helpOthers[idx]
      sigma2 = sigma2 + (label - mu)^2 - 2*(label - mu)*e_step_out$ms[e_step_out$video==video] +
        e_step_out$ms[e_step_out$video==video]^2 + e_step_out$s2s[e_step_out$video==video]
    }
    
    if ( useMu == 1 ){
      m_step_out$mus[m_step_out$worker==worker] = mu      
    } else {
      m_step_out$mus[m_step_out$worker==worker] = 0 
    }
    
    if ( sigmaType == "pretest" ){
      m_step_out$sigma2s[m_step_out$worker==worker] = 1/(mean(d$pretest[idxs]) + EPS)^(0.5)
    } else if ( sigmaType == "sigma" ){
      m_step_out$sigma2s[m_step_out$worker==worker] = sigma2 
    } else {
      m_step_out$sigma2s[m_step_out$worker==worker] = 1  
    }
  }
  return(m_step_out)
}

EStep <- function(d, videos, m_step_out, V, U, EPS){
  e_step_out=data.frame("ms" = double(length(videos)),
                        "s2s" = double(length(videos)),
                        "video" = videos)
  for ( video in videos ) {
    s2inv = 1/V^2
    m = U/V^2
    idxs = which(d$video == video)
    for ( idx in idxs ) {
      worker = d$worker[idx]
      label = d$helpOthers[idx]
      sigma2 = m_step_out$sigma2s[m_step_out$worker==worker] + EPS
      s2inv = s2inv + 1/sigma2
      m = m + (label - m_step_out$mus[m_step_out$worker==worker])/sigma2
    }
    s2 = s2inv^(-1)
    e_step_out$s2s[e_step_out$video == video] = s2
    m = s2 * m
    e_step_out$ms[e_step_out$video == video] = m
  } 
  return(e_step_out)
}

EM <- function( d, Niter, V, U, EPS, useMu, sigmaType, verbose, verbose_pval ){
  workers = unique(d$worker)
  videos = unique(d$video)
  m_step_out=data.frame("mus"=rep(0.0,length(workers)),
                        "sigma2s"=rep(1.0,length(workers)),
                        "worker"=workers)
  
  for (i in seq(1,Niter)) {
    e_step_out = EStep(d, videos, m_step_out, V, U, EPS)
    m_step_out = MStep(d, workers, e_step_out, m_step_out, useMu, sigmaType )

    cc=printResults(d, videos, e_step_out, i, verbose )  
  }
  
  if ( verbose_pval == 1 ){
    print(paste("Spearman =",
                format(cc$cc_spearman,digits=3), 
                "p-val =", format(cc$pval_spearman, digits=3),
                ", Pearson =",
                format(cc$cc_pearson,digits=3),
                "p-val =", format(cc$pval_pearson,digits=3)))    
  }
  
  return(e_step_out)
}

printResults <- function(d, videos, e_step_out, niter, verbose ){
  ys = c()
  yhats = c()
  for ( video in videos ){
    idxs = which(d$video == video)
    avgLearningGains = mean(d$posttest[idxs] - d$pretest[idxs])
    ys=c(ys,avgLearningGains)
    yhats=c(yhats,e_step_out$ms[e_step_out$video==video])
  }
  if ( verbose == 1 ){
    cc_spearman=cor(yhats, ys,method="spearman")
    cc_pearson=cor(yhats, ys,method="pearson")
    
    print(paste("Iteration:",niter,", Spearman =",
                format(cc_spearman,digits=3),
                ", Pearson =",
                format(cc_pearson,digits=3)))    
  }
  
  cc_spearman=cor.test(yhats, ys,method="spearman",
                         alternative = "two.sided")
  cc_pearson=cor.test(yhats, ys,method="pearson",
                        alternative = "two.sided")
  
  return( list("cc_spearman"=cc_spearman$estimate,
               "pval_spearman"=cc_spearman$p.value,
               "cc_pearson"=cc_pearson$estimate,
               "pval_pearson"=cc_pearson$p.value) )
}

print_crossValidation_results <- function( results_allTrials, cor_type, Kvals, nTrials ){
  
  aux=results_allTrials[results_allTrials$cor_type==cor_type,]
  setDT(aux)[,"average":=mean(value),by=model_type]
  setDT(aux)[,"standardError":=mean(value_sd/sqrt(Kvals)),by=model_type]
  aux=unique(aux,by="model_type")
  
  print(paste("Number of trials:",nTrials,", K-folds:",Kvals,
              ", Correlation type:",cor_type))
  print(aux[,names(aux) %in% c("model_type","average","standardError"),with=FALSE])
  
}

assign_k_folds <- function( N, K ){
  
  indx=seq(1,N)
  indx_rand=sample(indx)
  kfolds_aux=cut(seq(1,N),breaks=K,labels=FALSE)
  kfolds=double(N)
  for (i in seq(1,K)){
    kfolds[indx_rand[kfolds_aux==i]]=i  
  }
  
  return( kfolds )
}

check_video_assignment <- function( kfolds, d, minFolds ){
  
  workers=unique(d$worker)
  if ( length(workers) != length(kfolds) ){
    print("ERROR: kfolds and workers dont match!")
  }
  d_aux=data.table(d)
  d_aux$kfold<-NA
  for (w in workers){
    d_aux$kfold[d_aux$worker==w]=kfolds[workers==w]
  }
  setDT(d_aux)[,"Nfolds_video":=uniqueN(kfold),by=video]
  if ( min(d_aux$Nfolds_video) >= minFolds ){
    assignmentNotOk=0
  } else {
    assignmentNotOk=1
  }
  return(assignmentNotOk) 
}
