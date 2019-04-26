library(data.table)
rm(list=ls())
source("lib.R")

# --------------------------------------------------------
# Section 6.1: Are subjective ratings correlated with 
#              average learning gains?
# --------------------------------------------------------
cat("Section 6.1: Experiments","\n")

# Load dataset
d <- read.csv("~/Research/Ceibal/Fundacion/code/r/oer/CrowdsourcingExplanations6.csv")

# Compute learning gains
d$lg=d$posttest-d$pretest

# GRAM model parameters
## Prior distribution parameters
U = 3
V = 1000
EPS = 1e-1
Niter = 50
verbose = 0
verbose_pval = 1

# GRAM model estimation
for ( useMu in c(0,1) ){
  for ( sigmaType in c("pretest","sigma","") ) {
    print(paste("GRAM: useMu =",useMu,", sigmaType =",sigmaType))
    gram_scores=EM( d, Niter, V, U, EPS, useMu, sigmaType, verbose, verbose_pval )  
  }
}

# Compute correlation for unweighted average
dt=data.table(d)
setDT(dt)[,"unwAvg":=mean(helpOthers),by=video]
setDT(dt)[,"avgLearnGains":=mean(lg),by=video]
dt=unique(dt,by="video")
print("Unweighted average:")
cc_spearman=cor.test(dt$unwAvg,dt$avgLearnGains,method="spearman",
                     alternative = "two.sided")
cc_pearson=cor.test(dt$unwAvg,dt$avgLearnGains,method="pearson",
                     alternative = "two.sided")
print(paste("Spearman =",format(cc_spearman$estimate,digits=2),
            ", p-val =",format(cc_spearman$p.value,digits=2),
             ", Pearson =",format(cc_pearson$estimate,digits=2),
            ", p-val =",format(cc_pearson$p.value,digits=2)))

# --------------------------------------------------------
# Section 6.2: Do subjective ratings predict the average 
#              learning gains for new students?
# --------------------------------------------------------
cat("Section 6.2: Experiments","\n")

# Load dataset
d <- read.csv("~/Research/Ceibal/Fundacion/code/r/oer/CrowdsourcingExplanations6.csv")

# Compute learning gains
d$lg=d$posttest-d$pretest

# GRAM model parameters
useMu = 1 
sigmaType = 0 # ("pretest","sigma") 
## Prior distribution parameters
U = 3
V = 1000
EPS = 1e-1
Niter = 50
verbose = 0
verbose_pval = 0

# Cross-validation setup
workers=unique(d$worker)
Nworkers=length(workers)
Kvals=3
nTrials=50
Nmodels=3 # GRAM, UnwAvg, lm_testScores
minFolds=2
results_allTrials=data.frame("nTrials"=rep(seq(1,nTrials),2*Nmodels)) 
results_allTrials$cor_type<-c(rep("spearman",nTrials*Nmodels),
                              rep("pearson",nTrials*Nmodels))
results_allTrials$model_type<-c(rep("GRAM",nTrials),
                                rep("UnwAvg",nTrials),
                                rep("lm_testScores",nTrials))
results_allTrials$value<-NA
results_allTrials$value_sd<-NA


for ( cc in seq(1,nTrials) ){   # For each trial
  
  print(paste("Running trial",cc,"of",nTrials))
  
  # Keep track of the results of this trial
  results=data.frame("kfold"=double(Kvals))
  results$cor_spearGRAM<-NA
  results$cor_pearsGRAM<-NA
  results$cor_pearsUnwAvg<-NA
  results$cor_spearUnwAvg<-NA
  results$cor_pearsTestScores<-NA
  results$cor_spearTestScores<-NA
  
  # Check that all videos are in at least two folds to 
  # ensure that the videos in the test set are always 
  # in the training set. Otherwise re-run the k-fold 
  # assignment.
  assignmentNotOk=1
  tt=1
  while( assignmentNotOk==1 ){
    kfolds=assign_k_folds( Nworkers, Kvals )  
    assignmentNotOk=check_video_assignment(kfolds,d,minFolds)
    if (tt > 1){
      print(paste("Checking k-fold assignment",tt,"times."))  
    }
    tt=tt+1
  }
  
  for ( k in seq(1,Kvals)){  # For each k-fold
    
    # Workers assigned to the test set
    test_workers=workers[kfolds==k]
    
    # Define training set (all data except that of test workers)
    d_train=d[!(d$worker %in% test_workers),]
    
    print(paste("Running k-fold:",k,"of",Kvals,
                "with",dim(d_train)[1],"workers for training."))
    
    # Compute GRAM scores from training set
    gram_scores=EM( d_train, Niter, V, U, EPS, useMu, sigmaType, 
                    verbose, verbose_pval )  
    d_train$ms=NA
    for ( video in unique(d_train$video)){
      d_train$ms[d_train$video==video]=gram_scores$ms[gram_scores$video==video]
    }
    
    # Compute unweighted average scores from training set
    d_train=data.table(d_train)
    setDT(d_train)[,"ms_avg":=mean(helpOthers),by=video]
    d_train=data.frame(d_train)
    
    # Define testing set
    d_test=data.table(d[d$worker %in% test_workers,])
    setDT(d_test)[,"avgLG":=mean(lg),by=video]
    d_test=unique(d_test,by="video")
    
    # Predict using learning gains to take as baseline
    lm_testScores = lm( lg ~ video,
                        data=d_train )
    d_test$testScores = predict(lm_testScores,d_test,type="response")
    
    # Keep unique video scores per video in the training set
    d_train=data.table(d_train)
    d_train=unique(d_train,by="video")
    
    # Get computed GRAM and unweighted average scores for the test videos
    d_test$video_scoreGRAM<-NA
    d_test$video_scoreUnwAvg<-NA
    for ( vv in unique(d_test$video) ){
      d_test$video_scoreGRAM[d_test$video==vv]=d_train$ms[d_train$video==vv]
      d_test$video_scoreUnwAvg[d_test$video==vv]=d_train$ms_avg[d_train$video==vv]
    }
    
    # Compute correlations of each score type to the average learning gains per video
    results$kfold[k]=k
    results$cor_pearsGRAM[k]=cor(d_test$video_scoreGRAM,
                                 d_test$avgLG,method="pearson")
    results$cor_spearGRAM[k]=cor(d_test$video_scoreGRAM,
                                 d_test$avgLG,method="spearman")
    results$cor_pearsUnwAvg[k]=cor(d_test$video_scoreUnwAvg,
                                   d_test$avgLG,method="pearson")
    results$cor_spearUnwAvg[k]=cor(d_test$video_scoreUnwAvg,
                                   d_test$avgLG,method="spearman")
    results$cor_pearsTestScores[k]=cor(d_test$testScores,
                                       d_test$avgLG,method="pearson")
    results$cor_spearTestScores[k]=cor(d_test$testScores,
                                       d_test$avgLG,method="spearman")
  }
  
  # Mean across folds for trial cc
  # -- Spearman
  results_allTrials$value[results_allTrials$nTrials==cc & results_allTrials$cor_type=="spearman" &
                      results_allTrials$model_type == "GRAM"]=mean(results$cor_spearGRAM)
  results_allTrials$value[results_allTrials$nTrials==cc & results_allTrials$cor_type=="spearman" &
                      results_allTrials$model_type == "UnwAvg"]=mean(results$cor_spearUnwAvg)
  results_allTrials$value[results_allTrials$nTrials==cc & results_allTrials$cor_type=="spearman" & 
                      results_allTrials$model_type == "lm_testScores"]=mean(results$cor_spearTestScores)
  # -- Pearson   
  results_allTrials$value[results_allTrials$nTrials==cc & results_allTrials$cor_type=="pearson" & 
                      results_allTrials$model_type == "GRAM"]=mean(results$cor_pearsGRAM)
  results_allTrials$value[results_allTrials$nTrials==cc & results_allTrials$cor_type=="pearson" &
                      results_allTrials$model_type == "UnwAvg"]=mean(results$cor_pearsUnwAvg)
  results_allTrials$value[results_allTrials$nTrials==cc & results_allTrials$cor_type=="pearson" & 
                      results_allTrials$model_type == "lm_testScores"]=mean(results$cor_pearsTestScores)
  
  # Standard deviation across folds for trial cc
  # -- Spearman
  results_allTrials$value_sd[results_allTrials$nTrials==cc & results_allTrials$cor_type=="spearman" &
                         results_allTrials$model_type == "GRAM"]=sd(results$cor_spearGRAM)
  results_allTrials$value_sd[results_allTrials$nTrials==cc & results_allTrials$cor_type=="spearman" &
                         results_allTrials$model_type == "UnwAvg"]=sd(results$cor_spearUnwAvg)
  results_allTrials$value_sd[results_allTrials$nTrials==cc & results_allTrials$cor_type=="spearman" & 
                         results_allTrials$model_type == "lm_testScores"]=sd(results$cor_spearTestScores)
  # -- Pearson    
  results_allTrials$value_sd[results_allTrials$nTrials==cc & results_allTrials$cor_type=="pearson" & 
                         results_allTrials$model_type == "GRAM"]=sd(results$cor_pearsGRAM)
  results_allTrials$value_sd[results_allTrials$nTrials==cc & results_allTrials$cor_type=="pearson" &
                         results_allTrials$model_type == "UnwAvg"]=sd(results$cor_pearsUnwAvg)
  results_allTrials$value_sd[results_allTrials$nTrials==cc & results_allTrials$cor_type=="pearson" & 
                         results_allTrials$model_type == "lm_testScores"]=sd(results$cor_pearsTestScores)
  
}

cor_type="spearman"
print_crossValidation_results( results_allTrials, cor_type, Kvals, nTrials )

cor_type="pearson"
print_crossValidation_results( results_allTrials, cor_type, Kvals, nTrials )


