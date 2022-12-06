
#' testMethods
#'
#'Generates design matrix with user-specified grouping structure and runs user-specified group regularization methods
#'
#'
#' @param a Within group correlation (value from -1 to 1) for data generation
#' @param b Between group correlation (value from -1 to 1) for data generation
#' @param n Sample size for data generation
#' @param p Number of input variables for data generation
#' @param ng Number of groups for data generation
#' @param truegroup True grouping structure of input variables (p-dim vector with numerical true group labels)
#' @param withinGroupNoise Proportion of coefficients that switch from category 0 to category 1
#' @param chooseBetas True or false variable that indicates to take only last ng coefficients from each group to run methods on
#' @param methods Specified method ("lasso","ridge","group lasso","threshold group lasso","sparse group lasso","elastic net","group bridge")
#'
#' @return
#' @export
#'
#' @examples testMethods(a=0.9,b = 0.3,n=600,p=500,ng=10,methods=c("lasso","ridge","group lasso","threshold group lasso","elastic net","sparse group lasso","group bridge"))
#'
testMethods=function(a=0.9, b=0.5, n=500, p=500, ng=10, truegroup = sort(c(rep(1:(p/ng),ng))),withinGroupNoise=0,chooseBetas=FALSE,methods="lasso"){
  #LIBRARIES

  vp=as.vector(table(truegroup))

  SIGMA=matrix(a*b,p,p)
  diag(SIGMA)=1
  for (g in 1: ng) {
    pg=vp[g]
    sigma.tp=diag(pg)
    for (i in 1: vp[g]) {
      for (j in 1:(i-1)) {

        sigma.tp[i,j]=sigma.tp[j,i]=a
      }
    }
    SIGMA[which(truegroup==g),which(truegroup==g)]=sigma.tp
  }
  X = MASS::mvrnorm(n,mu=rep(0,p),Sigma=SIGMA);
  print(dim(X))

  true.betas = c()
  for (h in 1:ng){
    if (h %%2==0){
      temp = rep(1,p/ng)
    }else
      temp = rep(0,p/ng)
    true.betas=c(true.betas,temp)
  }
  switch.indices <- sample(seq(p),withinGroupNoise*p)#adding noise
  for (i in switch.indices){
    true.betas[i]=abs(true.betas[i]-1)
  }

  if (chooseBetas==TRUE){
    chosenBetas = c()
    Xsmall=matrix(0,nrow=n,ncol=(ng)^2)
    truegroup=sort(rep(1:ng,ng))
    count=0
    for (i in 1:p){
      for (j in 1:ng){
        if ((j)*(p/ng)-i >0 & (j)*(p/ng)-i < ng+1)
          count=count+1
        # chosenBetas=c(chosenBetas,true.betas[i])
        chosenBetas[count]=true.betas[i]
        Xsmall[,count]=X[,count]
      }
    }#for g in 1:ng
    true.betas=chosenBetas
    X=Xsmall
    print(length(true.betas))
    print(count)
  }#if chooseBetas is true
  Y = X %*% true.betas
  true_betas=true.betas
  #LASSO=======================================
  if ("lasso" %in% methods){
    print("Lasso=================================")
    cv_model <-glmnet::cv.glmnet(X,Y, alpha=1, nfolds=10)
    best_lambda <-cv_model$lambda.min
    best_lambda
    lasso_model <- glmnet::glmnet(X,Y,alpha=1,lambda=best_lambda)
    lasso.betas=coef(lasso_model)[2:length(coef(lasso_model))]
    cat("L2diff=",L2diff(true_betas,lasso.betas),'\n')
    cat("Sign match percentage =", SignMatch(true_betas,lasso.betas),'\n')
    lasso.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(lasso.betas)), reference =as.factor(true_betas))
    print(lasso.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = lasso.betas, weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Lasso")
  }#LASSO
  #RIDGE=======================================
  if ("ridge" %in% methods){
    print("Ridge=================================")
    cv_model2 <-glmnet::cv.glmnet(X,Y, alpha=0, nfolds=10)
    best_lambda <- cv_model2$lambda.min
    ridge_model <- glmnet::glmnet(X,Y,alpha=0,lambda=best_lambda)
    ridge.betas=coef(ridge_model)[2:length(coef(ridge_model))]
    cat("L2diff=",L2diff(true_betas,ridge.betas),'\n')
    cat("Sign match percentage =",SignMatch(true_betas,ridge.betas),'\n')
    rounded.ridge.betas <- betaRound(ridge.betas)
    ridge.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(ridge.betas)), reference =as.factor(true_betas))
    print(ridge.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = ridge.betas, weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Ridge")
  }#RIDGE
  #GRPLASSO=======================================
  if ("group lasso" %in% methods){
    print("Group Lasso=================================")
    cv_model <-gglasso::cv.gglasso(X,Y, group=truegroup, nfolds=10, loss="ls") #bruh this literally takes forever to run
    best_lambda <- cv_model$lambda.min
    grplasso_model <- gglasso::gglasso(X,Y,lambda=best_lambda,group=truegroup,loss="ls",intercept=F)
    grplasso.betas=coef(grplasso_model)[2:length(coef(grplasso_model))]
    cat("L2diff=",L2diff(true_betas,grplasso.betas),'\n')
    cat("Sign match percentage =",SignMatch(true_betas,grplasso.betas),'\n')
    rounded.grplasso.betas <- betaRound(grplasso.betas)
    grplasso.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(grplasso.betas)), reference =as.factor(true_betas))
    print(grplasso.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = grplasso.betas, weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Group Lasso")
  }#GRPLASSO

  #THRESHGRPLASSO=======================================
  if ("threshold group lasso" %in% methods){
    print("Threshold Group Lasso=================================")
    cv_model <-gglasso::cv.gglasso(X,Y, group=truegroup, nfolds=10, loss="ls")
    best_lambda <- cv_model$lambda.min
    grplasso_model <- gglasso::gglasso(X,Y,lambda=best_lambda,group=truegroup,loss="ls",intercept=F)
    grplasso.betas=coef(grplasso_model)[2:length(coef(grplasso_model))]
    thresh.betas <- thresh_within(grp.ind = truegroup, est =as.matrix(grplasso.betas,ncol=1),delta = 0.01)
    cat("L2diff=",L2diff(true_betas,thresh.betas),'\n')
    cat("Sign match percentage =",SignMatch(true_betas,thresh.betas),'\n')
    rounded.thresh.betas <- betaRound(thresh.betas)
    thresh.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(thresh.betas)), reference =as.factor(true_betas))
    print(thresh.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = thresh.betas, weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Thresholded Group Lasso")
  }#THRESHGRPLASSO
  #ELASTICNET=======================================
  if ("elastic net" %in% methods){
    print("Elastic Net=================================")
    df=data.frame(X,Y)
    names(df)=c(rep(1:dim(X)[2]),"seesaw")
    control <- caret::trainControl(method="repeatedcv",number=10,repeats=1, search="random",verboseIter=TRUE)
    elastic_model <- caret::train(seesaw~.,data=df,method="glmnet", preProcess=c("center","scale"),tuneLength=25,trControl=control)
    enet.betas <-stats::predict(elastic_model$finalModel,type="coefficients", s=elastic_model$bestTune$lambda)
    cat("L2diff=",L2diff(true_betas,enet.betas[-1]),'\n')
    cat("Sign match percentage =",SignMatch(true_betas,enet.betas[-1]),'\n')
    rounded.enet.betas <- betaRound(enet.betas[-1])
    enet.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(enet.betas[-1])), reference =as.factor(true_betas))
    print(enet.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = enet.betas[-1], weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Elastic Net")
  }#ELASTIC NET
  #SPARSEGROUPLASSO=======================================
  if ("sparse group lasso" %in% methods){
    print("Sparse Group Lasso=================================")
    sparse.fit <- sparsegl::cv.sparsegl(X,Y, group=truegroup,nfolds = 10)
    sparse.coef <- coef(sparse.fit)[-1]
    cat("L2diff=",L2diff(true_betas,sparse.coef),'\n')
    cat("Sign match percentage =",SignMatch(true_betas,sparse.coef),'\n')
    rounded.sparse.betas <- betaRound(sparse.coef)
    sparse.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(sparse.coef)), reference =as.factor(true_betas))
    print(sparse.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = sparse.coef, weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Sparse Group Lasso")
  }#SPARSE GROUP LASSO

  #GROUP BRIDGE=======================================
  if ("group bridge" %in% methods){
    print("Group Bridge=================================")
    gbridge <-grpreg::gBridge(X, Y, group=truegroup, family="gaussian", lambda=0.2, lambda.min={if (nrow(X) > ncol(X)) .001
      else .05}, lambda.max=20, alpha=1, eps=.001, delta=1e-7, max.iter=10000,
      gamma=0.5, warn=TRUE, returnX=FALSE, intercept=FALSE)
    group.bridge.betas= gbridge$beta[-1]
    cat("L2diff=",L2diff(true_betas,group.bridge.betas),'\n')
    cat("Sign match percentage =",SignMatch(true_betas,group.bridge.betas),'\n')
    rounded.gbridge.betas <- betaRound(group.bridge.betas)
    gbridge.confusion.matrix <- caret::confusionMatrix(data=as.factor(betaRound(group.bridge.betas)), reference =as.factor(true_betas))
    print(gbridge.confusion.matrix)
    PRROC_obj <- PRROC::roc.curve(scores.class0 = group.bridge.betas, weights.class0=true_betas, curve=TRUE)
    plot(PRROC_obj,main="Group Bridge")
  }#GROUPBRIDGE

  return()
}#testMethods


#' L2diff
#'
#' Calculates the l2 difference between two vectors of the same length
#'
#' @param a A vector
#' @param b A vector
#'
#' @return The l2 difference of the vectors
#' @export
#'
#' @examples L2diff(c(1,2,3,4),c(3,6,9,3))
#'
L2diff=function(a,b){
  norm(a-b,type="2")
}#L2 diff


#' Signmatch
#'
#' Takes pairs of elements between two vectors and sees whether or not both elements are in the same category (abs value less than 0.2, abs value greater than 0.2).
#' Calculates percentage of pairs that are in the same category
#'
#' @param a A vector
#' @param b A vector
#'
#' @return A matrix of True/False positive matches.
#' @export
#'
#' @examples SignMatch(c(1,2,3,4),c(3,6,9,3))
#'
SignMatch=function(a,b){
  matches=0
  for (i in 1:length(a)){
    if (((a[i] < 0.2 & a[i] > -0.2) & (b[i] < 0.2 & b[i] > -0.2))
        | ((a[i] > 0.2 | a[i] < -0.2) & (b[i] > 0.2 | b[i] < -0.2))) {
      matches = matches + 1
    }

  }
  return(matches/length(a))
} #signmatch

#' betaRound
#'
#' Rounds each element in the vector (to 0 or 1) based off of a threshold (abs value <0.2)
#'
#' @param a A vector
#'
#' @return rounded vector
#' @export
#'
#' @examples betaRound(c(1.30,0.02,-1.17))
#'
betaRound=function(a){
  b=c()
  for(i in 1:length(a)){
    if (a[i]<0.2 & a[i]>-0.2){
      b[i]=0
    }else
      b[i]=1
  }
  return(b)
}#beta round


#' thresh_within
#'
#' Runs thresholded group lasso on group lasso output. Created by: Sumanta Basu
#'
#' @param grp.ind a p-dim vector of group indices, e.g. 1, 1, 1, 2, 2, 3, 3, 3, 4, 4
#' @param est a p-dim vector (beta-hat) containing the group lasso solution
#' @param delta scalar, controlling the amount of threshold
#'
#' @return est.out
#' @export
#'
#' @examples thresh_within(grp.ind = c(0,0,1,1,2,2,3), est =as.matrix(c(1,2,3,4,7,6,5),ncol=1),delta = 0.01)
#'
thresh_within = function(grp.ind, est, delta){
  G = length(unique(grp.ind))
  p = nrow(est)
  l = ncol(est)
  delta.length = length(delta)

  if (delta.length == 1){
    for (g in 1:G){
      g.id = which(grp.ind == g)
      for (ll in 1:l){
        temp = sqrt(sum(est[g.id, ll]^2))
        if (temp)
          est[g.id, ll] = est[g.id, ll]*ifelse(abs(est[g.id, ll]) > delta*temp, 1, 0)
      }
    }
    est.out = est
  }

  else{
    est.D.abs = abs(est)
    for (g in 1:G){
      g.id = which(grp.ind == g)
      for (ll in 1:l){
        temp = sqrt(sum(est.D.abs[g.id, ll]^2))
        if (temp)
          est.D.abs[g.id, ll] = est.D.abs[g.id, ll]/temp
      }
    }

    est.out = array(0, c(p, l, delta.length))

    for (i in 1:delta.length){
      est.out[,,i] = est*(est.D.abs > delta[i])
    }

  }


  return(est.out)
}


