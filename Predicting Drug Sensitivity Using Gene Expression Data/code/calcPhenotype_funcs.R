
# below functions used linearridge described in the paper
calcPhenotype_ridge <- function(testExprData, trainingExprData, trainingPtype, batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE)
{
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting Ridge Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  rrModel <- linearRidge(Resp ~ ., data=trainFrame)
  if(printOutput) cat("Done\n\nCalculating predicted phenotype...");
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=testFrame)
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}


# below function used ridge from glmnet 
calcPhenotype_ridge_glmnet <- function(testExprData, trainingExprData, trainingPtype, 
                                       batchCorrect="eb", 
                                       powerTransformPhenotype=TRUE, 
                                       removeLowVaryingGenes=.2, 
                                       minNumSamples=10, 
                                       selection=-1, 
                                       printOutput=TRUE)
{
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting glmnet ridge Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  
  ridge_cv <- cv.glmnet(as.matrix(t(homData$train[keepRows, ])), trainingPtype, alpha = 0)
  best_lambda <- ridge_cv$lambda.min
  rrModel <- glmnet(as.matrix(t(homData$train[keepRows, ])), trainingPtype, alpha = 0, lambda = best_lambda)
  if(printOutput){
    cat("Done\n\nCalculating predicted phenotype...", "\n")
    cat("best lambda min is:", best_lambda)
  }
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newx=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newx=as.matrix(testFrame) )
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}

# below function used lasso from glmnet
# user can specify th print out the non-zero coefficients
calcPhenotype_lasso <- function(testExprData, trainingExprData, trainingPtype, 
                                       batchCorrect="eb", 
                                       powerTransformPhenotype=TRUE, 
                                       removeLowVaryingGenes=.2, 
                                       minNumSamples=10, 
                                       selection=-1, 
                                       printOutput=TRUE)
{
  set.seed(99)
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the lasso Regression model on our training data
  if(printOutput) cat("\nFitting glmnet lasso Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  
  ridge_cv <- cv.glmnet(as.matrix(t(homData$train[keepRows, ])), trainingPtype, alpha = 1)
  best_lambda <- ridge_cv$lambda.min
  rrModel <- glmnet(as.matrix(t(homData$train[keepRows, ])), trainingPtype, alpha = 1, lambda = best_lambda)
  if(printOutput){
    cat("Done\n\nCalculating predicted phenotype...", "\n")
    cat("best lambda min is:", best_lambda)
  }
  # print out non-zero coefficients
  myCoefs <- coef(rrModel, s="lambda.min");
  
  ## Asseble into a data.frame
  myResults <- data.frame(
    features = myCoefs@Dimnames[[1]][ which(myCoefs != 0 ) ], #intercept included
    coefs    = myCoefs              [ which(myCoefs != 0 ) ]  #intercept included
  )
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newx=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newx=as.matrix(testFrame) )
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(list(preds=preds, coeffs=myResults))
}


# below fuction used elestic net from glmnet
calcPhenotype_enet <- function(testExprData, trainingExprData, trainingPtype, 
                                batchCorrect="eb", 
                                powerTransformPhenotype=TRUE, 
                                removeLowVaryingGenes=.2, 
                                minNumSamples=10, 
                                selection=-1, 
                                printOutput=TRUE)
{
  set.seed(99)
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting glmnet enet Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  
  ridge_cv <- cv.glmnet(as.matrix(t(homData$train[keepRows, ])), trainingPtype, alpha = 0.5)
  best_lambda <- ridge_cv$lambda.min
  rrModel <- glmnet(as.matrix(t(homData$train[keepRows, ])), trainingPtype, alpha = 0.5, lambda = best_lambda)
  if(printOutput){
  cat("Done\n\nCalculating predicted phenotype...", "\n")
  cat("best lambda min is:", best_lambda)
  }
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newx=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newx=as.matrix(testFrame) )
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}




# below function used svm 
calcPhenotype_svm <- function(testExprData, trainingExprData, trainingPtype, 
                               batchCorrect="eb", 
                               powerTransformPhenotype=TRUE, 
                               removeLowVaryingGenes=.2, 
                               minNumSamples=10, 
                               selection=-1, 
                               printOutput=TRUE)
{
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting svm Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  
  # svm tunning parameter: gamma, cost
  svm.out <- tune.svm(Resp~., data = trainFrame)
  rrModel <- svm(Resp ~ ., data=trainFrame, 
                 type = 'eps-regression', 
                 kernel = 'linear',
                 gamma = svm.out$best.model$gamma, 
                 cost = svm.out$best.model$cost)
  
  if(printOutput){
    cat("Done\n\nCalculating predicted phenotype...")
    cat(" best gamma",svm.out$best.model$gamma, "best cost", svm.out$best.model$cost)
  }
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=as.matrix(testFrame) )
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}







# below function used ANN 
calcPhenotype_ann <- function(testExprData, trainingExprData, trainingPtype, 
                              batchCorrect="eb", 
                              powerTransformPhenotype=TRUE, 
                              removeLowVaryingGenes=.2, 
                              minNumSamples=10, 
                              selection=-1, 
                              printOutput=TRUE)
{
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting ANN Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  
  ann.out <- tune.nnet(Resp~., size = 1:10,
                       decay = 0.01,
                       data = trainFrame)
  
  
  rrModel <- nnet(Resp~.,size=ann.out$best.parameters$size, 
                  decay=0.01,
                  MaxNWts = 4e4,
                  linout=TRUE,
                  maxit = 5000,
                  data = trainFrame)
  
  if(printOutput){
    cat("Done\n\nCalculating predicted phenotype...", "/n")
    cat("best size is:", ann.out$best.parameters$size)
  }
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel,testFrame)
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}



# below function used PCR
calcPhenotype_pcr <- function(testExprData, trainingExprData, trainingPtype, 
                              batchCorrect="eb", 
                              powerTransformPhenotype=TRUE, 
                              removeLowVaryingGenes=.2, 
                              minNumSamples=10, 
                              selection=-1, 
                              printOutput=TRUE,
                              Ncomp)
{
  # Calculates a phenotype from gene expression microarray data, given a training set of expression data and known phenotype.
  #
  # Args:
  #   testExprData: The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
  #   trainingExprData: The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
  #   trainingPtype: The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
  #   batchCorrect: How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
  #   powerTransformPhenotype: Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
  #   removeLowVaryingGenes: What proportion of low varying genes should be removed? 20% be default
  #   minNumSamples: How many training and test samples are requried. Print an error if below this threshold
  #   selection: How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
  #
  # Returns:
  #   A vector of the estimated phenotype, in the same order as the columns of "testExprData".
  
  # check if the supplied data are of the correct classes
  if(class(testExprData) != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData) != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype) != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same lenght as the number of columns of the training expressin matrix.");
  
  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)
  
  # Do variable selection if specified. By default we remove 20% of least varying genes.
  # Otherwise, keep all genes.
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1)
  {
    keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
  }
  else 
    keepRows <- seq(1:nrow(homData$train))
  
  
  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    
    transForm <- powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  
  # create the PCR model on our training data
  if(printOutput) cat("\nFitting PCR Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  
  
  rrModel <- pcr(Resp~., 
                 data = trainFrame, 
                 scale = TRUE, 
                 validation = "CV")
  
  if(printOutput){
    cat("Done\n\nCalculating predicted phenotype...", "/n")
    validationplot(rrModel)
  }
  
  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  
  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test) == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel,testFrame, ncomp = Ncomp)
  }
  
  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");
  
  return(preds)
}








