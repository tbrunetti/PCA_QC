# if using half-rack VM, enviromental proxy setting is required
#Sys.setenv(http_proxy="http://cloud-proxy:3128")
#Sys.setenv(https_proxy="http://cloud-proxy:3128")
#prcomp is used instread of princomp here to perform PCA since:
# 1) princomp is limited to experiments where observations >> variables
# 2) prcomp can be used to easily visualize clusters using PCs and metadata

# counts_file is the FULL path to the counts matrix in CSV format--MUST HAVE HEADER AND ROW NAMES!!!
#             counts must be rounded to nearest integer
#             row names should be gene/transcript names
#             headers should be BIDs or sample ID names
# metadata_file is the FULL path to the metadata matrix in CSV format--MUST HAVE HEADER AND ROW NAMES!!!
#             row names should be BIDs or sample ID names
#             headers should be the name of the metadata measurement (i.e. sex, diagnosis, chrM contamination, etc...)
# to keep it straight, it makes sense to have n x m counts files and a m x u because the matrices can be multiplied since
# their inner components are the same dimensions ( nxm matrix [dotproduct] mxu matrix = nxu matrix)


# ggfortify and ggplot are required for plotting and clustering PCs to metadatas
library(factoextra)
library(ggfortify)
library(ggplot2)
library(data.table)
library(gplots)


################################################################################
##function: prep_data(counts_file, metadata_file, output_name)                ##
##input: results returned from output_file in addtion to counts_file path     ##
##       and metadata_file path                                               ##
##output: generates PDF of PC plot results (scree plot and PC1 vs PC2) and    ##    
##        calls function model_pcs(pca_matrix, metadata_file, output_name,    ##
##        cat_vars, num_vars, all_vars) -- added tryCatch error handling      ##
################################################################################ 
prep_data <- function(counts_file, metadata_file, output_name){
  print(counts_file)
  print(metadata_file)
  print(output_name)
  tryCatch({
    pdf(file=paste(output_name, '.pdf', sep = ''))
    varTypes <- read.csv(file=metadata_file,nrows=2,header = FALSE)
    
    #determines if variable is categorical or numerical/continuous
    varTypes <- varTypes[,colSums(is.na(varTypes))<nrow(varTypes)]
    cat_vars <- c()
    num_vars <- c()
    for(i in seq(1, dim(varTypes)[2])){
      if(toupper(trimws(varTypes[1, i], which = "both")) == "C"){
        cat_vars <- c(cat_vars, as.character(varTypes[2,i]))
      }
      else if(toupper(trimws(varTypes[1, i], which = "both")) == "N"){
        num_vars <- c(num_vars,as.character(varTypes[2,i]))
      }
      else{
        stop("invalid type")
      }
    }
    
    
    # read in gene expression table and metadata row.names needs to be set to the column number where the row names are listed
    # it is important to set check.names=FALSE so that R doesn't try to change any numerical names into odd characters
    #total_raw_counts <- read.table(counts_file, header=TRUE, row.names=432, check.names=FALSE, sep=",") ## check row.names solution
    total_raw_counts <- read.table(counts_file, header=TRUE, row.names = 1, check.names=FALSE, sep=",") ## check row.names solution
    metadata <- read.table(metadata_file, header=TRUE, skip=1, row.names=1, check.names=FALSE, sep=",")
  
    total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts) #converts table into dataframe
    metadata_dataframe <- as.data.frame(metadata) #converts table into dataframe
    
    print(total_raw_read_counts_dataframe)
    print(metadata_dataframe)
    #pre-filter step, remove any ROWS that have zero counts in the count matrix, it does not mean that a sample cannot
    # have zero counts, we are just removing the genes where no counts exist across all samples
    total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe[rowSums(total_raw_read_counts_dataframe)!=0, ]
    
    total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe + 1 #add pseudo-count, why??  Because the log(0) = Does not exist, so to address this issue add 1 to all counts
    transformed_dataframe <- log(total_raw_read_counts_dataframe) # take the log of the counts data, this will help normalize the data
  
    # remove BIDs that have missing or strange data
    #transformed_dataframe <- within(transformed_dataframe, rm('nothing'))
    #metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('nothing'), ]
    
    #remove rows that have a sum of zero
    transformed_dataframe <- transformed_dataframe[rowSums(transformed_dataframe)!=0, ]
    
    # sort dataframes.  Dataframes MUST be sorted by colnames in the counts matrix
    # and sorted by the rownames in the metadata matrix.  This ensures that the sample names properly match
    # each other between counts matix and metadata matrix.  Note, you will see the word 'TRUE' printed
    # if they are properly match, else you will see 'FALSE' in which case you need to sort
    counts_sorted <- transformed_dataframe[,order(colnames(transformed_dataframe))]
    metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]
    if(all(rownames(metadata_sorted)==colnames(counts_sorted))){
      print("TRUE")
    }else{
      print("FALSE")
      stop("Error:  Sorted matrices are incorrect")
    }
    
    print(counts_sorted)
    print(metadata_sorted)
    
    # Calculation of the prinipal components using prcomp
    # the t() function means to take the transpose of the counts matrix
    # scale.=TRUE means PCs will be based on corrlation matrix and not covariance matrix
    # corr is better if scales of variables are very different
    # center=TRUE uses the centers the data, use centroid/mean
    print("working")
    pca_matrix <- prcomp(t(counts_sorted), center=TRUE, scale. = TRUE)
    print("done")
    # plot the PCs againts how much variance each PC contributes to the data set
    # looking to incorporate the min number of PCs before "elbowing-effect" into the model unless
    # a PC is strongly correlated with a variable in the metadata set, in which case,
    # just regress out the variable rather than the PC
    screeplot(pca_matrix, type ="l")
    print(fviz_eig(pca_matrix))
    
    getMetadata_headers = colnames(metadata_sorted)
    print(pca_matrix$x)
    
    for(var in getMetadata_headers){
      print(var)
      if(var %in% cat_vars){
        metadata_sorted[,var] <- (as.factor(metadata_sorted[,var]))
        print(autoplot(prcomp(t(counts_sorted)), data = metadata_sorted, colour = var))
      }else if(var %in% num_vars){
        print(autoplot(prcomp(t(counts_sorted)), data = metadata_sorted, colour = var) + scale_color_gradientn(colours = terrain.colors(7)))
      }
    }
    
    write.csv(pca_matrix$x, file=paste(output_name, "_pca_matrix_output.csv", sep = ''))
    write.csv(metadata_sorted, file=paste(output_name, "_metadata_sorted_output.csv", sep=''))
    return(model_pcs(pca_matrix=pca_matrix, metadata_file=metadata_sorted, output_name=output_name, cat_vars = cat_vars, num_vars = num_vars, all_vars=getMetadata_headers))}
    ,
    error=function(error_message){
      message("Oops, prep_data(counts_file, metadata_file, output_name) encountered an error")
      message(error_message)
      return(NA)
    }
  )
}



# ***************************Function to correlate PCs with metadata******************************
################################################################################
##function: model_pcs(pca_matrix, metadata_file, output_name, cat_vars,       ##
##          num_vars, all_vars)                                               ##
##input: output from prep_data(counts_file, metadata_file, output_name)       ##
##output: PDF of heatmap with linear models for every variable listed in the  ##
##        metadata up to the first 10 PCs annotated with -log10(pvals) and    ##
##        adjusted-Rsq values -- added tryCatch error handling                ##
################################################################################ 
model_pcs <- function(pca_matrix, metadata_file, output_name, cat_vars, num_vars, all_vars){
  tryCatch({
    factors_affecting_pcs=list()
    # create heatmeap to visualize PC and metadata correlations
    # create a dataframe to store all the -log10 p-values and adjusted R squared vals for the visualization of PCA data
    pvalues <- data.frame()
    adjRsq <- data.frame()
    # iterate through all PCs and get -log10(p-value) and adjusted R-squared value to every PC correlated to each
    # metadata column listed
    
    # pick whichever value is smallest (10 PCs or total variables in metadata)
    total_pcs = 10
    if(length(all_vars)<10){
      total_pcs = length(all_vars)
    }
    
    for (pc in seq(1,dim(pca_matrix$rotation)[2])){
      print(pc)
      for (variable in seq(1,total_pcs)){
        print(all_vars[variable])
        if(all_vars[variable] %in% cat_vars){
          linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_file[, all_vars[variable]])))
        }
        else if(all_vars[variable] %in% num_vars){
          linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_file[,all_vars[variable]]))
        }
        else{
          stop("variable not listed as categorical or numerical.  Please correct.")
        }
        factors_affecting_pcs[[all_vars[variable]]][[as.character(pc)]]=list()
        factors_affecting_pcs[[all_vars[variable]]][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
        factors_affecting_pcs[[all_vars[variable]]][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
        pvalues[all_vars[variable], pc] <- unlist(factors_affecting_pcs[all_vars[variable]][[1]][[pc]][[2]])
        adjRsq[all_vars[variable], pc] <- unlist(factors_affecting_pcs[all_vars[variable]][[1]][[pc]][[1]])
        }
    }
    # get the row and column names match the p-values
    #print(factors_affecting_pcs)
    rownames(pvalues) <- all_vars
    print(pvalues)
    colnames(pvalues) <- unlist(lapply(seq(1,total_pcs),function(x) paste(c('PC',x),collapse='')))
    print(pvalues)
    
    # get the row and column names matching the adj R-sq values
    rownames(adjRsq) <- all_vars
    colnames(adjRsq) <- unlist(lapply(seq(1,total_pcs),function(x) paste(c('PC',x),collapse='')))
    
    # round all -log10(pvalue) in the dataframe to three decimal places
    is.num <- sapply(pvalues, is.numeric)
    pvalues[is.num] <- lapply(pvalues[is.num], round, 3)

    # create a heatmap of these values, value is -log10(p-val) and color is the adj R-sq value
    
    heatmap.2(as.matrix(adjRsq), cellnote=pvalues, notecol = "black", notecex = 0.5, cexRow = 0.3, dendrogram = "none", col=colorRampPalette(c("white", "yellow", "red"))(10))
    print("heatmap completed")
    dev.off()
    }
    ,
    error=function(error_message){
      message("An errror occurred while trying to calculate and build heatmap")
      message(error_message)
      return(NA)
    }
  )
}

#*********************************************End of function*********************************************

#---------------------------regress out variables of interest-------------------------------------------------
# for each variable that is to be regressed a linear model must be made and the residuals of the linear model
# must be extract and replace the old PC matrix
# regress out FlowcellBatch
#for (pc in seq(1,dim(pca_matrix$rotation)[2])){
#  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'FlowcellBatch'])))
#  pca_matrix$x[,pc]  <- linear_model$residuals
#}
# regress out Sex
#for (pc in seq(1,dim(pca_matrix$rotation)[2])){
#  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Sex'])))
#  pca_matrix$x[,pc]  <- linear_model$residuals
#}
# regress out 5 prime to 3 prime bias
#for (pc in seq(1,dim(pca_matrix$rotation)[2])){
#  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_MEDIAN_5PRIME_TO_3PRIME_BIAS']))
#  pca_matrix$x[,pc]  <- linear_model$residuals
#}
# regress out RIN
#for (pc in seq(1,dim(pca_matrix$rotation)[2])){
#  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'RIN']))
#  pca_matrix$x[,pc]  <- linear_model$residuals
#}
# regress out PMI
#for (pc in seq(1,dim(pca_matrix$rotation)[2])){
#  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'PMI']))
#  pca_matrix$x[,pc]  <- linear_model$residuals
#}

# call function again on regressed out variables
#model_pcs(pca_matrix)
# writes the variables that were regressed out, NOTE this must be changed manually!!!!!!
#mtext("vars regressed: FlowcellBatch, Sex, UF 5-3 bias, RIN, PMI", side=3, line=0)
# saves and closes newly created PDF
#dev.off()



################################################################################
##function: get_user_input()                                                  ##
##input:  None -- prompts user for location of counts_file                    ##
##output: Takes user input and proceeds to function                           ##
##        val_user_input(verification_1)                                      ##
################################################################################ 
get_user_input <- function(){
  cat("Full path to counts file (no whitespaces): ")
  count_file <- trimws(readLines("stdin", n=1), which="both")
  cat("You entered the following: ", count_file)
  cat('\n')
  cat('Is this correct? (Y/N/Q) ')
  verification_1 <- trimws(readLines('stdin', n=1), which="both")
  return(val_user_input(counts_file = count_file, verification_1=verification_1))
}
  
################################################################################
##function: val_user_input(verification_1)                                    ##
##input:  verification value returned from function get_user_input()          ##
##output: If count_file path is verified as correct, proceeds to function     ##
##        output <- get_metadata(); if user specifies no or n then recursively##
##        calls get_user_input() until user corretly validate file path or    ##
##        user quits program; if user validates metadata file then calls      ##
##        function output_file(counts_file, metadata_file); else if user does ##
##        not validate the metadata_file output<-get_metadata() is recursively##
##        called until file is validated or user exits program                ##
################################################################################ 
val_user_input <- function(counts_file, verification_1){
  if ((toupper(verification_1) == "YES") | (toupper(verification_1) == "Y")){
    output <- get_metadata()
    if (toupper(unlist(output[2]))=="YES" | toupper(unlist(output[2]))=='Y'){
      return(output_file(counts_file = counts_file, metadata_file=unlist(output[1])))
    }else if ((toupper(unlist(output[2])) == "NO") | (toupper(unlist(output[2])) == 'N')){
      output <- get_metadata()
      return(val_user_input(counts_file = counts_file, verification_1 = verification_1))
    }else{
      stop("Exiting Program.  Goodbye!")
    }
  }else if((toupper(verification_1) == "NO") | (toupper(verification_1) == 'N')){
    return(get_user_input())
  }else{
    stop("Exiting Program.  Goodbye!")
  }
}


################################################################################
##function: get_metadata()                                                    ##
##input:  None -- prompts user for location of metadata_file (properly        ## 
##        formatted)                                                          ##
##output: Takes user input and returns list with meta data file path and      ## 
##        verifcation                                                         ##
################################################################################ 
get_metadata <- function(){
  cat("Full path to meta data file (no whitespaces): ")
  metadata_file <- trimws(readLines('stdin', n=1), which="both")
  cat("You entered the following: ", metadata_file)
  cat("\n")
  cat("Is this correct? (Y/N/Q) ")
  verification_2 <- trimws(readLines('stdin', n=1), which="both")
  print(metadata_file)
  print(verification_2)
  return(list(metadata_file, verification_2))
}



################################################################################
##function: output_file()                                                     ##
##input: counts_file and metadata_file names and paths from function:         ##
##       val_user_input(verification_1)                                       ##
##output: If output_file is verified then calls function:                     ## 
##        prep_data(counts_file, metadata_file, output_name); else if         ##
##        verification is equal to no, then recursively calls output_file()   ##
##        until verification is yes or user quits (any verification           ##
##        character not equal to yes or no)                                   ##
################################################################################ 
output_file <- function(counts_file, metadata_file){
  cat("Name of output file (no whitespaces): ")
  output_file <- trimws(readLines("stdin", n=1), which="both")
  cat("You entered the following: ", output_file)
  cat('\n')
  cat('Is this correct? (Y/N/Q) ')
  verification_output <- trimws(readLines('stdin', n=1), which="both")
  print(verification_output)
  if ((toupper(verification_output) == "YES") | (toupper(verification_output) == 'Y')){
    return(prep_data(counts_file = counts_file, metadata_file=metadata_file, output_name=output_file))
  }else if ((toupper(verification_output) == "NO") | (toupper(verification_output) == 'N')){
    return(output_file())
  }else{
    stop("Exiting Program.  Goodbye!")
  }
}

get_user_input()





##-----------------------------------------------------END OF WORKING SCRIPT--------------------------------------------------------##
##-----------------------------------------EVERYTHING BELOW THIS LINE IS MANUAL VALIDATION------------------------------------------##

##validation##
validation_test <- function(){
  pca_matrix <- read.csv("github_repositories/PCA_QC/pca_matrix_output.csv", header = T, row.names = 1)
  metadata <- read.table("github_repositories/PCA_QC/metadata_sorted_output.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
  #categorical variables
  linear_model <- lm(pca_matrix$PC1 ~ na.omit(as.factor(metadata$sex)))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC2 ~ na.omit(as.factor(metadata$sex)))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC3 ~ na.omit(as.factor(metadata$sex)))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC1 ~ na.omit(as.factor(metadata$diagnosis)))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC2 ~ na.omit(as.factor(metadata$diagnosis)))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC3 ~ na.omit(as.factor(metadata$diagnosis)))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  #numeric variables
  linear_model <- lm(pca_matrix$PC1 ~ na.omit(metadata$contamination))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC2 ~ na.omit(metadata$contamination))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
  
  linear_model <- lm(pca_matrix$PC3 ~ na.omit(metadata$contamination))
  print(summary(linear_model)$adj.r.squared)
  print(-log10(anova(linear_model)$Pr[1]))
}

