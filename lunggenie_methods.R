# - check this exists - # 
mast="#########################################\n###   LungGENIE (LR-PCA imputation)  ###\n#########################################"

version = "\n\nversion 0.0.1"

contact = "\n\nJonathan L. Hess, PhD
Dept. of Psychiatry & Behavioral Sciences,
Norton College of Medicine,
SUNY Upstate Medical University
hessjo@upstate.edu"

mast_head = c(mast, version, contact)
message(mast_head)

# == load library
require(rtracklayer)
require(plyr)
require(readxl)
require(data.table)
require(rtracklayer)


# 1. Load pre-processed paired lung/blood data from GTEx donors
load_expr_data = function(path_to_data = NULL){
  # warning message
  if(is.null(path_to_data) == T){stop("Please provide a directory path")}
  # load data
  data_files = list.files(path_to_data, pattern=".Rdata", full.names=T)
  data_files = data_files[!grepl("Predictors|dtangle|Covariates|Scaling", data_files)]
  blood_expr = data.frame(t(readRDS(data_files[grepl("WholeBlood", data_files)])))
  lung_expr = data.frame(t(readRDS(data_files[!grepl("WholeBlood", data_files)])))
  # retain same samples
  common.sids = intersect(colnames(lung_expr), colnames(blood_expr))
  blood_expr <<- blood_expr[,colnames(blood_expr) %in% common.sids]
  lung_expr <<- lung_expr[,colnames(lung_expr) %in% common.sids]
}

## Quick convert ensembl gene ids to hgnc symbols
quick_convert_ensg_to_hgnc = function(gtf_path=NULL){
  if(is.null(gtf_path) == T){stop("Must supply path to GTF file compatible with data")}
  # convert blood gene expression ids to hgnc symbols
  convert = convert_ensg_to_hgnc(ensg = rownames(blood_expr), gtf_path = gtf_path)
  convert = convert[match(rownames(blood_expr), convert$gene_id), ]
  dat = data.frame(symbol = convert$gene_name, blood_expr)
  dat = dat[!is.na(dat$symbol),]
  dat = data.table(dat)
  dat = dat[,lapply(.SD, median),by=c("symbol")]
  dat = data.frame(dat)
  rownames(dat) = dat$symbol
  dat = dat[,!colnames(dat) %in% "symbol"]
  blood_expr <<- dat
  
  # convert lung gene expression ids to hgnc symbols
  convert = convert_ensg_to_hgnc(ensg = rownames(lung_expr), gtf_path = gtf_path)
  convert = convert[match(rownames(lung_expr), convert$gene_id), ]
  dat = data.frame(symbol = convert$gene_name, lung_expr)
  dat = dat[!is.na(dat$symbol),]
  dat = data.table(dat)
  dat = dat[,lapply(.SD, median),by=c("symbol")]
  dat = data.frame(dat)
  rownames(dat) = dat$symbol
  dat = dat[,!colnames(dat) %in% "symbol"]
  lung_expr <<- dat
}


 

# 2. load cross-validation performance 
load_cv_performance = function(){
  if(is.null(output_file_string)){stop("Please run cross-validation step to obtain gene-level prediction accuracies for LungGENIE")}
  return(data.frame(readRDS(output_file_string)))
}

# 3. principal component analysis in reference data
fit_pca = function(gene_list = NULL, autorun = FALSE){
  
  if(autorun == FALSE & is.null(gene_list) == T){stop("Gene list is expected!")}
  
  trained_pca_blood = blood_expr
  
  if(autorun == FALSE){
    if(is.null(gene_list)){stop("Please provide a list of genes that are present in new samples in order to run PCA correctly")}
    matched_genes = intersect(rownames(trained_pca_blood), gene_list)
    if(length(matched_genes) < 1){"No genes common to reference and new samples! Please check that IDs are in ENSEMBL gene ID format."}
    trained_pca_blood = trained_pca_blood[rownames(trained_pca_blood) %in% matched_genes, ]
  }
  
  return_these_genes = rownames(trained_pca_blood)
  
  pca_blood = prcomp(t(trained_pca_blood))
  
  return(list(pca = pca_blood, genes = return_these_genes))
}

# 4. predict PCA in new data
predict_pca = function(dat = NULL, pca_model = NULL, mean_imputation = FALSE){
  
  if(is.null(dat)){stop("Please supply data frame for new samples")}
  if(is.null(pca_model)){stop("Please supply fitted PCA model")}
  if(class(pca_model) != "list"){stop("Expecting list object for fitted PCA model")}
  
  new_samples = dat[,colnames(dat) %in% pca_model$genes]
  
  if(mean_imputation == TRUE){
    
    means = colMeans(new_samples,na.rm=TRUE)
    for(n in 1:ncol(new_samples)){
      new_samples[is.na(new_samples[,n]),n] = means[[n]]
    }
    
  }
  
  if(mean_imputation == FALSE){
    
    na_detected = colSums(is.na(new_samples))
    na_detected = na_detected[na_detected > 0]
    if(length(na_detected) > 0){stop("Warning! Missing values (NAs) detected in new samples. Switch mean_imputation to TRUE to resolve this issue.")}
    
  }
  
  preds = predict(pca_model$pca, newdata = new_samples)
  return(preds)
}


# 5. Fit LR-PCA prediction model to GTEx counts, then apply the resulting weights to the PCA scores derived from new samples
fit_lr_weights_in_gtex = function(pca_model  = NULL, tissue = NULL, gene_list = NULL, n_comps = 20){
  if(is.null(tissue)){stop("Please specify a tissue based on nomenclature available in the cross-validation table")}
  if(is.null(pca_model)){stop("Argument required for blood_PCA")}
  if(is.null(gene_list)){stop("Please specify gene IDs to include in the training process.")}
  # TODO: Linear regression: lung transcript ~ blood PCA
  Y = data.frame(t(lung_expr)) # select normalized GTEx lung counts
  # filter Y by genes that are well predicted from 5-fold cross-validation
  cv_performance = gene_list
  matched_genes = intersect(colnames(Y), cv_performance)
  if(length(matched_genes) < 1){stop("Unexpected issue! No matching gene IDs between CV performance file and normalized GTEx counts.")}
  common_genes = intersect(matched_genes, pca_model$genes)
  if(length(matched_genes) < 1){stop("Unexpected issue! No genes available to run LR.")}
  message("\rTraining LR models for: ", length(matched_genes), " genes")
  
  Y = Y[,colnames(Y) %in% matched_genes, ]
  X = data.frame(pca_model$pca$x[,1:n_comps]) # use PCA model derived from GTEx paired blood-lung data
  fit = lm(as.matrix(Y) ~ ., data = X) # fit a LR model per gene
  return(fit) # return model
}

fit_lr_weights_in_gtex_age_sex = function(pca_model  = NULL, tissue = NULL, gene_list = NULL, n_comps = 20){
  if(is.null(tissue)){stop("Please specify a tissue based on nomenclature available in the cross-validation table")}
  if(is.null(pca_model)){stop("Argument required for blood_PCA")}
  if(is.null(gene_list)){stop("Please specify gene IDs to include in the training process.")}
  # TODO: Linear regression: lung transcript ~ blood PCA
  Y = data.frame(t(lung_expr)) # select normalized GTEx lung counts
  # filter Y by genes that are well predicted from 5-fold cross-validation
  cv_performance = gene_list
  matched_genes = intersect(colnames(Y), cv_performance)
  if(length(matched_genes) < 1){stop("Unexpected issue! No matching gene IDs between CV performance file and normalized GTEx counts.")}
  common_genes = intersect(matched_genes, pca_model$genes)
  if(length(matched_genes) < 1){stop("Unexpected issue! No genes available to run LR.")}
  message("\rTraining LR models for: ", length(matched_genes), " genes")
  
  Y = Y[,colnames(Y) %in% matched_genes, ]
  X = data.frame(pca_model$pca$x[,1:n_comps]) # use PCA model derived from GTEx paired blood-lung data
  X = data.frame(load_predictors[,!colnames(load_predictors) %in% "SampleID"], X) # bind age, sex, and blood PCs
  fit = lm(as.matrix(Y) ~ ., data = X) # fit a LR model per gene
  return(fit) # return model
}

# 6. apply lung transcriptome prediction model to blood PCA
impute_gxp = function(pca_in_new_sample = NULL, trained_model = NULL, scale = TRUE){
  
  if(is.null(trained_model)){stop("Please specify object containing pre-trained models")}
  
  if(is.null(pca_in_new_sample)){stop("Please provide blood-based transcriptome PCs for new samples")}
  
  predict_vals = predict(trained_model, data.frame(pca_in_new_sample))
  
  predict_ci = predict(trained_model, data.frame(pca_in_new_sample), interval='confidence', level = 0.95)
  
  if(scale == TRUE){
    
    predict_vals = scale(predict_vals)
  }
  
  return(data.frame(predict_vals))
  
  imputed_ci <<- predict_ci
}

impute_gxp_with_age_sex = function(pca_in_new_sample = NULL, trained_model = NULL, scale = TRUE){
  
  if(is.null(trained_model)){stop("Please specify object containing pre-trained models")}
  
  if(is.null(pca_in_new_sample)){stop("Please provide blood-based transcriptome PCs for new samples")}
  
  predict_vals = predict(trained_model, data.frame(pca_in_new_sample))
  
  if(scale == TRUE){
    
    predict_vals = scale(predict_vals)
  }
  
  return(data.frame(predict_vals, check.names = F))
}


# additional functions:

# > convert hgnc symbols to ensembl gene ids (gencode v26)
convert_hgnc_to_ensg = function(hgnc = NULL, gtf_path = NULL){
  
  # check if parameters are filled in
  if(is.null(gtf_path)){stop("Please specify path to the gtf file in the LungGENIE provided repository")}
  if(is.null(hgnc)){stop("Please specify a vector of HGNC gene symbols")}
  
  # load GTF file
  gtf = data.frame(fread(gtf_path))
  gtf = gtf[,colnames(gtf) %in% c("gene_id", "gene_name")]
  gtf$gene_name = gsub("[-]", ".", gtf$gene_name)
  
  hgnc = gsub("[-]", ".", hgnc)
  
  symbols = data.frame(gene_name = hgnc)
  convert = merge(symbols, gtf, by='gene_name')
  return(convert)
}

convert_ensg_to_hgnc = function(ensg = NULL, gtf_path = NULL){
  
  # check if parameters are filled in
  if(is.null(gtf_path)){stop("Please specify path to the gtf file in the LungGENIE provided repository")}
  if(is.null(ensg)){stop("Please specify a vector of ensembl gene ids")}
  
  # load GTF file
  gtf = data.frame(fread(gtf_path))
  gtf = gtf[,colnames(gtf) %in% c("gene_id", "gene_name")]
  gtf$gene_name = gsub("[-]", ".", gtf$gene_name)
  
  ensg = gsub("[-]", ".", ensg)
  
  symbols = data.frame(gene_id = ensg)
  convert = merge(symbols, gtf, by='gene_id')
  return(convert)
}



# re-train PCA models using specified gene lists present in target sample:
retrain_gtex = function(gene_list = NULL, output = "", set_seed = T, seed = 123, ncomps = 20, prop_for_test_set = 0.1, n_folds = 5){
  
  if(base::exists("lung_expr") == FALSE | base::exists("blood_expr") == FALSE){stop("Warning! Paired blood-lung RNA-seq data have not been imported properly. Please run `load_expr_data` first.")}
  if(is.null(n_folds) == T){stop("Specify the number of cross-validation folds")}
  if(output == ""){warning("No output directory specified!")}
  
  # make output file path into string for loading CV performance
  output_file_string <<- paste(output, "/", tissue, "_LungGENIE_retrain-", ncomps, ".Rdata", sep="")
  
  use_n_comps = ncomps
  
  lung_counts = lung_expr
  blood_counts = blood_expr
  blood_counts = blood_counts[rownames(blood_counts) %in% gene_list, ]
  blood_counts = blood_counts[!rowSums(is.na(blood_counts)), ]
  
  # optional: define test set
  if(prop_for_test_set > 0){
    test.set = sample(colnames(lung_counts), ncol(lung_counts)*prop_for_test_set)
    
    test.lung = lung_counts[,colnames(lung_counts) %in% test.set]
    test.blood = blood_counts[,colnames(blood_counts) %in% test.set]
    
    # remove test set from cross-validation samples
    lung_counts = lung_counts[,!colnames(lung_counts) %in% test.set]
    blood_counts = blood_counts[,!colnames(blood_counts) %in% test.set]
  }
  
  # set up fold ids for cross-validation
  if(set_seed == T){seed = seed;
  
  # create folds for cross-validation
  create.folds = function(nfolds = 5, ids = NULL){
    set.seed(seed)
    if(is.null(ids)){stop("Expecting IDs column for creating folds")}
    sample(dplyr::ntile(ids, nfolds)) # random folds
  }
  }
  if(set_seed == F){
    if(set_seed == T){
      # create folds for cross-validation
      create.folds = function(nfolds = 5, ids = NULL){
        if(is.null(ids)){stop("Expecting IDs column for creating folds")}
        sample(dplyr::ntile(ids, nfolds)) # random folds
      }
    }
  }
  folds = create.folds(nfolds = n_folds, ids=colnames(lung_counts))
  
  allmods = list();
  pred.test.expr = list();
  for(fold in 1:max(folds)){
    
    cat("\rFold: ", fold);cat("\n")
    
    # define training/testing splits
    x_train = t(blood_counts[,which(folds != fold)])
    y_train = t(lung_counts[,which(folds != fold)])
    
    x_test = t(blood_counts[,!colnames(blood_counts) %in% rownames(x_train)])
    y_test = t(lung_counts[,!colnames(lung_counts) %in% rownames(y_train)])
    
    # - Run PCA on blood gene expression data
    pca.blood.train = prcomp(x_train)
    n.blood.comps = pca.blood.train$sdev^2 / sum(pca.blood.train$sdev^2)
    
    # apply PCA model to test set
    pca.in.test = predict(pca.blood.train, x_test)
    
    # format PCA predictor matrices based on chosen number of components 
    pca.TRAIN = data.frame(pca.blood.train$x[,1:use_n_comps])
    pca.TEST = data.frame(pca.in.test[,1:use_n_comps])
    
    
    # apply PCA model to predict gene expression in lung
    fitted.model = lm(as.matrix(y_train)  ~ ., data = pca.TRAIN)
    pred.train = predict(fitted.model)
    train.cors = lapply(1:ncol(pred.train), function(x) cor(pred.train[,x], y_train[,x]))
    names(train.cors) = colnames(pred.train)
    train.cors = ldply(train.cors, .id='gene')
    colnames(train.cors)[2] = 'train_cor'
    train.cors$train_rsq = train.cors$train_cor^2
    
    # apply PCA model to test set, evaluate accuracy
    pred.test = predict(fitted.model, newdata = pca.TEST)
    test.cors = lapply(1:ncol(pred.test), function(x) cor(pred.test[,x], y_test[,x]))
    names(test.cors) = colnames(pred.test)
    test.cors = ldply(test.cors, .id='gene')
    colnames(test.cors)[2] = 'test_cor'
    test.cors$test_rsq = test.cors$test_cor^2
    test.cors$N = nrow(pca.TEST)
    allmods[[fold]] = data.frame(cbind(train.cors, test.cors[,-1]))
    
    # save predicted expr for the validation fold
    pred.test.expr[[fold]] = data.frame(sampleid = rownames(pred.test), pred.test, check.names = F)
  }
  # combine all validation expr data into single d.f.
  names(pred.test.expr) = 1:max(folds)
  validation_expr <<- ldply(pred.test.expr, .id='fold')
  
  
  stats = ldply(allmods)
  
  # deal with missing values (NAs)
  stats$test_cor[is.na(stats$test_cor)] = 0
  stats$train_cor[is.na(stats$train_cor)] = 0
  stats$test_rsq[is.na(stats$test_rsq)] = 0
  stats$train_rsq[is.na(stats$train_rsq)] = 0
  
  # transform correlation coefficients 
  stats$fisherZ = atanh(stats$test_cor)*sqrt(stats$N-3)
  stats = data.table(stats)
  
  perf_train = stats[,list(Cor = tanh(mean(atanh(test_cor))),
                           train.cor = tanh(mean(atanh(train_cor))),
                           SE = sd(test_cor,  na.rm=TRUE)/sqrt(length(test_cor)),
                           SD = sd(test_cor, na.rm=T),
                           Rsq=mean(test_cor)^2,
                           FisherZ = sum(fisherZ),
                           sample.size = mean(N),
                           Zscore = sum(fisherZ) / sqrt(length(unique(folds)))),by=c("gene")]
  
  perf_train$pval = 2*pnorm(abs(perf_train$Zscore), lower.tail = FALSE) # calculate two-tailed p-value based on combined z-score
  perf_train$fdr = p.adjust(perf_train$pval, "fdr")
  
  if(prop_for_test_set > 0){
    # Run final test
    # define training/testing splits
    x_train = t(blood_counts)
    y_train = t(lung_counts)
    
    x_test = t(test.blood)
    y_test = t(test.lung)
    
    # - Run PCA on blood gene expression data
    pca.blood.train = prcomp(x_train)
    n.blood.comps = pca.blood.train$sdev^2 / sum(pca.blood.train$sdev^2)
    # use_n_comps = min(which(cumsum(n.blood.comps) >= pca.prop.min))
    
    
    # apply PCA model to test set
    pca.in.test = predict(pca.blood.train, x_test)
    
    # format PCA predictor matrices based on chosen number of components 
    pca.TRAIN = data.frame(pca.blood.train$x[,1:use_n_comps])
    pca.TEST = data.frame(pca.in.test[,1:use_n_comps])
    
    
    # apply PCA model to predict gene expression in lung
    fitted.model = lm(as.matrix(y_train)  ~ ., data = pca.TRAIN)
    pred.train = predict(fitted.model)
    train.cors = lapply(1:ncol(pred.train), function(x) cor(pred.train[,x], y_train[,x]))
    names(train.cors) = colnames(pred.train)
    train.cors = ldply(train.cors, .id='gene')
    colnames(train.cors)[2] = 'train_cor'
    train.cors$train_rsq = train.cors$train_cor^2
    
    # apply PCA model to test set, evaluate accuracy
    pred.test = predict(fitted.model, newdata = pca.TEST)
    test.cors = lapply(1:ncol(pred.test), function(x) cor.test(pred.test[,x], y_test[,x]))
    test.cors.est = unlist(lapply(test.cors, function(x) x$estimate))
    test.cors.pval = unlist(lapply(test.cors, function(x) x$p.value))
    names(test.cors.est) = colnames(pred.test)
    test.cors = ldply(test.cors.est, .id='gene')
    colnames(test.cors)[2] = 'test_cor'
    test.cors$test_pval = unlist(test.cors.pval)
    test.cors$test_rsq = test.cors$test_cor^2
    test.cors$N = nrow(pca.TEST)
    final.model = data.frame(cbind(train.cors, test.cors[,-1]))
    colnames(final.model)[-1] = paste("final_", colnames(final.model)[-1], sep="")
    
    
    colnames(perf_train)[-c(1)] = paste("cv_", colnames(perf_train)[-c(1)], sep="")
    all.models.merge = merge(perf_train, final.model, by='gene')
    all.models.merge$final_test_fdr = p.adjust(all.models.merge$final_test_pval ,"fdr")
    
  } else {
    all.models.merge = perf_train
  }
  
  all.models.merge$tissue = tissue
  
  # save output as .Rdata file
  saveRDS(all.models.merge, file=paste(output, "/", tissue, "_LungGENIE_retrain-",ncomps,".Rdata", sep=""))
  
  perf <<- all.models.merge
  
}

 