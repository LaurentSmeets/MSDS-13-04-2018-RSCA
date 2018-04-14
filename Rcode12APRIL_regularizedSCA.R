 library(RegularizedSCA)
 library(devtools)
 install_github("ZhengguoGu/RegularizedSCA")
 
 #Herring data: 2 blocks
 names(Herring)
 ChemPhy <- pre_process(Herring$Herring_ChemPhy)
 Sensory <- pre_process(Herring$Herring_Sensory)
 herring_data <- cbind(ChemPhy, Sensory)
 num_var <- cbind(dim(ChemPhy)[2], dim(Sensory)[2])
 
 #model selection: nr of components; status (common/distinctive)
 #1. vaf to decide nr of components
 vaf <- VAF(DATA = herring_data, Jk = num_var, R = 10)
 summary(vaf)
 #2. discosca to decide status
 discoresult <- DISCOsca(DATA = herring_data, R = 4, Jk = num_var)
 summary(discoresult)
 #3. pca-gca to decide number of components and status
 pca_gca(DATA = herring_data, Jk = num_var)
 
 #ANALYSIS WITH PRE-DEFINED STRUCTURE
 targetmatrix <- matrix(c(1, 1, 1, 1, 1, 0, 0, 1), nrow = 2, ncol = 4)
 maxLasso <- maxLGlasso(DATA = herring_data, num_var, R = 4)$Lasso
 set.seed(115)
 #introduce sparseness: lasso penalty
 #!!!!!!
 results_cvS <- cv_structuredSCA(DATA = herring_data, Jk = num_var, R = 4,
                Target = targetmatrix, Position = c(1, 2, 3, 4), LassoSequence = seq(from = 0.0000001, to = maxLasso, length.out = 200))
 #!!!!!!!
 plot(results_cvS)
 results_cvS$LassoRegion
 set.seed(115)
 result_str <- structuredSCA(DATA = herring_data, Jk = num_var, R = 4, Target =
                                 targetmatrix, Position = c(1, 2, 3, 4), LASSO = 0.8922256)
 final_comLoadingS <- undoShrinkage(DATA = herring_data, R = 4, Phat =
                                        result_str$Pmatrix)
 summary(final_comLoadingS)
 set.seed(115)
 #!!!!!!!
 system.time(results_cv <- cv_sparseSCA(DATA = herring_data, Jk = num_var, R = 4))
 #!!!!!!!
 summary(results_cv)
 set.seed(115)
 final_results <- sparseSCA(herring_data, num_var, R = 4, LASSO = 1.503148,
                              GROUPLASSO = 0.3655355, NRSTART = 20)
 final_Loading <- undoShrinkage(herring_data, R = 4, Phat = final_results$Pmatrix)
 summary(final_Loading)
 