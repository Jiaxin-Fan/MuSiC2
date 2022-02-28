# Utility Function
#
# Author: Jiaxin Fan
###################################################
#' @title MuSiC2_Deconvolution
#'
#' @description This function is used to deconvolve bulk RNA-seq data using single-cell reference generated under a different condition.
#' @param bulk.eset ExpressionSet for bulk data
#' @param sc.eset ExpressionSet for single cell data
#' @param condition character, the phenoData of bulk dataset used for indicating clinical conditions;
#' @param control character, the clinical condition of samples in bulk data that is the same as the clinical condition of samples in the single cell data;
#' @param case character, the clinical condition of samples in bulk data that is different from the clincial condition of samples in the single cell data;
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character, the phenoData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by the single cell dataset;
#' @param expr_low numeric, cutoff for defining lowly expressed genes in bulk data. Genes with mean expression across samples in bulk data < expr_low will be excluded from cell-type-specific DE gene detection. Default is 20;
#' @param prop_r numeric, cutoff for defining rare cell types. Cell types with mean proportion across samples in bulk data < prop_r will be characterized as rare cell types. Otherwise, will be characterized as common cell types. Default is 0.1;
#' @param eps_c numeric, convergence cutoff for common cell types. The cell type proportion estimate is converged if absolute relative change of proportion estimates for the current iteration against the previous iteration < eps_c. Default is 0.05;
#' @param eps_r numeric, convergence cutoff for rare cell types. The cell type proportion estimate is converged if absolute change of proportion estimates for the current iteration against the previous iteration < eps_r. Default is 0.01;
#' @param n_resample numeric, number of resamples used for detecting cell-type-specific DE genes. Default is 20;
#' @param sample_prop numeric, proportion of samples to be randomly sampled without replacement under each condition at each resampling. Default is 0.5;
#' @param cutoff_expr numeric, cutoff for defining lowly expressed genes over resamples. Genes with average cell-type-specific expression calculated over all resamples in the lower cutoff_expr quantile are excluded from cell-type-specific DE gene detection. Default is 0.05;
#' @param cutoff_c numeric, cutoff for defining cell-type-specific DE genes for common cell types. Genes with the value of statistic, defined as the absolute value of the ratio of the mean and standard deviation of the log fold change over all resamples, in the upper cutoff_c quantile are considered as cell-type-specific DE genes. Default is 0.05;
#' @param cutoff_r numeric, cutoff for defining cell-type-specific DE genes for rare cell types. Genes with the value of statistic, defined as the absolute value of the ratio of the mean and standard deviation of the log fold change over all resamples, in the upper cutoff_r quantile are considered as cell-type-specific DE genes. Default is 0.01;
#' @param maxiter numeric, maximum number of iterations. Default is 200;
#' @param markers vector or list of gene names. Default as NULL, i.e., use all genes that provided by both bulk and single cell datasets;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param centered logic, substract avg of Y and D;
#' @param normalize logic, divide Y and D by their standard deviation;
#' @return If MuSiC2 converges, return:
#' \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {n.iter: numeric, number of iterations.}
#'    \item {DE.genes: vector, cell-type-specific DE genes being removed.}
#'    }
#'  Or if MuSiC2 does not converge, return:
#'  \itemize{
#'    \item {Est.prop: matrix, cell type proportion estimates.}
#'    \item {convergence: logical, whether MuSiC2 converged or not.}
#'    \item {id.not.converge: vector, sample ids that failed to converge.}
#'    }
#' @seealso
#' \code{\link[MuSiC:music_prop]{music_prop}}
#' @export
#' @import MuSiC nnls xbioc

music2_prop = function(bulk.eset, sc.eset, condition, control,case, clusters, samples, select.ct, expr_low=20, prop_r=0.1, eps_c=0.05, eps_r=0.01, n_resample=20, sample_prop=0.5,cutoff_expr=0.05, cutoff_c=0.05, cutoff_r=0.01, maxiter = 200, markers = NULL, cell_size = NULL, ct.cov = FALSE, centered = FALSE, normalize = FALSE){
  gene_all = intersect(rownames(bulk.eset),rownames(sc.eset))
  bulk.eset = bulk.eset[gene_all,]
  sc.iter.eset = sc.eset

  # remove lowly expressed genes from DE analysis: i.e., gene with average expression < expr_low
  expr = exprs(bulk.eset)
  expr = apply(expr,1,mean)
  exp_genel = names(expr[expr>=expr_low])

  # separate heterozygous bulk samples based on their disease status
  control_sample = bulk.eset[,sampleNames(bulk.eset)[pData(bulk.eset)[ , condition] == control]]
  CASE_sample = bulk.eset[,sampleNames(bulk.eset)[pData(bulk.eset)[ , condition] == case]]

  # Step 1: cell type deconvolution, set initial value
  # estimate cell type proportion for controls using music.
  # tutorial for basic music see https://xuranw.github.io/MuSiC/articles/MuSiC.html
  prop_control = music_prop(bulk.eset = control_sample, sc.eset = sc.eset,
                            clusters = clusters, samples = samples, select.ct = select.ct,
                            markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                            nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_case_fix = NULL
  prop_case_ini = music_prop(bulk.eset = CASE_sample, sc.eset = sc.eset,
                              clusters = clusters, samples = samples, select.ct = select.ct,
                              markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                              nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize, verbose = F)$Est.prop.weighted
  prop_CASE = prop_case_ini
  prop_all = rbind(prop_control,prop_CASE)

  # start iteration
  iter=1
  ncell = length(select.ct)
  id_conv = NULL
  while(iter <= maxiter){
    # print(iter)
    # step 2: identify cell-type-specific DE genes
    # calculate mean/ sd of log fold change as an indicator of differential expression
    LOGFC = NULL
    MOD0=MOD1=matrix(0L, nrow = ncell, ncol = length(exp_genel))

    for(i in 1:n_resample){
      id_h = sample(colnames(control_sample),round(ncol(control_sample)*sample_prop))
      control_s = exprs(control_sample[exp_genel, colnames(control_sample) %in% id_h])
      prop_h = prop_control[colnames(control_s),]

      mod0 = apply(control_s, 1, function(x){
        mod = nnls(prop_h,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD0=MOD0+mod0

      id_d = sample(colnames(CASE_sample),round(ncol(CASE_sample)*sample_prop))
      case_s = exprs(CASE_sample[exp_genel, colnames(CASE_sample) %in% id_d])
      prop_d = prop_CASE[colnames(case_s),]

      mod1 = apply(case_s, 1, function(x){
        mod = nnls(prop_d,x)
        if(mod$mode==1){
          return(mod$x)
        }else{
          return(rep(0,ncell))
        }
      })
      MOD1=MOD1+mod1
      LOGFC = rbind(LOGFC,log1p(mod1)-log1p(mod0))
    }

    rcv=NULL
    for(i in 1:ncell){
      s = LOGFC[seq(from=i,to=nrow(LOGFC), by=ncell),]
      rcv = rbind(rcv,apply(s,2,function(x){ifelse(mean(x)==0, 0, mean(x)/sd(x))}))
    }
    abs_rcv_logfc = abs(rcv)
    MOD0 = MOD0/n_resample
    MOD1 = MOD1/n_resample
    rownames(MOD0)=rownames(MOD1)=rownames(abs_rcv_logfc)=select.ct

    # average cell type proportion
    mex = apply(prop_all,2,mean)
    lr = NULL
    for(celltype in select.ct){
      m = mex[celltype]
      rh = MOD0[celltype,]
      rd = MOD1[celltype,]
      # for genes with average expression within lower cutoff_expr for both conditions are removed from cell-type-specific DE genes detection
      llr = unique(intersect(names(rd[rd <= quantile(rd,prob=cutoff_expr)]),
                             names(rh[rh <= quantile(rh,prob=cutoff_expr)])))
      x = abs_rcv_logfc[celltype,]
      x = x[!names(x) %in% llr]
      # select genes with large mean/cv of log fc as DE
      if(m >= prop_r){
        lr = c(lr, names(x[x >= quantile(x,prob=1-cutoff_c)]))
      }else{
        lr = c(lr, names(x[x >= quantile(x,prob=1-cutoff_r)]))
      }
    }
    lr = unique(lr)

    # step 3: update sc gene list
    # remove identified DE genes from sc rna-seq data
    l = setdiff(gene_all,lr)
    sc.iter.eset = sc.eset[l,]

    # step 1: update cell type proportion based on new gene list
    if(length(id_conv)>0){
      case_sample = CASE_sample[ , !colnames(CASE_sample) %in% id_conv]
    }else{
      case_sample = CASE_sample
    }

    prop_case = music_prop(bulk.eset = case_sample, sc.eset = sc.iter.eset,
                           clusters=clusters, samples=samples, select.ct=select.ct,
                           markers = markers, cell_size = cell_size, ct.cov = ct.cov, iter.max = 1000,
                           nu = 0.0001, eps = 0.01, centered = centered, normalize = normalize,verbose = F)$Est.prop.weighted

    prop_CASE = rbind(prop_case,prop_case_fix)

    if(length(id_conv)==1){
      rownames(prop_CASE) = c(rownames(prop_case),id_conv)
    }

    prop_all = rbind(prop_control,prop_CASE)

    # check convergence, by cell type
    prop_case=prop_case[rownames(prop_case_ini),]
    pc = abs(prop_case-prop_case_ini)
    conv = pc
    conv[,] = 1
    # use difference if rare cell type
    conv[prop_case_ini <= prop_r] = ifelse(pc[prop_case_ini <= prop_r] < eps_r, 0, 1)
    # use percent change if common cell type
    pc[prop_case_ini > prop_r] = pc[prop_case_ini>prop_r]/prop_case_ini[prop_case_ini>prop_r]
    conv[prop_case_ini > prop_r] = ifelse(pc[prop_case_ini > prop_r] < eps_c,0,1)
    convf = apply(conv,1,function(x){all(x==0)})

    # if an id converged, not updating anymore
    all_converge=FALSE
    id_conv = c(id_conv,names(convf[convf==TRUE]))
    prop_case_ini = prop_CASE[!rownames(prop_CASE) %in% id_conv,]
    prop_case_fix = prop_CASE[rownames(prop_CASE) %in% id_conv,]

    # if all converged or if only one subjects not converging--> music2 converged
    if(is.vector(prop_case_ini)){
      all_converge=TRUE
      break
    }else if(nrow(prop_case_ini)==0){
      all_converge=TRUE
      break
    }
    iter=iter+1
  }
  # return estimated proportion
  if(all_converge){
    return(list('Est.prop' = prop_all,'convergence'=TRUE,'n.iter'=iter,'DE.genes'=lr))}
  else{
    return(list('Est.prop' = prop_all,'convergence'=FALSE,'id.not.converge'=rownames(prop_case_ini)))}
}

