setwd("/novo2019/Temp_files/Rtest/sem")
library(lavaan)
library(semPlot)
library(corrplot)
library(psych)

#######################################################
#the value of the modification index (mi in the output) 
#is the expected decrease in the model χ2
#model fit and modification using modification index
#######################################################
fitted_model_search <- function(model_formula, data, bootstrap = FALSE){
  fit0 <- sem(model_formula, data=data)
  pval <- fitMeasures(fit0, c('pvalue'))
  if(pval > 0.05){
    cat(model_formula, "\n")
    return(list(fit.res = fit0, model = model_formula))
  }else{
    mi <- modindices(fit0)
    sel.mi <- mi[mi$op == '~~', ]
    max.mi.formula <- paste(sel.mi[which.max(sel.mi$mi), ][1:3], collapse = ' ')
    model_formula <- paste(model_formula, max.mi.formula, sep = "; ")
    temp.fit <- sem(model_formula, data=data, bootstrap = bootstrap)
    pval <- fitMeasures(temp.fit, c('pvalue'))
    if(pval > 0.05){
      cat(model_formula, "\n")
      return(list(fit.res = temp.fit, model = model_formula))
    }else{
      return(fitted_model_search(model_formula, data, bootstrap = bootstrap))
    }
  }
}

plot_matrix <- function(matrix_toplot){
  corrplot::corrplot(matrix_toplot, is.corr = FALSE,
                     type = 'lower',
                     order = "original",
                     tl.col='black', tl.cex=.75)
}

obs_vars_extract <- function(model_formula){
  return(gsub(' |\"', "", unlist(strsplit(unlist(strsplit(model_formula, '=~'))[-1],'\\+'))))
}
######################################################################
#                KMO test and Bartlett's test                        #
# fi_is_fit function is used for test whether the observed variabels #
# matrix is fit for factor analysis, two indicators, MSA (>0.6 is    #
# acceptable for fa) pvalue < 0.05 shows variables are correlated    #
######################################################################

########################################################################
#https://www.statisticshowto.datasciencecentral.com/kaiser-meyer-olkin/#
#KMO values between 0.8 and 1 indicate the sampling is adequate.	   #
#KMO values less than 0.6 indicate the sampling is not adequate 	   #
#and that remedial action should be taken. Some authors put this 	   #
#value at 0.5, so use your own judgment for values between 0.5 		   #
#and 0.6.															   #
#KMO Values close to zero means that there are large partial           #
#correlations compared to the sum of correlations. In other words,	   #
#there are widespread correlations which are a large problem for 	   #
#factor analysis.													   #
########################################################################

fa_is_fit <- function(measurement_model, data, n = 100){
  obs_vars <- obs_vars_extract(measurement_model)
  myidx <- (colnames(data) %in% obs_vars)
  select_dt <- data[, myidx]
  corr.dt <- psych::corr.test(select_dt)$r
  msa <- psych::KMO(corr.dt)$MSA
  cor.bar <- psych::cortest.bartlett(corr.dt, n = n)
  res <- data.frame(MSA = msa, Chisq = cor.bar$chisq, pvalue = cor.bar$p.value, df = cor.bar$df)
  return(res)
}

################################################################
#TEST DATA IN lavaan
################################################################
data(HolzingerSwineford1939)

### initial model format
initial_model_form <- ' visual  =~ x1 + x2 + x3;  
                        textual =~ x4 + x5 + x6; 
                        speed   =~ x7 + x8 + x9'

fa_is_fit('visual  =~ x1 + x2 + x3', data = HolzingerSwineford1939, n = 301)
fa_is_fit('textual =~ x4 + x5 + x6', data = HolzingerSwineford1939, n = 301)
fa_is_fit('speed   =~ x7 + x8 + x9', data = HolzingerSwineford1939, n = 301)



####best fitted model
#meanstructure = T， get the intercept
fit.res <- fitted_model_search(model_formula = initial_model_form, data=HolzingerSwineford1939)
#bootstrap can efficiently test the significance of coefficiencts
fit <- sem(model = fit.res$model, data = HolzingerSwineford1939, se = 'bootstrap')
#######
res <- summary(fit, fit.measures = T, standardized=T, rsq = T)

est_tab <- res$PE
est_tab$exo <- NULL
est_tab$std.nox <- NULL
############################
# parameterEstimates(fit) ##
# standardizedsolution(fit)#

# varTable(fit)
# parTable(fit)
# parameterTable(fit)
# fitMeasures(fit)
md.index <- fitMeasures(fit, c('chisq', 'df', 'pvalue', 'cfi', 'aic', 'bic', 'rmsea', 'rmsea.pvalue', 'srmr', 'gfi', 'mfi'))

out.index <- data.frame(index = names(md.index), value = md.index)

semPaths(fit, 
         layout = 'circle',
         what = 'stand', whatLabels = 'stand', # label the standardized path coefficient
         intStyle = 'single', # tree style, single or multi
         nCharNodes = 100, nCharEdges = 100,  #max lable string width for nodes and edges
         sizeMan = 8, sizeLat = 10, # size of manifested variable (observed) and latent variables
         #shapeMan = '', shapeLat = '',
         structural = T, # if only show the structure model(latent vars)
         residuals = F, # if show the residual covariance error matrix
         combineGroups = TRUE, ### if multiple groups, combine these figures to one 
         color = list('man' = 'lightblue', 'lat' = 'lightgreen'),
         edge.color = 'grey10',
         nDigits = 3,
         edge.label.cex = 1, 
         layoutSplit = F,
         posCol = 'red', negCol = 'blue',
         curvePivot = TRUE,
         edge.label.position = 0.5,
         bg = F,
)

semPaths(fit, 
         layout = 'tree',
         what = 'stand', whatLabels = 'stand', # label the standardized path coefficient
         intStyle = 'single', # tree style, single or multi
         nCharNodes = 100, nCharEdges = 100,  #max lable string width for nodes and edges
         sizeMan = 8, sizeLat = 10, # size of manifested variable (observed) and latent variables
         #shapeMan = '', shapeLat = '',
         structural = F, # if only show the structure model(latent vars)
         residuals = F, # if show the residual covariance error matrix
         combineGroups = TRUE, ### if multiple groups, combine these figures to one 
         color = list('man' = 'lightblue', 'lat' = 'lightgreen'),
         edge.color = 'grey10',
         nDigits = 3,
         edge.label.cex = 1, 
         layoutSplit = F,
         posCol = 'red', negCol = 'blue',
         edge.label.position = 0.5,
         bg = F,
)

#original model-implied covariance matrix
fitted(fit)
#standardized model-implied covariance matrix 
inspect(fit, what="cor.all")
#sample (observed) covariance matrix
lavCor(fit)
#bmisfit between observed- and model-implied matrices 
resid(fit, "cor")
#plot misfit residual covariance matirx
plot_matrix(resid(fit, type="cor")$cov)
