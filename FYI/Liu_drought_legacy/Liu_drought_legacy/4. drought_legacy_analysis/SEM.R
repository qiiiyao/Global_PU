
library(lavaan)
library(semPlot)
#读取数据

data <-read.csv('K:\\Liu_drought_legacy\\sem\\sem_variables.csv')
names(data)
#删除缺失值
data <- na.omit(data)

#标准化变量
data_scaled <- scale(data)

mod3 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss
          ST.anomaly  ~  WT.anomaly'

fit3 <- sem(mod3, data = data_scaled)
fitMeasures(fit3, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))
semPaths(fit3,what="std",
         layout="tree",
         fade=F,
         residuals=F,
         nCharNodes=0, rotation = 2,fontsize =50
)




mod1 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss'
fit <- sem(mod1, data = data_scaled)


mod2 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly+SM.loss
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss'
fit2 <- sem(mod2, data = data_scaled)



mod3 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss
          ST.anomaly  ~  WT.anomaly'

fit3 <- sem(mod3, data = data_scaled)



mod4 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly+SM.loss
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss
          ST.anomaly  ~  WT.anomaly'

fit4 <- sem(mod4, data = data_scaled)


mod5 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss
          ST.anomaly  ~  WT.anomaly
         ST.anomaly  ~  SM.anomaly'
fit5 <- sem(mod5, data = data_scaled)



mod6 <- 'SOS_change~GPP.anomaly+EOS.anomaly+SM.anomaly+ST.anomaly+WT.anomaly+SM.loss
          GPP.anomaly~SM.loss
          EOS.anomaly~SM.loss
          SM.anomaly~SM.loss
          ST.anomaly~SM.loss                                                                                                                                                                                              
          WT.anomaly~SM.loss
          ST.anomaly  ~  WT.anomaly
          ST.anomaly  ~  SM.anomaly'
fit6 <- sem(mod6, data = data_scaled)
fitMeasures(fit, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))
fitMeasures(fit2, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))
fitMeasures(fit3, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))
fitMeasures(fit4, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))
fitMeasures(fit5, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))
fitMeasures(fit6, c('chisq', 'df', 'pvalue', 'cfi', 'gfi','rmsea', 'srmr', 'AIC'))


semPaths(fit3,what="std",
         layout="tree",
         fade=F,
         residuals=F,
         nCharNodes=0, rotation = 2,fontsize =50
)


semPaths(fit4,what="std",
         layout="tree",
         fade=F,
         residuals=F,
         nCharNodes=0, rotation = 2,fontsize =50
)


semPaths(fit5,what="std",
         layout="tree",
         fade=F,
         residuals=F,
         nCharNodes=0, rotation = 2,fontsize =50
)

semPaths(fit6,what="std",
         layout="tree",
         fade=F,
         residuals=F,
         nCharNodes=0, rotation = 2,fontsize =50
)


