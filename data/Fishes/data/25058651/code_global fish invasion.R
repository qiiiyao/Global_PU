library(tidyverse)
library(conflicted)
library(FishPhyloMaker)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ape)
library(picante)
library(lme4)
library(lmerTest)
library(performance)
library(nlme)
library(INLA)
library(sjPlot)
library(sf)
library(spData)
library(ggsci)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggthemes)
library(scico)
library(circlize)
library(piecewiseSEM)
library(brms)

#读入数据
drainage_basins = read.csv(file.choose(), sep = ";")
occ_drainage = read.csv(file.choose(), sep = ";")

#1.The generation of valid species names
occ_drainage$X2.Species.Name.in.Source = gsub("[.]", "_", occ_drainage$X2.Species.Name.in.Source)
occ_drainage$X6.Fishbase.Valid.Species.Name = gsub("[.]", "_", occ_drainage$X6.Fishbase.Valid.Species.Name)

species_list = unique(occ_drainage$X6.Fishbase.Valid.Species.Name)
all_taxa_names = FishTaxaMaker(data = species_list, allow.manual.insert = TRUE)
Cyprinidae
Cypriniformes
Loricariidae
Siluriformes
Pseudopimelodidae
Siluriformes
Petromyzontidae
Petromyzontiformes
Nemacheilidae
Cypriniformes
belontiidae
perciformes
Loricariidae
Siluriformes
Trichomycteridae
Siluriformes
Poeciliidae
Cyprinodontiformes

#2.Further cleaning of the valid species names file Taxon_data_FishPhyloMaker
all_taxa_names$Taxon_data_FishPhyloMaker$o  =  gsub("/.*", replacement = " ", all_taxa_names$Taxon_data_FishPhyloMaker$o) # removing slash in Perciformes
which(is.na(all_taxa_names$Taxon_data_FishPhyloMaker$s))#查找species为NA的行
all_taxa_names$Taxon_data_FishPhyloMaker = all_taxa_names$Taxon_data_FishPhyloMaker[-810,]#去除species为NA的行，这个很特殊，都搞定了怎么还有一个NA，remove它


#3.Details of the process for determining valid species names
#（1）提取All_info_fishbase 里 valid_names的61个重复物种名
all_taxa_names$All_info_fishbase[all_taxa_names$All_info_fishbase$valid_names%>%duplicated(),]%>%drop_na(valid_names)%>%pull(valid_names)%>%write.csv(file="duplicated names.csv")

# (2) 提取没有被辨别出来的9个物种名，这些名字非NA，但又没有出现在后面的Taxon_data_FishPhyloMaker的物种名里面
all_taxa_names$All_info_fishbase[!(gsub(" ","_",all_taxa_names$All_info_fishbase$valid_names)%in%all_taxa_names$Taxon_data_FishPhyloMaker$s),]%>%drop_na(valid_names)%>%pull(valid_names)%>%write.csv(file="not identified names.csv")

#（3）被辨别出来为无效物种名NA的9个物种，这些种通过FishTaxaMaker(allow.manual.insert=T), 给科名和目名，成功变为有效物种名。提取出物种名及给定的科和目名。
all_taxa_names$Taxon_data_FishPhyloMaker[all_taxa_names$Taxon_data_FishPhyloMaker$s %in% gsub(" ","_",all_taxa_names$Species_not_in_Fishbase),]%>%write.csv(file="na to valid names.csv")


#4.Generate valid species names for the species not identified in 3.(2), and integrate them with the Taxon_data_FishPhyloMaker.
not_identified_names = all_taxa_names$All_info_fishbase[!(gsub(" ","_",all_taxa_names$All_info_fishbase$valid_names)%in%all_taxa_names$Taxon_data_FishPhyloMaker$s),]%>%drop_na(valid_names)%>%pull(valid_names)
not_identified_names = gsub(" ","_",not_identified_names)
not_identified_taxa = FishTaxaMaker(not_identified_names)
not_identified_df = not_identified_taxa$Taxon_data_FishPhyloMaker#没有被辨别的有效物种名数据框
generated_df = all_taxa_names$Taxon_data_FishPhyloMaker#上面自动产生的有效物种名数据框

df_for_phylomaker = rbind(not_identified_df,generated_df)#共14892个物种

#5.Construction of a global fish phylogenetic tree
phylo_all_spp  =  FishPhyloMaker(data = df_for_phylomaker, 
                                insert.base.node = TRUE, 
                                return.insertions = TRUE)

#最终构建的树包含了14708个物种，因为有184个物种虽然是有效名但是FishPhyloMaker没有成功插入这些种到fish tree of life骨架树中。
phylo_all_spp$Insertions_data%>%filter(insertions=="Not_inserted")%>%pull(s)%>%write.csv(file="not inserted to fish tree.csv")

#提取出谱系树
phylo = phylo_all_spp$Phylogeny

#6.Calculation of pairwise phylogenetic distances
library(ape)
distance_all = cophenetic.phylo(phylo)
distance_all = distance_all/max(distance_all)
tibble(distance_all)

#7.Matching occurrence data with phylogenetic data

#更新数据中物种名为valid_names,使得与phylo的物种名保持一致
phylo.valid.names.df = all_taxa_names$All_info_fishbase%>%select(1:2)%>%tibble()
phylo.valid.names.df$valid_names = gsub(" ","_",phylo.valid.names.df$valid_names)
phylo.valid.names.df_na = phylo.valid.names.df%>%filter(is.na(valid_names))%>%mutate(valid_names=user_spp)#有效物种名有9个没有辨识出来为NA,建树时已经加上，这里也要加上
phylo.valid.names.df_residue = phylo.valid.names.df%>%filter(!is.na(valid_names))
phylo.valid.names.df.final = rbind(phylo.valid.names.df_residue,phylo.valid.names.df_na)

data.occurrence = left_join(occ_drainage,phylo.valid.names.df.final,by=c("X6.Fishbase.Valid.Species.Name"="user_spp"))%>%tibble()#有效物种名加入到数据中了

data.all = left_join(data.occurrence,drainage_basins)%>%select(X2.Country,everything())#把发生数据和流域信息数据结合

#排除没有exotic种的国家
country.used = data.all%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X2.Country)%>%unique()#143个国家中只有121个国家有外来鱼类
data.used = data.all%>%filter(X2.Country%in%country.used)%>%
  group_by(X1.Basin.Name)%>%
  distinct(valid_names,.keep_all = T)%>%
  select(c(1,2,4,9))%>%
  arrange(X2.Country,desc(X3.Native.Exotic.Status),X1.Basin.Name)%>%
  ungroup()#只用有外来鱼类的国家.😄每个流域内去除有效物种名后的重复物种.包含了14803个物种

#排除不在phylo中的物种,共183个物种,其中包括在4个国家被定义为exotic的6个物种(这6个种都是流域间转移种，因为检查发现一个国家内其他流域定义它们为本地种)
species.exclude = (data.used$valid_names%>%unique)[!((data.used$valid_names%>%unique) %in% phylo$tip.label)]
data.exclude = data.used%>%filter(valid_names%in%species.exclude)
data.exclude%>%summarise(n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.status=n_distinct(X3.Native.Exotic.Status),n.species=n_distinct(valid_names))#总结排除的物种
data.exclude%>%filter(X3.Native.Exotic.Status=="exotic")#总结4个国家被定义为exotic的6个物种

data.used = data.used%>%filter(!valid_names%in%species.exclude)
data.used$valid_names%>%unique()%>%length()#包含了14620个物种
data.used%>%print(n=100)

#纠正一个外来状态的错误，中国的Cyprinus_carpio鲤鱼被定义为exotic,这是不对的，将中国内鲤鱼的Native.Exotic.Status改为native,共54个数据
data.used = data.used%>%mutate(X3.Native.Exotic.Status=if_else(X2.Country=="China"&valid_names=="Cyprinus_carpio","native",X3.Native.Exotic.Status))

#--------------------Analyze the introduced alien species

#8.Conduct analysis using data on introduced alien species,excluding translocated species

#找出一个国家中流域间转移的外来种，即在一个流域内被定义为exotic但在其他流域是本地种
data.distinct = data.used %>% group_by(X2.Country,X3.Native.Exotic.Status) %>%
  distinct(valid_names,.keep_all = T)#一个国家内外来和本地分别独一无二的种
data.distinct = data.distinct%>%ungroup()
data.distinct%>%print(n=100)

data.dup = data.distinct %>% group_by(X2.Country) %>% 
  mutate(dup=duplicated(valid_names)) %>% 
  filter(dup)#一个国家内，外来和本地重复的种，这些种应该都是流域间转移种
data.dup = data.dup%>%unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)
data.dup%>%print(n=100)

#从data.used中将这些流域间转移种的属性exotic改为native,从而在国外引入计算exotic-native谱系距离时不考虑这些种为exotic
data.replace = data.used%>%
  unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  filter(new.col %in% data.dup$new.col)%>%
  mutate(X3.Native.Exotic.Status="native")
data.replace

data.residue = data.used%>%
  unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  filter(!(new.col %in% data.dup$new.col))
data.residue

data.used.combine = rbind(data.residue,data.replace)%>%
  arrange(X2.Country,X3.Native.Exotic.Status,X1.Basin.Name)%>%
  mutate(new.col=NULL)
data.used.combine

data.used.combine%>%group_by(X2.Country,X3.Native.Exotic.Status)%>%
  distinct(valid_names,.keep_all = T)%>%ungroup()%>%
  group_by(X2.Country)%>%mutate(dup=duplicated(valid_names))%>%filter(dup)#验证一下，确实国家内没有外来和本地重复了，流域间转移的外来种已全部更改为native

country.used.final = data.used.combine%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X2.Country)%>%unique()#改变后只有118个国家有exotic了
data.used.final = data.used.combine%>%filter(X2.Country%in%country.used.final)#过滤出这些有外来物种的国家
data.used.final = data.used.final%>%arrange(X2.Country,X1.Basin.Name,valid_names)
data.used.final#最终使用的dataframe
data.used.final%>%print(n=200)

#总结这个最终使用的dataframe
data.used.final%>%summarise(n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))
data.used.final%>%group_by(X3.Native.Exotic.Status)%>%summarise(n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))

#9.Generate a data frame combining basin and species in pairs
data.used.final_add.occur = data.used.final%>%
  mutate(occurrence=1,.after=valid_names)%>%
  arrange(X2.Country,X1.Basin.Name,X3.Native.Exotic.Status)#添加一列occurrence

data.expand = data.used.final%>%
  unite(new.col,X3.Native.Exotic.Status,valid_names,sep="/")%>%
  group_by(X2.Country)%>%
  tidyr::expand(X1.Basin.Name,new.col)%>%
  separate(new.col,into=c("X3.Native.Exotic.Status","valid_names"),sep="/")#产生basin和species两两结合的数据。expand()好啊，会把数据框两列中独一无二的数据两两结合

data.pair = left_join(data.expand,data.used.final_add.occur)%>%
  replace_na(list(occurrence=0))%>%
  arrange(X2.Country,X1.Basin.Name,X3.Native.Exotic.Status,valid_names)

data.pair%>%print(n=200)

#10.Calculate the phylogenetic distances between all alien species and all native species in each basin within a country
phylodist.function = function(data){
  
  data.out = NULL
  
  for(i in 1:length(unique(data$X1.Basin.Name))){
    
    basin = filter(data,X1.Basin.Name==unique(data$X1.Basin.Name)[i])
    
    exo = filter(basin,X3.Native.Exotic.Status=="exotic")%>%pull(valid_names)
    
    nat = filter(basin,X3.Native.Exotic.Status=="native"&occurrence!=0)%>%pull(valid_names)
    
    data.out.basin = list()
    
    for(j in 1:length(exo)){
      
      dist = distance_all[exo[j],nat]
      
      data.out.basin[[j]] = c(filter(basin,valid_names==exo[j])[1,1], filter(basin,valid_names==exo[j])[1,2], exo[j], filter(basin,valid_names==exo[j])[1,5],mean(dist,na.rm=T),min(dist,na.rm=T))
    }
    
   data.out = append(data.out,data.out.basin) 
    
  }
  
  data.final = as.data.frame(do.call(rbind,data.out))
  names(data.final) = c("country","basin","exotic_species","occurrence","mpd","mntd")
  data.final
}

data.distance = phylodist.function(data.pair)

data.distance.final = data.distance%>%mutate(country=unlist(country),basin=unlist(basin),exotic_species=unlist(exotic_species),occurrence=unlist(occurrence),mpd=unlist(mpd),mntd=unlist(mntd))%>%tibble()#把每个成分的列表转为向量
data.distance.final = data.distance.final%>%filter(!is.na(mpd))#有三行mpd为NaN,mntd为Inf,排除它们。Wadi.Libya的Wadi.Nashu就3个国外引入种，没有本地种😌😌

#11. Statistical tests
glmer_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd)
tab_model(glmer_mpd)

#手动画图
get.data_mpd = get_model_data(glmer_mpd,type="pred",terms="mpd [all]")
tibble(get.data_mpd)
(plot.mpd = 
  ggplot(data=get.data_mpd,aes(x=x,y=predicted))+
  geom_line(color="#00468BFF",linewidth=1)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#00468BFF",alpha=0.2)+
  annotate(geom="text",x=c(0.527,0.506,0.486,0.545,0.561,0.490),y=c(0.5,0.455,0.405,0.355,0.305,0.255),label=c("italic(β)[mpd]==-5.27","italic(z)==-15.61","italic(P)<2e-16","italic(R^2)[marginal]==0.02","italic(R^2)[conditional]==0.81","italic(n)==61090"),parse=T,size=3.5)+
  annotate(geom="text",x=0.5,y=0.6,label="Exotic fish species",fontface="bold")+
  labs(x=NULL,y="Occurrence probability")+
  scale_y_continuous(labels=function(x) sprintf("%.2f",x))+
  theme_classic()
)

glmer_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd)
tab_model(glmer_mntd)

#手动画图
get.data_mntd = get_model_data(glmer_mntd,type="pred",terms="mntd [all]")
tibble(get.data_mntd)
(plot.mntd = 
  ggplot(data=get.data_mntd,aes(x=x,y=predicted))+
  geom_line(color="#00468BFF",linewidth=1)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#00468BFF",alpha=0.2)+
  annotate(geom="text",x=c(0.527,0.506,0.486,0.545,0.545,0.490),y=c(0.215,0.195,0.175,0.155,0.135,0.114),label=c("italic(β)[mntd]==-4.37","italic(z)==-23.14","italic(P)<2e-16","italic(R^2)[marginal]==0.04","italic(R^2)[conditional]==0.80","italic(n)==61090"),parse=T,size=3.5)+
  annotate(geom="text",x=0.5,y=0.2615,label="Exotic fish species",fontface="bold")+
  labs(x=NULL,y=NULL)+
  scale_y_continuous(labels=function(x) sprintf("%.2f",x))+
  theme_classic()
)

#合并两个回归图,图片大小7.53*4.05
ggarrange(
  plot.mpd,plot.mntd,
  labels=c("(a)","(b)"),
  hjust=-1,vjust=1.2
)

#--------------------Analyze the translocated species

#12.Conduct analysis using data on translocated species,excluding alien species。

#找出国外引入的外来种
data.foreign = data.used.final%>%filter(X3.Native.Exotic.Status=="exotic")%>%
unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)

#从data.used中将这些国外引入的外来种的属性exotic改为native,从而在流域间转移计算exotic-native谱系距离时不考虑这些种为exotic
data.replace.translocation = data.used%>%
  unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  filter(new.col %in% data.foreign$new.col)%>%
  mutate(X3.Native.Exotic.Status="native")
data.replace.translocation

data.residue.translocation = data.used%>%
  unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  filter(!(new.col %in% data.foreign$new.col))
data.residue.translocation

data.used.combine.translocation = rbind(data.residue.translocation,data.replace.translocation)%>%
  arrange(X2.Country,X3.Native.Exotic.Status,X1.Basin.Name)%>%
  mutate(new.col=NULL)
data.used.combine.translocation

country.used.final.translocation = data.used.combine.translocation%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X2.Country)%>%unique()#改变后只有53个国家有转移exotic
data.used.final.translocation = data.used.combine.translocation%>%filter(X2.Country%in%country.used.final.translocation)#过滤出这些有转移外来物种的国家
data.used.final.translocation = data.used.final.translocation%>%arrange(X2.Country,X1.Basin.Name,valid_names)
data.used.final.translocation%>%print(n=200)

#总结这个最终使用的dataframe
data.used.final.translocation%>%summarise(n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))
data.used.final.translocation%>%group_by(X3.Native.Exotic.Status)%>%summarise(n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))

#13.Generate a data frame combining basin and species in pairs
data.used.final.translocation_add.occur = data.used.final.translocation%>%
  mutate(occurrence=1,.after=valid_names)%>%
  arrange(X2.Country,X1.Basin.Name,X3.Native.Exotic.Status)#添加一列occurrence

data.expand.translocation = data.used.final.translocation%>%
  unite(new.col,X3.Native.Exotic.Status,valid_names,sep="/")%>%
  group_by(X2.Country)%>%
  tidyr::expand(X1.Basin.Name,new.col)%>%
  separate(new.col,into=c("X3.Native.Exotic.Status","valid_names"),sep="/")#产生basin和species两两结合的数据。expand()好啊，会把数据框两列中独一无二的数据两两结合

data.pair.translocation = left_join(data.expand.translocation,data.used.final.translocation_add.occur)%>%
  replace_na(list(occurrence=0))%>%
  arrange(X2.Country,X1.Basin.Name,X3.Native.Exotic.Status,valid_names)

data.pair.translocation%>%print(n=200)

#14.Calculate the phylogenetic distances between all translocated species and all native species in each basin within a country
data.distance.translocation = phylodist.function(data.pair.translocation)

data.distance.translocation.final = data.distance.translocation%>%mutate(country=unlist(country),basin=unlist(basin),exotic_species=unlist(exotic_species),occurrence=unlist(occurrence),mpd=unlist(mpd),mntd=unlist(mntd))%>%tibble()#把每个成分的列表转为向量
data.distance.translocation.final = data.distance.translocation.final%>%filter(!is.na(mpd))#有1行mpd为NaN,mntd为Inf,排除它们
data.distance.translocation.final = data.distance.translocation.final%>%filter(mntd!=0)#这个很重要，这是忽略每个国家内转移外来种出去的流域和外来种的谱系距离（数据两两配对后，流域内同种native发生，exotic没有发生的流域，即是转移外来种出去的流域。因为exotic-native同种，所以计算出来的外来种与这个流域的最近谱系距离为0，忽略mntd=0 这些行的数据即可），因为只是要算潜在接收外来种的流域与外来种的距离

#15. Statistical tests
glmer_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd.translocation)
tab_model(glmer_mpd.translocation)

#手动画图
get.data_mpd.translocation = get_model_data(glmer_mpd.translocation,type="pred",terms="mpd [all]")
tibble(get.data_mpd.translocation)
(plot.mpd.translocation = 
    ggplot(data=get.data_mpd.translocation,aes(x=x,y=predicted))+
    geom_line(color="#f47920",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#f47920",alpha=0.2)+
    annotate(geom="text",x=c(0.528,0.493,0.524,0.546,0.56,0.486),y=c(0.095,0.085,0.075,0.065,0.055,0.045),label=c("italic(β)[mpd]==-3.81","italic(z)==-7.33","italic(P)==2.27e-13","italic(R^2)[marginal]==0.01","italic(R^2)[conditional]==0.86","italic(n)==62951"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.11,label="Translocated fish species",fontface="bold")+
    labs(x="Nonnative-native MPD",y="Occurrence probability")+
    theme_classic()
)

glmer_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd.translocation)
tab_model(glmer_mntd.translocation)

#手动画图
get.data_mntd.translocation = get_model_data(glmer_mntd.translocation,type="pred",terms="mntd [all]")
tibble(get.data_mntd.translocation)
(plot.mntd.translocation = 
    ggplot(data=get.data_mntd.translocation,aes(x=x,y=predicted))+
    geom_line(color="#f47920",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#f47920",alpha=0.2)+
    annotate(geom="text",x=c(0.527,0.506,0.486,0.545,0.561,0.490),y=c(0.047,0.042,0.037,0.032,0.027,0.022),label=c("italic(β)[mntd]==-3.55","italic(z)==-13.66","italic(P)<2e-16","italic(R^2)[marginal]==0.03","italic(R^2)[conditional]==0.85","italic(n)==62951"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.0542,label="Translocated fish species",fontface="bold")+
    labs(x="Nonnative-native MNTD",y=NULL)+
    theme_classic()
)

#合并两个回归图,图片大小7.53*4.05
ggarrange(
  plot.mpd.translocation,plot.mntd.translocation,
  labels=c("(c)","(d)"),
  hjust=-0.5,vjust=1.2
)

#合并国外引入和流域间转移的回归图.大小6.9*6.28.命名logistic_globe，使用😄。
ggarrange(
  plot.mpd,plot.mntd,
  plot.mpd.translocation,plot.mntd.translocation,
  labels=c("a","b","c","d"),
  ncol=2,nrow=2,
  font.label=list(size=16),
  hjust=-1,vjust=1.2
)


#16.Create a global sampling map, pie charts, and bar charts

#(1)地图(原始数据)。展示全球外来种的分布格局和占比情况，所有国家和流域都包括。

library(sf)
library(rnaturalearth)
W  =  ne_countries(scale = 50, returnclass = "sf")#世界地图

data.all1 = data.all%>%mutate(X3.Native.Exotic.Status=if_else(X2.Country=="China"&valid_names=="Cyprinus_carpio","native",X3.Native.Exotic.Status))#用原始数据data.all看流域及外来物种数，但需要修正下中国的鲤鱼为本地种(统计用数据data.used 已修正)，方便后面正确画出有外来鱼类流域的比例饼图，所以用data.all1
data.map1 = data.all1%>%select(c(1,2,10,14,15,16))%>%distinct(X2.Country,X1.Basin.Name,X3.Ecoregion,.keep_all = T)
map2 = ggplot(W)+#大小9.25*3.5
  geom_sf(fill="grey50",linewidth=0.1,color="white")+
  coord_sf(ylim = c(-50, 80))+
  geom_point(data=data.map1,aes(x=X7.Median.Longitude,y=X8.Median.Latitude,size=X9.Surface.Area,fill=X3.Ecoregion),shape=21,stroke=0.1,color="grey50")+
  scale_size(range=c(0.5,10))+
  scale_fill_npg(name="Biogeographical realms")+
  theme_map(base_size=11)+
  theme(legend.position=c(-0.01,0.10),legend.background = element_blank())+
  guides(size="none",fill=guide_legend(override.aes=list(size=3)))

map2

#(2)总结地图,画饼图。插入到地图(原始数据)中
#全球
data.all1%>%pull(X2.Country)%>%unique%>%length()#国家个数143
data.all1%>%pull(X1.Basin.Name)%>%unique%>%length()#流域数3119
data.all1%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length#有外来物种的流域数1719  55.11%
data.used.final%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#有国外引入外来物种的流域数1581  50.69%
data.used.final.translocation%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#有流域间转移外来物种的流域数603  19.33%

data.pie = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
           counts=c(1719,3119-1719))%>%arrange(desc(Basins))
data.pie = mutate(data.pie,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=9,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines")) 

#古北区
data.all1%>%filter(X3.Ecoregion=="Palearctic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区流域数1033
data.all1%>%filter(X3.Ecoregion=="Palearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区有外来物种的流域数649  62.83%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Palearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区有国外引入外来物种的流域数598  57.89%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Palearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区有流域间转移外来物种的流域数219  21.20%

data.pie.palearctic = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                     counts=c(649,1033-649))%>%arrange(desc(Basins))
data.pie.palearctic = mutate(data.pie.palearctic,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.palearctic,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=11,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines"))

#非洲热带区
data.all1%>%filter(X3.Ecoregion=="Afrotropic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区流域数282
data.all1%>%filter(X3.Ecoregion=="Afrotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区有外来物种的流域数100  35.46%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Afrotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区有国外引入外来物种的流域数92  32.62%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Afrotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区有流域间转移外来物种的流域数35  12.41%

data.pie.afrotropic = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                                counts=c(100,282-100))%>%arrange(desc(Basins))
data.pie.afrotropic = mutate(data.pie.afrotropic,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.afrotropic,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=11,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines"))

#新热带区
data.all1%>%filter(X3.Ecoregion=="Neotropic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区流域数426
data.all1%>%filter(X3.Ecoregion=="Neotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区有外来物种的流域数169  39.67%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Neotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区有国外引入外来物种的流域数152  35.68%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Neotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区有流域间转移外来物种的流域数57  13.38%

data.pie.neotropic = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                                counts=c(169,426-169))%>%arrange(desc(Basins))
data.pie.neotropic = mutate(data.pie.neotropic,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.neotropic,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=11,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines"))

#澳新区
data.all1%>%filter(X3.Ecoregion=="Australasia")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区流域数725
data.all1%>%filter(X3.Ecoregion=="Australasia"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区有外来物种的流域数397  54.76%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Australasia"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区有国外引入外来物种的流域数388  53.52%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Australasia"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区有流域间转移外来物种的流域数93  12.83%

data.pie.australasia = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                                counts=c(397,725-397))%>%arrange(desc(Basins))
data.pie.australasia = mutate(data.pie.australasia,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.australasia,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=11,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines"))

#东洋区
data.all1%>%filter(X3.Ecoregion=="Indo-Malay")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区流域数340
data.all1%>%filter(X3.Ecoregion=="Indo-Malay"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区有外来物种的流域数234  68.82%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Indo-Malay"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区有国外引入外来物种的流域数220  64.71%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Indo-Malay"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区有流域间转移外来物种的流域数48  14.12%

data.pie.indomalay = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                                counts=c(234,340-234))%>%arrange(desc(Basins))
data.pie.indomalay = mutate(data.pie.indomalay,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.indomalay,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=11,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines"))

#新北区
data.all1%>%filter(X3.Ecoregion=="Nearctic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区流域数296
data.all1%>%filter(X3.Ecoregion=="Nearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区有外来物种的流域数164  55.41%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Nearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区有国外引入外来物种的流域数125  42.23%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Nearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区有流域间转移外来物种的流域数151  51.01%

data.pie.nearctic = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                               counts=c(164,296-164))%>%arrange(desc(Basins))
data.pie.nearctic = mutate(data.pie.nearctic,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.nearctic,aes(x="",y=prop,fill=Basins))+#大小15*6.28
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=15,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"),labels=c("River basins with non-native fishes","River basins without non-native fishes"))+
  guides(fill=guide_legend(keywidth=4,keyheight=4))+
  theme_void()+
  theme(legend.text=element_text(size=33))+
  theme(plot.margin = unit(c(0,0.2,0,-1.8),units="lines"))

#大洋区
data.all1%>%filter(X3.Ecoregion=="Oceania")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区流域数17
data.all1%>%filter(X3.Ecoregion=="Oceania"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区有外来物种的流域数6  35.29%
data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Oceania"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区有国外引入外来物种的流域数6  35.29%
data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Oceania"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区有流域间转移外来物种的流域数0  0%

data.pie.oceania = data.frame(Basins=c("Drainage basins with exotic fishes","Drainage basins without exotic fishes"),
                                counts=c(6,17-6))%>%arrange(desc(Basins))
data.pie.oceania = mutate(data.pie.oceania,prop=round(100*counts/sum(counts),2),position=cumsum(prop)-0.5*prop)

ggplot(data.pie.oceania,aes(x="",y=prop,fill=Basins))+#大小4*4
  geom_col(color="white",linewidth=3)+
  geom_text(aes(y=position,label=paste(prop,"%")),size=11,fontface="bold")+
  coord_polar(theta="y")+
  scale_fill_manual(name="",values=c("#fdb933","grey50"))+
  guides(fill="none")+
  theme_void()+
  theme(plot.margin = unit(c(-2.3,-2.3,-2.3,-2.3),units="lines"))

#(3) 柱状图，与取样地图组合。全球和每个生物地理区被外来鱼类和转移鱼类入侵的流域比例.
#全球
global.basin.number = data.all1%>%pull(X1.Basin.Name)%>%unique%>%length()#总流域数3119
global.basin.nonnative.number = data.all1%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length#有外来物种的流域数1719  55.11%
global.basin.alien.number = data.used.final%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#有国外引入外来物种的流域数1581  50.69%
global.basin.translocation.number = data.used.final.translocation%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#有流域间转移外来物种的流域数603  19.33%

#古北区
Palearctic.basin.number = data.all1%>%filter(X3.Ecoregion=="Palearctic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区流域数1033
Palearctic.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Palearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区有外来物种的流域数649  62.83%
Palearctic.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Palearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区有国外引入外来物种的流域数598  57.89%
Palearctic.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Palearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#古北区有流域间转移外来物种的流域数219  21.20%

#非洲热带区
Afrotropic.basin.number = data.all1%>%filter(X3.Ecoregion=="Afrotropic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区流域数282
Afrotropic.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Afrotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区有外来物种的流域数100  35.46%
Afrotropic.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Afrotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区有国外引入外来物种的流域数92  32.62%
Afrotropic.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Afrotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#非洲热带区有流域间转移外来物种的流域数35  12.41%

#新热带区
Neotropic.basin.number = data.all1%>%filter(X3.Ecoregion=="Neotropic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区流域数426
Neotropic.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Neotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区有外来物种的流域数169  39.67%
Neotropic.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Neotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区有国外引入外来物种的流域数152  35.68%
Neotropic.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Neotropic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新热带区有流域间转移外来物种的流域数57  13.38%

#澳新区
Australasia.basin.number = data.all1%>%filter(X3.Ecoregion=="Australasia")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区流域数725
Australasia.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Australasia"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区有外来物种的流域数397  54.76%
Australasia.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Australasia"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区有国外引入外来物种的流域数388  53.52%
Australasia.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Australasia"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#澳新区有流域间转移外来物种的流域数93  12.83%

#东洋区
Indo_Malay.basin.number = data.all1%>%filter(X3.Ecoregion=="Indo-Malay")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区流域数340
Indo_Malay.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Indo-Malay"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区有外来物种的流域数234  68.82%
Indo_Malay.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Indo-Malay"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区有国外引入外来物种的流域数220  64.71%
Indo_Malay.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Indo-Malay"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#东洋区有流域间转移外来物种的流域数48  14.12%

#新北区
Nearctic.basin.number = data.all1%>%filter(X3.Ecoregion=="Nearctic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区流域数296
Nearctic.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Nearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区有外来物种的流域数164  55.41%
Nearctic.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Nearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区有国外引入外来物种的流域数125  42.23%
Nearctic.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Nearctic"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#新北区有流域间转移外来物种的流域数151  51.01%

#大洋区
Oceania.basin.number = data.all1%>%filter(X3.Ecoregion=="Oceania")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区流域数17
Oceania.basin.nonnative.number = data.all1%>%filter(X3.Ecoregion=="Oceania"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区有外来物种的流域数6  35.29%
Oceania.basin.alien.number = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Oceania"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区有国外引入外来物种的流域数6  35.29%
Oceania.basin.translocation.number = data.used.final.translocation%>%left_join(drainage_basins)%>%filter(X3.Ecoregion=="Oceania"&X3.Native.Exotic.Status=="exotic")%>%pull(X1.Basin.Name)%>%unique%>%length()#大洋区有流域间转移外来物种的流域数0  0%

data.number = rbind(Global=c(global.basin.number,global.basin.nonnative.number,global.basin.alien.number,global.basin.translocation.number),
                   Palearctic=c(Palearctic.basin.number,Palearctic.basin.nonnative.number,Palearctic.basin.alien.number,Palearctic.basin.translocation.number),
                   Afrotropic=c(Afrotropic.basin.number,Afrotropic.basin.nonnative.number,Afrotropic.basin.alien.number,Afrotropic.basin.translocation.number),
                   Neotropic=c(Neotropic.basin.number,Neotropic.basin.nonnative.number,Neotropic.basin.alien.number,Neotropic.basin.translocation.number),
                   Australasia=c(Australasia.basin.number,Australasia.basin.nonnative.number,Australasia.basin.alien.number,Australasia.basin.translocation.number),
                   Indo_Malay=c(Indo_Malay.basin.number,Indo_Malay.basin.nonnative.number,Indo_Malay.basin.alien.number,Indo_Malay.basin.translocation.number),
                   Nearctic=c(Nearctic.basin.number,Nearctic.basin.nonnative.number,Nearctic.basin.alien.number,Nearctic.basin.translocation.number),
                   Oceania=c(Oceania.basin.number,Oceania.basin.nonnative.number,Oceania.basin.alien.number,Oceania.basin.translocation.number)
)

colnames(data.number) = c("total.number","nonnative.number","alien.number","translocation.number")         
data.number = data.number%>%as.data.frame%>%rownames_to_column(var="region")
data.number = data.number%>%mutate(nonnative.ratio=nonnative.number/total.number,alien.ratio=alien.number/total.number,translocation.ratio=translocation.number/total.number)

data.ratio = data.number%>%
  select(c(1,7:8))%>%
  pivot_longer(-1,names_to = "origin",values_to = "percentage")%>%
  filter(region!="Global")

data.ratio = data.ratio%>%mutate(origin=fct_relevel(origin,"translocation.ratio"))
data.ratio[14,3] = 0.003#仅仅为了画图显示出0值

ggplot(data=data.ratio,aes(y=region,x=percentage*100,fill=origin))+#大小3.5*4.5 portrait
  geom_col(position=position_dodge(),width=0.8)+
  xlab("Percentage of river basins (%)")+
  xlim(0,70)+
  #scale_fill_lancet(name=NULL,labels=c("colonization by alien fishes","colonization by translocated fishes"))+
  scale_fill_manual(name=NULL,values=c("#f47920","#00468BFF"),labels=c("Colonization by translocated fishes","Colonization by exotic fishes"))+
  scale_y_discrete(name=NULL,limits=data.ratio%>%filter(origin=="alien.ratio")%>%arrange(percentage)%>%pull(region))+
  theme_test()+
  theme(axis.text.y=element_text(size=10))+
  theme(legend.position=c(0.65,0.05),legend.background=element_blank(),legend.text = element_text(size=8))+
  guides(fill=guide_legend(keywidth=0.7,keyheight=0.7,reverse=T))

#17.Draw a phylogenetic tree

#以是否被引入过分组所有物种
species.out = phylo_all_spp$Insertions_data%>%filter(insertions=="Not_inserted")%>%pull(s)
phylo.status = data.occurrence%>%filter(!valid_names %in% species.out)%>%arrange(X3.Native.Exotic.Status)%>%distinct(valid_names,.keep_all=T)
grp = list(Native=phylo.status%>%filter(X3.Native.Exotic.Status=="native")%>%pull(valid_names),Exotic=phylo.status%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(valid_names))
group.phylo = groupOTU(phylo,grp,"Status")

ggtree = 
  ggtree(group.phylo,aes(color=Status),layout="circular",size=0.3)+
  scale_color_manual(name="Status of introduction",values=c("grey70","#00ae9d"),limits=c("Native","Exotic"),labels=c("Always native","Non-native species"))+
  theme(legend.position=c(0.45,0.57),legend.background=element_blank())+
  guides(color=guide_legend(keyheight=1,ncol=1,order=1,override.aes=list(linewidth=0.7)))

#---加物种对应的主要科热图，ggtree::gheatmap()
heatmap.family = phylo_all_spp$Insertions_data%>%filter(insertions!="Not_inserted")
main.family = heatmap.family%>%group_by(f)%>%summarise(n=length(s))%>%arrange(desc(n))%>%slice(1:8)%>%pull(f)
heatmap.main.family = heatmap.family%>%mutate(f=if_else(f %in% main.family,f,"NA"))
heatmap.main.family = heatmap.main.family[,1:2]%>%remove_rownames()%>%column_to_rownames(var="s")

ggtree.heatmap = #大小6.28*6.28
  gheatmap(ggtree,heatmap.main.family,offset=0,width=0.04,colnames=F,color = NA)+
  scale_fill_manual(name="Main family",limits=c(main.family,"NA"),values=c(pal_lancet()(8)[1:8],"white"),labels=c(main.family,""))+
  theme(legend.position=c(0.57,0.5),legend.background=element_blank())+
  guides(fill=guide_legend(keyheight=0.7,ncol=3,nrow=4,order=2))
ggtree.heatmap

#---采用ggtreeExtra包,加外圈热图和柱状图，使用😄！
#第一圈，热图，物种所属的科
heatmap.family = phylo_all_spp$Insertions_data%>%filter(insertions!="Not_inserted")
main.family1 = heatmap.family%>%group_by(f)%>%summarise(n=length(s))%>%arrange(desc(n))%>%slice(1:9)%>%pull(f)#主要的9个科
ggtreeExtra.main.family = heatmap.family%>%mutate(f=if_else(f %in% main.family1,f,"Others"))
ggtreeExtra.main.family%>%tibble

p1 = ggtree+
    ggtreeExtra::geom_fruit(
    data=ggtreeExtra.main.family,
    geom=geom_tile,
    aes(x="",y=s,fill=f),#x=""使得圈在同一位置😄，否则就根据x分配位置了。y是tree中tiplabel即种名😄
    width=25,
    offset=0.025)+
    scale_fill_manual(name="Main families",limits=c(main.family1,"Others"),values=c(pal_lancet()(9)[c(1:7)],"#fcaf17","#6c4c49","grey80"))+#c(pal_simpsons()(16)[c(12:16,4:7,3)])
    theme(legend.position="right")+
    guides(fill=guide_legend(keyheight=0.5,keywidth=0.5,ncol=2,order=2))

#第二圈，热图，物种所属的目
heatmap.order = phylo_all_spp$Insertions_data%>%filter(insertions!="Not_inserted")
main.order = heatmap.order%>%group_by(o)%>%summarise(n=length(s))%>%arrange(desc(n))%>%slice(1:9)%>%pull(o)#主要的9个科
ggtreeExtra.main.order = heatmap.order%>%mutate(o=if_else(o %in% main.order,o,"Others"))
ggtreeExtra.main.order%>%tibble

p2 = p1+
  ggnewscale::new_scale_fill()+
  geom_fruit(
  data=ggtreeExtra.main.order,
  geom=geom_tile,
  aes(x="",y=s,fill=o),
  width=25,
  offset=0.04)+
  scale_fill_manual(name="Main orders",limits=c(main.order,"Others"),values=c(pal_frontiers()(9)[c(1,8:9,2:7)],"grey80"))+
  guides(fill=guide_legend(keyheight=0.5,keywidth=0.5,ncol=2,order=3))

#给第二圈外层加线，通过压缩tile得到线，其他方式不好得到
line.ecoregion = data.frame(x=1,valid_names=phylo.ecoregion$valid_names%>%unique())
p3 = p2+
  geom_fruit(
    data=line.ecoregion,
    geom=geom_tile,
    aes(x=x,y=valid_names),
    fill=NA,color="black",
    pwidth=0.01,
    offset=0.04)

#第三圈，热图，物种所在的生物地理区
phylo.ecoregion = #创造每个ecoregion里对应的物种
  data.all%>%
  filter(!valid_names %in% species.out)%>%
  select(c(1:2,4),9:10)%>%
  filter(X3.Ecoregion!="Oceania")%>%
  group_by(X3.Ecoregion)%>%
  reframe(valid_names=unique(valid_names))#reframe相当于summarise,只是总结结果为多行时新版dplyr用这个reframe代替了summarise

p4 = p3+
  new_scale_fill()+
  geom_fruit(
    data=phylo.ecoregion,
    geom=geom_tile,
    aes(x=X3.Ecoregion,y=valid_names,fill=X3.Ecoregion),
    color=NA,
    pwidth=0.3,
    offset=0.018)+
  scale_fill_npg(name="Biogeographical realms")+
  guides(fill=guide_legend(keyheight=0.5,keywidth=0.5,ncol=2,order=4))

#给第三圈外层加线，通过压缩tile得到线
p5 = p4+#大小8.15*6.15
  geom_fruit(
    data=line.ecoregion,
    geom=geom_tile,
    aes(x=x,y=valid_names),
    fill=NA,color="black",
    pwidth=0.01,
    offset=0.03)

p5#最终使用

#18.The relationship between the occurrence likelihood and phylogenetic distance in different biogeographic regions. Alien species.

data.distance.final_ecoregion = data.distance.final%>%left_join(drainage_basins,by=c("country"="X2.Country","basin"="X1.Basin.Name"))%>%select(1:7)

#古北区
glmer.palearctic_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic"))
summary(glmer.palearctic_mpd)
#(plot.palearctic.mpd = ggplot(get_model_data(glmer.palearctic_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="Probability of exotic occurrence")+annotate("text",x=c(0.5,0.5),y=c(0.75,0.65),label=c("italic(P)<2e-16","italic(n)==14196"),parse=T))#画出的不平滑，因为x不是均匀分布的，下面自己手动画
newdata.palearctic.mpd = data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic")%>%pull(mpd)
predict.data.palearctic.mpd = predict(glmer.palearctic_mpd,type="response",re.form=NA,newdata=data.frame(mpd=seq(min(newdata.palearctic.mpd),max(newdata.palearctic.mpd),length=500)))#注意mpd=必须要
predict.data.palearctic.mpd = data.frame(mpd=seq(min(newdata.palearctic.mpd),max(newdata.palearctic.mpd),length=500),predicted=predict.data.palearctic.mpd)
(plot.palearctic.mpd = ggplot(predict.data.palearctic.mpd,aes(mpd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[7])+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.5,0.5),y=c(0.75,0.65),label=c("italic(P)<2e-16","italic(n)==14196"),parse=T))
                             
glmer.palearctic_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic"))
summary(glmer.palearctic_mntd)
#(plot.palearctic.mntd = ggplot(get_model_data(glmer.palearctic_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.5,0.5),y=c(0.52,0.45),label=c("italic(P)<2e-16","italic(n)==14196"),parse=T))#画出的不平滑，因为x不是均匀分布的，下面自己手动画
newdata.palearctic.mntd = data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic")%>%pull(mntd)
predict.data.palearctic.mntd = predict(glmer.palearctic_mntd,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.palearctic.mntd),max(newdata.palearctic.mntd),length=500)))#注意mntd=必须要
predict.data.palearctic.mntd = data.frame(mntd=seq(min(newdata.palearctic.mntd),max(newdata.palearctic.mntd),length=500),predicted=predict.data.palearctic.mntd)
(plot.palearctic.mntd = ggplot(predict.data.palearctic.mntd,aes(mntd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[7])+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.5,0.5),y=c(0.52,0.45),label=c("italic(P)<2e-16","italic(n)==14196"),parse=T))

#非洲热带区
glmer.afrotropic_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Afrotropic"))
summary(glmer.afrotropic_mpd)
(plot.afrotropic.mpd = ggplot(get_model_data(glmer.afrotropic_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2,color=pal_npg()(7)[1])+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.255,0.25),y=c(0.13,0.123),label=c("italic(P)==0.406","italic(n)==1722"),parse=T))

glmer.afrotropic_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Afrotropic"))
summary(glmer.afrotropic_mntd)
(plot.afrotropic.mntd = ggplot(get_model_data(glmer.afrotropic_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2,color=pal_npg()(7)[1])+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.255,0.25),y=c(0.126,0.117),label=c("italic(P)==0.089","italic(n)==1722"),parse=T))

#新热带区
glmer.neotropic_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Neotropic"))
summary(glmer.neotropic_mpd)
#(plot.neotropic.mpd = ggplot(get_model_data(glmer.neotropic_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="Probability of exotic occurrence")+annotate("text",x=c(0.51,0.5),y=c(0.05,0.044),label=c("italic(P)==0.009","italic(n)==6277"),parse=T))
newdata.neotropic.mpd = data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Neotropic")%>%pull(mpd)
predict.data.neotropic.mpd = predict(glmer.neotropic_mpd,type="response",re.form=NA,newdata=data.frame(mpd=seq(min(newdata.neotropic.mpd),max(newdata.neotropic.mpd),length=500)))#注意mpd=必须要
predict.data.neotropic.mpd = data.frame(mpd=seq(min(newdata.neotropic.mpd),max(newdata.neotropic.mpd),length=500),predicted=predict.data.neotropic.mpd)
(plot.neotropic.mpd = ggplot(predict.data.neotropic.mpd,aes(mpd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[5])+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.51,0.5),y=c(0.05,0.044),label=c("italic(P)==0.009","italic(n)==6277"),parse=T))

glmer.neotropic_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Neotropic"))
summary(glmer.neotropic_mntd)
#(plot.neotropic.mntd = ggplot(get_model_data(glmer.neotropic_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.56,0.5),y=c(0.084,0.074),label=c("italic(P)==3.37e-7","italic(n)==6277"),parse=T))
newdata.neotropic.mntd = data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Neotropic")%>%pull(mntd)
predict.data.neotropic.mntd = predict(glmer.neotropic_mntd,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.neotropic.mntd),max(newdata.neotropic.mntd),length=500)))#注意mntd=必须要
predict.data.neotropic.mntd = data.frame(mntd=seq(min(newdata.neotropic.mntd),max(newdata.neotropic.mntd),length=500),predicted=predict.data.neotropic.mntd)
(plot.neotropic.mntd = ggplot(predict.data.neotropic.mntd,aes(mntd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[5])+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.56,0.5),y=c(0.084,0.074),label=c("italic(P)==3.37e-7","italic(n)==6277"),parse=T))

#澳新区
glmer.australasia_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Australasia"))
summary(glmer.australasia_mpd)
(plot.australasia.mpd = ggplot(get_model_data(glmer.australasia_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[2])+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.41,0.4),y=c(0.051,0.044),label=c("italic(P)==4.57e-7","italic(n)==15391"),parse=T))

glmer.australasia_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Australasia"))
summary(glmer.australasia_mntd)
(plot.australasia.mntd = ggplot(get_model_data(glmer.australasia_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[2])+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.301,0.3),y=c(0.06,0.052),label=c("italic(P)==0.0001","italic(n)==15391"),parse=T))

#东洋区
glmer.indomalay_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Indo-Malay"))
summary(glmer.indomalay_mpd)
(plot.indomalay.mpd = ggplot(get_model_data(glmer.indomalay_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[3])+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.303,0.3),y=c(0.158,0.145),label=c("italic(P)==0.018","italic(n)==4388"),parse=T))

glmer.indomalay_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Indo-Malay"))
summary(glmer.indomalay_mntd)
(plot.indomalay.mntd = ggplot(get_model_data(glmer.indomalay_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[3])+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.305,0.3),y=c(0.18,0.16),label=c("italic(P)<2e-16","italic(n)==4388"),parse=T))

#新北区
glmer.nearctic_mpd = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Nearctic"))
summary(glmer.nearctic_mpd)
(plot.nearctic.mpd = ggplot(get_model_data(glmer.nearctic_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[4])+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.403,0.4),y=c(0.04,0.035),label=c("italic(P)<2e-16","italic(n)==19065"),parse=T))

glmer.nearctic_mntd = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Nearctic"))
summary(glmer.nearctic_mntd)
(plot.nearctic.mntd = ggplot(get_model_data(glmer.nearctic_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[4])+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.4,0.4),y=c(0.006,0.0053),label=c("italic(P)<2e-16","italic(n)==19065"),parse=T))

#大洋区，样本少，自变量range也非常小，基本无法拟合。不单独呈现大洋区
glmer.oceania_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Oceania"))#只有一个国家波利尼西亚，country不能作为随机效应
summary(glmer.oceania_mpd)
plot_model(glmer.oceania_mpd,type="pred",terms="mpd [all]")

glmer.oceania_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final_ecoregion%>%filter(X3.Ecoregion=="Oceania"))#只有一个国家波利尼西亚，country不能作为随机效应,拟合不出来，最近谱系距离非常接近，样本量也少
summary(glmer.oceania_mntd)
plot_model(glmer.oceania_mntd,type="pred",terms="mntd [all]")

#19.The relationship between the occurrence likelihood and phylogenetic distance in different biogeographic regions. Translocated species.

data.distance.translocation.final_ecoregion = data.distance.translocation.final%>%left_join(drainage_basins,by=c("country"="X2.Country","basin"="X1.Basin.Name"))%>%select(1:7)

#古北区
glmer.palearctic_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic"))
summary(glmer.palearctic_mpd.translocation)
#(plot.palearctic.mpd_translocation = ggplot(get_model_data(glmer.palearctic_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="MPD",y="Probability of exotic occurrence")+annotate("text",x=c(0.5,0.5),y=c(0.17,0.15),label=c("italic(P)==0.020","italic(n)==8080"),parse=T))#不平滑，因x分布不均匀，下面自己手动画
newdata.palearctic.mpd.translocation = data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic")%>%pull(mpd)
predict.data.palearctic.mpd.translocation = predict(glmer.palearctic_mpd.translocation,type="response",re.form=NA,newdata=data.frame(mpd=seq(min(newdata.palearctic.mpd.translocation),max(newdata.palearctic.mpd.translocation),length=500)))#注意mpd=必须要
predict.data.palearctic.mpd.translocation = data.frame(mpd=seq(min(newdata.palearctic.mpd.translocation),max(newdata.palearctic.mpd.translocation),length=500),predicted=predict.data.palearctic.mpd.translocation)
(plot.palearctic.mpd_translocation = ggplot(predict.data.palearctic.mpd.translocation,aes(mpd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[7])+theme_classic()+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.5,0.5),y=c(0.17,0.15),label=c("italic(P)==0.020","italic(n)==8080"),parse=T))

glmer.palearctic_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic"))
summary(glmer.palearctic_mntd.translocation)
#(plot.palearctic.mntd_translocation = ggplot(get_model_data(glmer.palearctic_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="MNTD",y="")+annotate("text",x=c(0.56,0.5),y=c(0.15,0.13),label=c("italic(P)==3.82e-11","italic(n)==8080"),parse=T))#不平滑，因x分布不均匀，下面自己手动画
newdata.palearctic.mntd.translocation = data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Palearctic")%>%pull(mntd)
predict.data.palearctic.mntd.translocation = predict(glmer.palearctic_mntd.translocation,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.palearctic.mntd.translocation),max(newdata.palearctic.mntd.translocation),length=500)))#注意mntd=必须要
predict.data.palearctic.mntd.translocation = data.frame(mntd=seq(min(newdata.palearctic.mntd.translocation),max(newdata.palearctic.mntd.translocation),length=500),predicted=predict.data.palearctic.mntd.translocation)
(plot.palearctic.mntd_translocation = ggplot(predict.data.palearctic.mntd.translocation,aes(mntd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[7])+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.56,0.5),y=c(0.15,0.13),label=c("italic(P)==3.82e-11","italic(n)==8080"),parse=T))

#非洲热带区
glmer.afrotropic_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Afrotropic"))
summary(glmer.afrotropic_mpd.translocation)
(plot.afrotropic.mpd_translocation = ggplot(get_model_data(glmer.afrotropic_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2,color=pal_npg()(7)[1])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.255,0.25),y=c(0.00016,0.00014),label=c("italic(P)==0.470","italic(n)==356"),parse=T))

glmer.afrotropic_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Afrotropic"))
summary(glmer.afrotropic_mntd.translocation)
(plot.afrotropic.mntd_translocation = ggplot(get_model_data(glmer.afrotropic_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[1])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.265,0.25),y=c(0.00071,0.00061),label=c("italic(P)==0.002","italic(n)==356"),parse=T))

#新热带区
glmer.neotropic_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Neotropic"))
summary(glmer.neotropic_mpd.translocation)
(plot.neotropic.mpd_translocation = ggplot(get_model_data(glmer.neotropic_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2,color=pal_npg()(7)[5])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.354,0.35),y=c(3e-4,2.7e-4),label=c("italic(P)==0.243","italic(n)==4134"),parse=T))

glmer.neotropic_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Neotropic"))
summary(glmer.neotropic_mntd.translocation)
(plot.neotropic.mntd_translocation = ggplot(get_model_data(glmer.neotropic_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[5])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.285,0.25),y=c(0.0019,0.0016),label=c("italic(P)==4.08e-6","italic(n)==4134"),parse=T))

#澳新区
glmer.australasia_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Australasia"))
summary(glmer.australasia_mpd.translocation)
(plot.australasia.mpd_translocation = ggplot(get_model_data(glmer.australasia_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2,color=pal_npg()(7)[2])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.305,0.3),y=c(5.3e-8,4.7e-8),label=c("italic(P)==0.516","italic(n)==7998"),parse=T))

glmer.australasia_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Australasia"))
summary(glmer.australasia_mntd.translocation)
(plot.australasia.mntd_translocation = ggplot(get_model_data(glmer.australasia_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[2])+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.256,0.25),y=c(0.007,0.006),label=c("italic(P)==0.0002","italic(n)==7998"),parse=T))

#东洋区
glmer.indomalay_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Indo-Malay"))
summary(glmer.indomalay_mpd.translocation)
(plot.indomalay.mpd_translocation = ggplot(get_model_data(glmer.indomalay_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2,color=pal_npg()(7)[3])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.305,0.3),y=c(2e-4,1.85e-4),label=c("italic(P)==0.584","italic(n)==1869"),parse=T))

glmer.indomalay_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Indo-Malay"))
summary(glmer.indomalay_mntd.translocation)
(plot.indomalay.mntd_translocation = ggplot(get_model_data(glmer.indomalay_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[3])+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.256,0.25),y=c(0.00062,0.00055),label=c("italic(P)==0.003","italic(n)==1869"),parse=T))

#新北区
glmer.nearctic_mpd.translocation = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Nearctic"))
summary(glmer.nearctic_mpd.translocation)
#(plot.nearctic.mpd_translocation = ggplot(get_model_data(glmer.nearctic_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+scale_y_continuous(label=scientific)+labs(x="MPD",y="Probability of exotic occurrence")+annotate("text",x=c(0.445,0.4),y=c(4e-3,3.4e-3),label=c("italic(P)==5.09e-12","italic(n)==40514"),parse=T))
newdata.nearctic.mpd.translocation = data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Nearctic")%>%pull(mpd)
predict.data.nearctic.mpd.translocation = predict(glmer.nearctic_mpd.translocation,type="response",re.form=NA,newdata=data.frame(mpd=seq(min(newdata.nearctic.mpd.translocation),max(newdata.nearctic.mpd.translocation),length=500)))#注意mpd=必须要
predict.data.nearctic.mpd.translocation = data.frame(mpd=seq(min(newdata.nearctic.mpd.translocation),max(newdata.nearctic.mpd.translocation),length=500),predicted=predict.data.nearctic.mpd.translocation)
(plot.nearctic.mpd_translocation = ggplot(predict.data.nearctic.mpd.translocation,aes(mpd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[4])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.445,0.4),y=c(4e-3,3.4e-3),label=c("italic(P)==5.09e-12","italic(n)==40514"),parse=T))

glmer.nearctic_mntd.translocation = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Nearctic"))
summary(glmer.nearctic_mntd.translocation)
#(plot.nearctic.mntd_translocation = ggplot(get_model_data(glmer.nearctic_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+scale_y_continuous(label=scientific)+labs(x="MNTD",y="")+annotate("text",x=c(0.405,0.4),y=c(0.0013,0.00113),label=c("italic(P)<2e-16","italic(n)==40514"),parse=T))
newdata.nearctic.mntd.translocation = data.distance.translocation.final_ecoregion%>%filter(X3.Ecoregion=="Nearctic")%>%pull(mntd)
predict.data.nearctic.mntd.translocation = predict(glmer.nearctic_mntd.translocation,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.nearctic.mntd.translocation),max(newdata.nearctic.mntd.translocation),length=500)))#注意mntd=必须要
predict.data.nearctic.mntd.translocation = data.frame(mntd=seq(min(newdata.nearctic.mntd.translocation),max(newdata.nearctic.mntd.translocation),length=500),predicted=predict.data.nearctic.mntd.translocation)
(plot.nearctic.mntd_translocation = ggplot(predict.data.nearctic.mntd.translocation,aes(mntd,predicted))+geom_line(linewidth=2,color=pal_npg()(7)[4])+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.405,0.4),y=c(0.0013,0.00113),label=c("italic(P)<2e-16","italic(n)==40514"),parse=T))

#大洋区,就一个国家，并且不存在流域间转移的外来种，没有拟合模型

#合并各个生物地理区两个尺度的四个logistic图,大小6.7*6.28
ggarrange(#古北区
  plot.palearctic.mpd,plot.palearctic.mntd,
  plot.palearctic.mpd_translocation,plot.palearctic.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#非洲热带区
  plot.afrotropic.mpd,plot.afrotropic.mntd,
  plot.afrotropic.mpd_translocation,plot.afrotropic.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#新热带区
  plot.neotropic.mpd,plot.neotropic.mntd,
  plot.neotropic.mpd_translocation,plot.neotropic.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#澳新区
  plot.australasia.mpd,plot.australasia.mntd,
  plot.australasia.mpd_translocation,plot.australasia.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#东洋区
  plot.indomalay.mpd,plot.indomalay.mntd,
  plot.indomalay.mpd_translocation,plot.indomalay.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#新北区
  plot.nearctic.mpd,plot.nearctic.mntd,
  plot.nearctic.mpd_translocation,plot.nearctic.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

#20.Chord diagram depicting the relationships between alien species and countries in each biogeographic region
data.used.final_chord = data.used.final%>%left_join(drainage_basins)%>%filter(X3.Native.Exotic.Status=="exotic")%>%distinct(X2.Country,valid_names,.keep_all = T)

#古北区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Palearctic")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Palearctic")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Oncorhynchus_mykiss被引入到最多国家，引入到 32 countries
data.used.final_chord%>%filter(X3.Ecoregion=="Palearctic")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#Spain是引入外来种最多的国家，引入了 27 species

#非洲热带区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Afrotropic")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Afrotropic")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Micropterus_salmoides被引入到最多国家，引入到 8 countries
data.used.final_chord%>%filter(X3.Ecoregion=="Afrotropic")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#South Africa是引入外来种最多的国家，引入了 22 species

#新热带区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Neotropic")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Neotropic")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Oreochromis_mossambicus被引入到最多国家，引入到 9 countries
data.used.final_chord%>%filter(X3.Ecoregion=="Neotropic")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#Brazil是引入外来种最多的国家，引入了 32 species

#澳新区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Australasia")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Australasia")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Carassius_auratus被引入到最多国家，引入到 4c ountries
data.used.final_chord%>%filter(X3.Ecoregion=="Australasia")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#Australia是引入外来种最多的国家，引入了 24 species

#东洋区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Indo-Malay")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Indo-Malay")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Oreochromis_mossambicus被引入到最多国家，引入到 13 countries
data.used.final_chord%>%filter(X3.Ecoregion=="Indo-Malay")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#Nepal是引入外来种最多的国家，引入了 24 species

#新北区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Nearctic")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Nearctic")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Salmo_trutta被引入到最多国家，引入到 3 countries
data.used.final_chord%>%filter(X3.Ecoregion=="Nearctic")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#United States是引入外来种最多的国家，引入了 84 species

#大洋区 大小5*5
chordDiagram(data.used.final_chord%>%filter(X3.Ecoregion=="Oceania")%>%select(c(4,1)),annotationTrack = c("grid","axis"))
data.used.final_chord%>%filter(X3.Ecoregion=="Oceania")%>%select(c(4,1))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Oreochromis_mossambicus被引入到最多国家，引入到 1 countries
data.used.final_chord%>%filter(X3.Ecoregion=="Oceania")%>%select(c(4,1))%>%group_by(X2.Country)%>%summarise(n=n())%>%arrange(desc(n))#French Polynesia是引入外来种最多的国家，引入了 3 species


#21.The relationship between the occurrence likelihood and phylogenetic distance in different countries. Alien species.

#外来物种最多的国家，国外引入+流域间转移
data.all1%>%filter(X3.Native.Exotic.Status=="exotic")%>%group_by(X2.Country)%>%distinct(valid_names,.keep_all = T)%>%summarise(n=length(valid_names))%>%arrange(desc(n))

#按照所有外来物种最多的6个国家，选择的国家如下：United States,Canada,Brazil,Russia,Mexico,China

#查看国外引入外来物种最多的国家，仅是查看，选择展示的国家仍按所有外来物种的多少的前6名
data.used.final%>%filter(X3.Native.Exotic.Status=="exotic")%>%group_by(X2.Country)%>%distinct(valid_names,.keep_all = T)%>%summarise(n=length(valid_names))%>%arrange(desc(n))

#United States
glmer.US_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final%>%filter(country=="United States"))
summary(glmer.US_mpd)
(plot.US.mpd = ggplot(get_model_data(glmer.US_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.325,0.3),y=c(0.25,0.22),label=c("italic(P)==5.05E-13","italic(n)==16296"),parse=T))

glmer.US_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.final%>%filter(country=="United States"))
summary(glmer.US_mntd)
(plot.US.mntd = ggplot(get_model_data(glmer.US_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.252,0.25),y=c(0.016,0.014),label=c("italic(P)<2e-16","italic(n)==16296"),parse=T))

#Canada
glmer.Canada_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Canada"))
summary(glmer.Canada_mpd)
(plot.Canada.mpd = ggplot(get_model_data(glmer.Canada_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.51,0.5),y=c(3e-4,2.7e-4),label=c("italic(P)==0.055","italic(n)==2160"),parse=T))

glmer.Canada_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Canada"))
summary(glmer.Canada_mntd)
#(plot.Canada.mntd = ggplot(get_model_data(glmer.Canada_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=1)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.51,0.5),y=c(6e-5,5.4e-5),label=c("italic(P)==0.009","italic(n)==2160"),parse=T))
newdata.Canada.mntd = data.distance.final%>%filter(country=="Canada")%>%pull(mntd)
predict.data.Canada.mntd = predict(glmer.Canada_mntd,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.Canada.mntd),max(newdata.Canada.mntd),length=500)))#注意mntd=必须要
predict.data.Canada.mntd = data.frame(mntd=seq(min(newdata.Canada.mntd),max(newdata.Canada.mntd),length=500),predicted=predict.data.Canada.mntd)
(plot.Canada.mntd = ggplot(predict.data.Canada.mntd,aes(mntd,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.51,0.5),y=c(6e-5,5.4e-5),label=c("italic(P)==0.009","italic(n)==2160"),parse=T))

#Brazil
glmer.Brazil_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Brazil"))
summary(glmer.Brazil_mpd)
#(plot.Brazil.mpd = ggplot(get_model_data(glmer.Brazil_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=1,linetype=2)+theme_classic()+labs(x=NULL,y="Probability of exotic occurrence")+annotate("text",x=c(0.605,0.6),y=c(0.005,0.0045),label=c("italic(P)==0.316","italic(n)==2944"),parse=T))
newdata.Brazil.mpd = data.distance.final%>%filter(country=="Brazil")%>%pull(mpd)
predict.data.Brazil.mpd = predict(glmer.Brazil_mpd,type="response",re.form=NA,newdata=data.frame(mpd=seq(min(newdata.Brazil.mpd),max(newdata.Brazil.mpd),length=500)))#注意mpd=必须要
predict.data.Brazil.mpd = data.frame(mpd=seq(min(newdata.Brazil.mpd),max(newdata.Brazil.mpd),length=500),predicted=predict.data.Brazil.mpd)
(plot.Brazil.mpd = ggplot(predict.data.Brazil.mpd,aes(mpd,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.605,0.6),y=c(0.005,0.0045),label=c("italic(P)==0.316","italic(n)==2944"),parse=T))

glmer.Brazil_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Brazil"))
summary(glmer.Brazil_mntd)
#(plot.Brazil.mntd = ggplot(get_model_data(glmer.Brazil_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=1)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.505,0.5),y=c(0.0073,0.0065),label=c("italic(P)==0.047","italic(n)==2944"),parse=T))
newdata.Brazil.mntd = data.distance.final%>%filter(country=="Brazil")%>%pull(mntd)
predict.data.Brazil.mntd = predict(glmer.Brazil_mntd,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.Brazil.mntd),max(newdata.Brazil.mntd),length=500)))#注意mntd=必须要
predict.data.Brazil.mntd = data.frame(mntd=seq(min(newdata.Brazil.mntd),max(newdata.Brazil.mntd),length=500),predicted=predict.data.Brazil.mntd)
(plot.Brazil.mntd = ggplot(predict.data.Brazil.mntd,aes(mntd,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.505,0.5),y=c(0.0073,0.0065),label=c("italic(P)==0.047","italic(n)==2944"),parse=T))

#Russia
glmer.Russia_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Russia"))
summary(glmer.Russia_mpd)
(plot.Russia.mpd = ggplot(get_model_data(glmer.Russia_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.505,0.5),y=c(2.6e-5,2.35e-5),label=c("italic(P)==0.155","italic(n)==2678"),parse=T))

glmer.Russia_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Russia"))
summary(glmer.Russia_mntd)
(plot.Russia.mntd = ggplot(get_model_data(glmer.Russia_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.505,0.5),y=c(1.45e-5,1.35e-5),label=c("italic(P)==0.401","italic(n)==2678"),parse=T))

#Mexico
glmer.Mexico_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Mexico"))
summary(glmer.Mexico_mpd)
(plot.Mexico.mpd = ggplot(get_model_data(glmer.Mexico_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.355,0.35),y=c(0.015,0.0135),label=c("italic(P)==0.003","italic(n)==2407"),parse=T))

glmer.Mexico_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="Mexico"))
summary(glmer.Mexico_mntd)
(plot.Mexico.mntd = ggplot(get_model_data(glmer.Mexico_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.255,0.25),y=c(0.0255,0.023),label=c("italic(P)==2e-6","italic(n)==2407"),parse=T))

#China
glmer.China_mpd = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="China"))
summary(glmer.China_mpd)
(plot.China.mpd = ggplot(get_model_data(glmer.China_mpd,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x=NULL,y="Occurrence probability")+annotate("text",x=c(0.265,0.25),y=c(0.0225,0.02),label=c("italic(P)==0.189","italic(n)==864"),parse=T))

glmer.China_mntd = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final%>%filter(country=="China"))
summary(glmer.China_mntd)
(plot.China.mntd = ggplot(get_model_data(glmer.China_mntd,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x=NULL,y="")+annotate("text",x=c(0.265,0.25),y=c(0.025,0.0225),label=c("italic(P)==0.169","italic(n)==864"),parse=T))

#22.The relationship between the occurrence likelihood and phylogenetic distance in different countries. Translocated species.

#按照以上标准选择的国家，United States,Canada,Brazil,Russia,Mexico,China

#查看流域间转移外来物种最多的国家，仅是查看，选择展示的国家仍按所有外来物种的多少的前6名
data.used.final.translocation%>%filter(X3.Native.Exotic.Status=="exotic")%>%group_by(X2.Country)%>%distinct(valid_names,.keep_all = T)%>%summarise(n=length(valid_names))%>%arrange(desc(n))
#注意因为前期data.used数据就排除了不在谱系树中的物种，data.exclude包括了在4个国家被定义为exotic的6个物种(这6个种都是流域间转移种），展示流域间外来物种数时需要加上这些物种

#United States
glmer.US_mpd.translocation = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="United States"))
summary(glmer.US_mpd.translocation)
(plot.US.mpd_translocation = ggplot(get_model_data(glmer.US_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.54,0.5),y=c(0.0075,0.0066),label=c("italic(P)==1.54e-7","italic(n)==37491"),parse=T))

glmer.US_mntd.translocation = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="United States"))
summary(glmer.US_mntd.translocation)
(plot.US.mntd_translocation = ggplot(get_model_data(glmer.US_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.52,0.5),y=c(0.0028,0.0025),label=c("italic(P)==3.3e-16","italic(n)==37491"),parse=T))

#Canada
glmer.Canada_mpd.translocation = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final%>%filter(country=="Canada"))
summary(glmer.Canada_mpd.translocation)
(plot.Canada.mpd_translocation = ggplot(get_model_data(glmer.Canada_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.31,0.3),y=c(0.0002340,0.0002330),label=c("italic(P)==0.981","italic(n)==2468"),parse=T))

glmer.Canada_mntd.translocation = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final%>%filter(country=="Canada"))
summary(glmer.Canada_mntd.translocation)
(plot.Canada.mntd_translocation = ggplot(get_model_data(glmer.Canada_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.25,0.25),y=c(0.0002368,0.0002364),label=c("italic(P)==0.980","italic(n)==2468"),parse=T))

#Brazil
glmer.Brazil_mpd.translocation = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="Brazil"))
summary(glmer.Brazil_mpd.translocation)
(plot.Brazil.mpd_translocation = ggplot(get_model_data(glmer.Brazil_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.35,0.35),y=c(0.00235,0.0022),label=c("italic(P)==0.316","italic(n)==2334"),parse=T))

glmer.Brazil_mntd.translocation = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="Brazil"))
summary(glmer.Brazil_mntd.translocation)
#(plot.Brazil.mntd_translocation = ggplot(get_model_data(glmer.Brazil_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=1)+theme_classic()+labs(x="MNTD",y="Probability of exotic occurrence")+annotate("text",x=c(0.25,0.25),y=c(0.006,0.0054),label=c("italic(P)==0.016","italic(n)==2334"),parse=T))
newdata.Brazil.mntd.translocation = data.distance.translocation.final%>%filter(country=="Brazil")%>%pull(mntd)
predict.data.Brazil.mntd.translocation = predict(glmer.Brazil_mntd.translocation,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.Brazil.mntd.translocation),max(newdata.Brazil.mntd.translocation),length=500)))#注意mntd=必须要
predict.data.Brazil.mntd.translocation = data.frame(mntd=seq(min(newdata.Brazil.mntd.translocation),max(newdata.Brazil.mntd.translocation),length=500),predicted=predict.data.Brazil.mntd.translocation)
(plot.Brazil.mntd_translocation = ggplot(predict.data.Brazil.mntd.translocation,aes(mntd,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.25,0.25),y=c(0.006,0.0054),label=c("italic(P)==0.016","italic(n)==2334"),parse=T))

#Russia
glmer.Russia_mpd.translocation = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="Russia"))
summary(glmer.Russia_mpd.translocation)
(plot.Russia.mpd_translocation = ggplot(get_model_data(glmer.Russia_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+scale_y_continuous(label=scientific)+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.355,0.35),y=c(0.00013,0.000115),label=c("italic(P)==0.085","italic(n)==2716"),parse=T))

glmer.Russia_mntd.translocation = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="Russia"))
summary(glmer.Russia_mntd.translocation)
#(plot.Russia.mntd_translocation = ggplot(get_model_data(glmer.Russia_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=1)+theme_classic()+labs(x="MNTD",y="")+annotate("text",x=c(0.5,0.5),y=c(0.9e-4,0.8e-4),label=c("italic(P)==0.048","italic(n)==2716"),parse=T))
newdata.Russia.mntd.translocation = data.distance.translocation.final%>%filter(country=="Russia")%>%pull(mntd)
predict.data.Russia.mntd.translocation = predict(glmer.Russia_mntd.translocation,type="response",re.form=NA,newdata=data.frame(mntd=seq(min(newdata.Russia.mntd.translocation),max(newdata.Russia.mntd.translocation),length=500)))#注意mntd=必须要
predict.data.Russia.mntd.translocation = data.frame(mntd=seq(min(newdata.Russia.mntd.translocation),max(newdata.Russia.mntd.translocation),length=500),predicted=predict.data.Russia.mntd.translocation)
(plot.Russia.mntd_translocation = ggplot(predict.data.Russia.mntd.translocation,aes(mntd,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.5,0.5),y=c(0.9e-4,0.8e-4),label=c("italic(P)==0.048","italic(n)==2716"),parse=T))

#Mexico
glmer.Mexico_mpd.translocation = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="Mexico"))
summary(glmer.Mexico_mpd.translocation)
(plot.Mexico.mpd_translocation = ggplot(get_model_data(glmer.Mexico_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.355,0.35),y=c(0.013,0.0118),label=c("italic(P)==0.096","italic(n)==2071"),parse=T))

glmer.Mexico_mntd.translocation = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final%>%filter(country=="Mexico"))
summary(glmer.Mexico_mntd.translocation)
(plot.Mexico.mntd_translocation = ggplot(get_model_data(glmer.Mexico_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.265,0.25),y=c(0.025,0.022),label=c("italic(P)==1.1e-6","italic(n)==2071"),parse=T))

#China
glmer.China_mpd.translocation = glmer(occurrence~mpd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final%>%filter(country=="China"))
summary(glmer.China_mpd.translocation)
(plot.China.mpd_translocation = ggplot(get_model_data(glmer.China_mpd.translocation,type="pred",terms="mpd [all]"),aes(x,predicted))+geom_line(linewidth=2,linetype=2)+theme_classic()+labs(x="Nonnative-native MPD",y="Occurrence probability")+annotate("text",x=c(0.257,0.25),y=c(0.009,0.00845),label=c("italic(P)==0.635","italic(n)==2367"),parse=T))

glmer.China_mntd.translocation = glmer(occurrence~mntd+(1|basin)+(1|exotic_species),family=binomial(link = "logit"),data=data.distance.translocation.final%>%filter(country=="China"))
summary(glmer.China_mntd.translocation)
(plot.China.mntd_translocation = ggplot(get_model_data(glmer.China_mntd.translocation,type="pred",terms="mntd [all]"),aes(x,predicted))+geom_line(linewidth=2)+theme_classic()+labs(x="Nonnative-native MNTD",y="")+annotate("text",x=c(0.255,0.25),y=c(0.015,0.013),label=c("italic(P)==0.002","italic(n)==2367"),parse=T))

#合并6个国家两个尺度的四个logistic图，大小6.7*6.28

ggarrange(#美国
  plot.US.mpd,plot.US.mntd,
  plot.US.mpd_translocation,plot.US.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#加拿大
  plot.Canada.mpd,plot.Canada.mntd,
  plot.Canada.mpd_translocation,plot.Canada.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#巴西
  plot.Brazil.mpd,plot.Brazil.mntd,
  plot.Brazil.mpd_translocation,plot.Brazil.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#俄罗斯
  plot.Russia.mpd,plot.Russia.mntd,
  plot.Russia.mpd_translocation,plot.Russia.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#墨西哥
  plot.Mexico.mpd,plot.Mexico.mntd,
  plot.Mexico.mpd_translocation,plot.Mexico.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

ggarrange(#中国
  plot.China.mpd,plot.China.mntd,
  plot.China.mpd_translocation,plot.China.mntd_translocation,
  labels=c("(i)","(ii)","(iii)","(iv)"),
  ncol=2,nrow=2,
  hjust=0,vjust=1.2,
  align="v"
)

#23.Chord diagram depicting the relationships between alien species and river basins in each country

#按照以上标准选择的国家，United States,Canada,Brazil,Russia,Mexico,China

data.used.final_chord_country = data.used.final%>%filter(X3.Native.Exotic.Status=="exotic")

#美国 大小5*5
chordDiagram(data.used.final_chord_country%>%filter(X2.Country=="United States")%>%select(c(4,2)),annotationTrack = c("grid","axis"))
data.used.final_chord_country%>%filter(X2.Country=="United States")%>%select(c(4,2))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Cyprinus_carpio被引入到最多流域，引入到 90 basins
data.used.final_chord_country%>%filter(X2.Country=="United States")%>%select(c(4,2))%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))#Caloosahatchee是引入外来种最多的流域，引入了 35 species

#加拿大 大小5*5
chordDiagram(data.used.final_chord_country%>%filter(X2.Country=="Canada")%>%select(c(4,2)),annotationTrack = c("grid","axis"))
data.used.final_chord_country%>%filter(X2.Country=="Canada")%>%select(c(4,2))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Salmo_trutta被引入到最多流域，引入到 6 basins
data.used.final_chord_country%>%filter(X2.Country=="Canada")%>%select(c(4,2))%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))#Saint.Laurent是引入外来种最多的流域，引入了 24 species

#巴西 大小5*5
chordDiagram(data.used.final_chord_country%>%filter(X2.Country=="Brazil")%>%select(c(4,2)),annotationTrack = c("grid","axis"))
data.used.final_chord_country%>%filter(X2.Country=="Brazil")%>%select(c(4,2))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Poecilia_reticulata被引入到最多流域，引入到 28 basins
data.used.final_chord_country%>%filter(X2.Country=="Brazil")%>%select(c(4,2))%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))#Paraiba.do.Sul是引入外来种最多的流域，引入了 21 species

#俄罗斯 大小5*5
chordDiagram(data.used.final_chord_country%>%filter(X2.Country=="Russia")%>%select(c(4,2)),annotationTrack = c("grid","axis"))
data.used.final_chord_country%>%filter(X2.Country=="Russia")%>%select(c(4,2))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Ctenopharyngodon_idella被引入到最多流域，引入到 5 basins
data.used.final_chord_country%>%filter(X2.Country=="Russia")%>%select(c(4,2))%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))#Kuban是引入外来种最多的流域，引入了 15 species

#墨西哥 大小5*5
chordDiagram(data.used.final_chord_country%>%filter(X2.Country=="Mexico")%>%select(c(4,2)),annotationTrack = c("grid","axis"))
data.used.final_chord_country%>%filter(X2.Country=="Mexico")%>%select(c(4,2))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Oreochromis_aureus被引入到最多流域，引入到 14 basins
data.used.final_chord_country%>%filter(X2.Country=="Mexico")%>%select(c(4,2))%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))#Yaqui是引入外来种最多的流域，引入了 14 species

#中国 大小5*5
chordDiagram(data.used.final_chord_country%>%filter(X2.Country=="China")%>%select(c(4,2)),annotationTrack = c("grid","axis"))
data.used.final_chord_country%>%filter(X2.Country=="China")%>%select(c(4,2))%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))#外来种Oreochromis_mossambicus被引入到最多流域，引入到 11 basins
data.used.final_chord_country%>%filter(X2.Country=="China")%>%select(c(4,2))%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))#Nandu是引入外来种最多的流域，引入了 6 species

#24.Calculate invasion metrics using the raw data

#(1).data.used1不排除没有外来物种的国家以及不在谱系树中的物种，即是原始数据。用以计算出国外引入的外来物种

#找出一个国家中流域间转移的外来种，即在一个流域内被定义为exotic但在其他流域是本地种
data.used1 = data.all1%>%
  group_by(X1.Basin.Name)%>%
  distinct(valid_names,.keep_all = T)%>%
  select(c(1,2,4,9))%>%
  arrange(X2.Country,desc(X3.Native.Exotic.Status),X1.Basin.Name)%>%
  ungroup()

data.distinct1 = data.used1%>%group_by(X2.Country,X3.Native.Exotic.Status)%>%
  distinct(valid_names,.keep_all = T)#一个国家内外来和本地分别独一无二的种
data.distinct1 = data.distinct1%>%ungroup()
data.distinct1%>%print(n=100)

data.dup1 = data.distinct1%>%group_by(X2.Country)%>%mutate(dup=duplicated(valid_names))%>%filter(dup)#一个国家内，外来和本地重复的种，这些种应该都是流域间转移种
data.dup1 = data.dup1%>%unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)
data.dup1%>%print(n=100)

#从data.used1中将这些流域间转移种的属性exotic改为native,从而在国外引入计算exotic-native谱系距离时不考虑这些种为exotic
data.used.final1 = 
  data.used1%>%unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  mutate(X3.Native.Exotic.Status=if_else(new.col %in% data.dup1$new.col,"native",X3.Native.Exotic.Status))%>%
  mutate(new.col=NULL)

#(2).data.used.final.translocation1是不排除没有外来物种的国家以及不在谱系树中的物种，即是原始数据。

#找出国外引入的外来种
data.foreign1 = data.used.final1%>%filter(X3.Native.Exotic.Status=="exotic")%>%
  unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)

#从data.used1中将这些国外引入的外来种的属性exotic改为native,从而在流域间转移计算exotic-native谱系距离时不考虑这些种为exotic
data.used.final.translocation1 = 
  data.used1%>%unite(new.col,X2.Country,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  mutate(X3.Native.Exotic.Status=if_else(new.col %in% data.foreign1$new.col,"native",X3.Native.Exotic.Status))%>%
  mutate(new.col=NULL)

#采用以上原始数据计算以下指标：

#最被入侵的流域
data.used1%>%filter(X3.Native.Exotic.Status=="exotic")%>%group_by(X1.Basin.Name)%>%summarise(n=n())%>%arrange(desc(n))

#入侵流域最多的物种
data.used1%>%filter(X3.Native.Exotic.Status=="exotic")%>%group_by(valid_names)%>%summarise(n=n())%>%arrange(desc(n))

#入侵其他国家最多的物种
data.used.final1%>%filter(X3.Native.Exotic.Status=="exotic")%>%
  distinct(X2.Country,valid_names)%>%
  group_by(valid_names)%>%
  summarise(n=n())%>%arrange(desc(n))

#记录的外来物种总数
data.used1%>%filter(X3.Native.Exotic.Status=="exotic")%>%distinct(valid_names)

#25.Using INLA to account for phylogenetic non-independence among species, test the relationship between phylogenetic relatedness and occurrence likelihood on a global scale.

#--------国外引入的外来物种

#从全球谱系树中提取数据中涉及的物种,构建谱系协方差矩阵，标准化，计算逆矩阵
library(ape);library(INLA)
phylo1 = keep.tip(phylo,data.distance.final$exotic_species%>%unique)
phyloMatrix = vcv(phylo1)
phyloMatrix = phyloMatrix/max(phyloMatrix)
inverse_phyloMatrix = inla.as.sparse(solve(phyloMatrix))

idx.dataframe = data.frame(exotic_species=rownames(inverse_phyloMatrix),exotic_species.idx=1:nrow(inverse_phyloMatrix))
data.distance.final.inla = full_join(data.distance.final,idx.dataframe)

inla_mpd = inla(occurrence~mpd+
              f(country,model="iid")+
              f(country:basin,model="iid")+
              f(exotic_species.idx,model="generic0",Cmatrix=inverse_phyloMatrix),
              family="binomial",data=data.distance.final.inla,
              control.compute=list(dic=T,waic=T,cpo=T))
summary(inla_mpd)

inla_mntd = inla(occurrence~mntd+
                 f(country,model="iid")+
                 f(country:basin,model="iid")+
                 f(exotic_species.idx,model="generic0",Cmatrix=inverse_phyloMatrix),
               family="binomial",data=data.distance.final.inla,
               control.compute=list(dic=T,waic=T,cpo=T))
summary(inla_mntd)

#--------流域间转移的外来物种

#从全球谱系树中提取数据中涉及的物种,构建谱系协方差矩阵，标准化，计算逆矩阵
phylo2 = keep.tip(phylo,data.distance.translocation.final$exotic_species%>%unique)
phyloMatrix1 = vcv(phylo2)
phyloMatrix1 = phyloMatrix1/max(phyloMatrix1)
inverse_phyloMatrix1 = inla.as.sparse(solve(phyloMatrix1))

idx.dataframe1 = data.frame(exotic_species=rownames(inverse_phyloMatrix1),exotic_species.idx=1:nrow(inverse_phyloMatrix1))
data.distance.translocation.final.inla = full_join(data.distance.translocation.final,idx.dataframe1)

inla_mpd.translocation = inla(occurrence~mpd+
                        f(country,model="iid")+
                        f(country:basin,model="iid")+
                        f(exotic_species.idx,model="generic0",Cmatrix=inverse_phyloMatrix1),
                        family="binomial",data=data.distance.translocation.final.inla,
                        control.compute=list(dic=T,waic=T,cpo=T))
summary(inla_mpd.translocation)

inla_mntd.translocation = inla(occurrence~mntd+
                               f(country,model="iid")+
                               f(country:basin,model="iid")+
                               f(exotic_species.idx,model="generic0",Cmatrix=inverse_phyloMatrix1),
                             family="binomial",data=data.distance.translocation.final.inla,
                             control.compute=list(dic=T,waic=T,cpo=T))
summary(inla_mntd.translocation)

#26.Conduct Structural Equation Modeling (SEM), considering the effects of local species richness and phylogenetic diversity in each basin, to test the relationship between phylogenetic distance and occurrence likelihood.

#---------------------国外引入的外来物种

#-计算本地丰富度
data.distance.final.richness = 
  data.used.final%>%filter(X3.Native.Exotic.Status=="native")%>%#Wadi.Libya的Wadi.Nashu就3个国外引入种，没有本地种😌😌,所以筛选掉了。
  group_by(X1.Basin.Name)%>%summarise(richness=n())%>%
  rename(basin=X1.Basin.Name)%>%
  right_join(data.distance.final)

#-计算本地mpd，mntd
library(picante)
com.matrix = data.used.final%>%filter(X3.Native.Exotic.Status=="native")%>%select(2,4)%>%mutate(occur=1)%>%
  pivot_wider(names_from=valid_names,values_from = occur)
com.matrix[is.na(com.matrix)] = 0
com.matrix = com.matrix%>%column_to_rownames(var="X1.Basin.Name")
com.matrix%>%view
native.mpd = mpd(com.matrix,distance_all)#NA是因为basin中就1个本地种，没有两两间的谱系距离
native.mpd.df = data.frame(basin=row.names(com.matrix),native.mpd)
native.mntd = mntd(com.matrix,distance_all)
native.mntd.df = data.frame(basin=row.names(com.matrix),native.mntd)

#-整合谱系多样性到丰富度数据库
data.distance.final.diversity = 
  data.distance.final.richness%>%
  left_join(native.mpd.df)%>%
  left_join(native.mntd.df)

data.distance.final.diversity = drop_na(data.distance.final.diversity)#去除NA
data.distance.final.diversity = data.distance.final.diversity%>%#先标准化连续自变量，使得系数可比。注意这得到的不是标准意义的std，没有标准y,但是系数也是可比的。
                               mutate(richness.scale=scale(richness)%>%as.numeric(),mpd.scale=scale(mpd)%>%as.numeric(),
                                      mntd.scale=scale(mntd)%>%as.numeric(),native.mpd.scale=scale(native.mpd)%>%as.numeric(),
                                      native.mntd.scale=scale(native.mntd)%>%as.numeric()) 
data.distance.final.diversity = as.data.frame(data.distance.final.diversity)#不转化为dataframe,在psem中会出现问题

#--开展SEM,基于brms
model1 = bf(occurrence~mpd.scale+richness.scale+native.mpd.scale+(1|country/basin)+(1|exotic_species),family=bernoulli)
model2 = bf(mpd.scale~richness.scale+native.mpd.scale+(1|country/basin)+(1|exotic_species),family=gaussian)
brm.sem = brm(model1+model2+set_rescor(FALSE),data=data.distance.final.diversity,chains=4)
summary(brm.sem)
#pp_check(brm.sem,resp="mpdscale",ndraws=100)
bayes_R2(brm.sem)#只给conditional R2,及误差和置信区间
performance::r2_bayes(brm.sem)#会给conditional and marginal R2，及置信区间

#画出效应图，插入Omnigraffle组图 大小4.44*2.18
sem.effect.data = data.frame(coef=c(-0.53,0.58,0.31,-0.06),effect=c("Direct","Direct","Direct","Indirect"),predictor=c("MPD","Richness","Native MPD","Native MPD "))
sem.effect.data = sem.effect.data%>%mutate(predictor=factor(predictor,levels = c("MPD","Richness","Native MPD","Native MPD ")))
ggplot(data=sem.effect.data,aes(x=predictor,y=coef))+
  geom_col(width=0.5,fill=NA,color="black",linewidth=0.05)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y="Standardized effects")+
  theme_test()+
  theme(plot.margin = unit(c(0.05,0.05,0.05,0.05),units="lines"))

model3 = bf(occurrence~mntd.scale+richness.scale+native.mntd.scale+(1|country/basin)+(1|exotic_species),family=bernoulli)
model4 = bf(mntd.scale~richness.scale+native.mntd.scale+(1|country/basin)+(1|exotic_species),family=gaussian)
brm.sem2 = brm(model3+model4+set_rescor(FALSE),data=data.distance.final.diversity,chains=4)
summary(brm.sem2)
bayes_R2(brm.sem2)#只给conditional R2,及误差和置信区间

#画出效应图，插入Omnigraffle组图 大小4.44*2.18
sem.effect.data2 = data.frame(coef=c(-0.72,0.39,-0.72,0.04,-0.06),effect=c("Direct","Direct","Direct","Indirect","Indirect"),predictor=c("MNTD","Richness","Native MNTD","Richness ","Native MNTD "))
sem.effect.data2 = sem.effect.data2%>%mutate(predictor=factor(predictor,levels = c("MNTD","Richness","Native MNTD","Richness ","Native MNTD ")))
ggplot(data=sem.effect.data2,aes(x=predictor,y=coef))+
  geom_col(width=0.5,fill=NA,color="black",linewidth=0.05)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y="Standardized effects")+
  theme_test()+
  theme(plot.margin = unit(c(0.05,0.05,0.05,0.05),units="lines"))

#检查共线性
glmer_mpd.diversity = glmer(occurrence~mpd.scale+richness.scale+native.mpd.scale+(1|country/basin)+(1|exotic_species),family=binomial,data=data.distance.final.diversity)
summary(glmer_mpd.diversity)
car::vif(glmer_mpd.diversity)

glmer_mntd.diversity = glmer(occurrence~mntd.scale+richness.scale+native.mntd.scale+(1|country/basin)+(1|exotic_species),family=binomial,data=data.distance.final.diversity)
summary(glmer_mntd.diversity)
car::vif(glmer_mntd.diversity)

#---------------------------流域间转移的外来物种

#-计算本地丰富度
data.distance.translocation.final.richness = 
  data.used.final.translocation%>%filter(X3.Native.Exotic.Status=="native")%>%
  group_by(X1.Basin.Name)%>%summarise(richness=n())%>%
  rename(basin=X1.Basin.Name)%>%
  right_join(data.distance.translocation.final)

#-计算本地mpd，mntd
com.matrix.translocation = data.used.final.translocation%>%filter(X3.Native.Exotic.Status=="native")%>%select(2,4)%>%mutate(occur=1)%>%
  pivot_wider(names_from=valid_names,values_from = occur)
com.matrix.translocation[is.na(com.matrix.translocation)] = 0
com.matrix.translocation = com.matrix.translocation%>%column_to_rownames(var="X1.Basin.Name")
#com.matrix.translocation%>%view
native.mpd.translocation = mpd(com.matrix.translocation,distance_all)#NA是因为basin中就1个本地种，没有两两间的谱系距离
native.mpd.df.translocation = data.frame(basin=row.names(com.matrix.translocation),native.mpd=native.mpd.translocation)
native.mntd.translocation = mntd(com.matrix.translocation,distance_all)
native.mntd.df.translocation = data.frame(basin=row.names(com.matrix.translocation),native.mntd=native.mntd.translocation)

#-整合谱系多样性到丰富度数据库
data.distance.translocation.final.diversity = 
  data.distance.translocation.final.richness%>%
  left_join(native.mpd.df.translocation)%>%
  left_join(native.mntd.df.translocation)

data.distance.translocation.final.diversity = drop_na(data.distance.translocation.final.diversity)#去除NA
data.distance.translocation.final.diversity = data.distance.translocation.final.diversity%>%#先标准化连续自变量，使得系数可比。注意这得到的不是标准意义的std，没有标准y,但是系数也是可比的。
  mutate(richness.scale=scale(richness)%>%as.numeric(),mpd.scale=scale(mpd)%>%as.numeric(),
         mntd.scale=scale(mntd)%>%as.numeric(),native.mpd.scale=scale(native.mpd)%>%as.numeric(),
         native.mntd.scale=scale(native.mntd)%>%as.numeric()) 
data.distance.translocation.final.diversity = as.data.frame(data.distance.translocation.final.diversity)#不转化为dataframe,在psem中会出现问题

#开展SEM,基于brms
model5 = bf(occurrence~mpd.scale+richness.scale+native.mpd.scale+(1|country/basin)+(1|exotic_species),family=bernoulli)
model6 = bf(mpd.scale~richness.scale+native.mpd.scale+(1|country/basin)+(1|exotic_species),family=gaussian)
brm.sem.translocation = brm(model5+model6+set_rescor(FALSE),data=data.distance.translocation.final.diversity,chains=4)
summary(brm.sem.translocation)
#pp_check(brm.sem.translocation,resp="mpdscale",ndraws=100)
bayes_R2(brm.sem.translocation)#遇到内存不足，删除上面的brm.sem, brm.sem2释放内存，可以正常给结果！

#画出效应图，插入Omnigraffle组图 大小4.44*2.18
sem.effect.data.trans = data.frame(coef=c(-0.51,0.49,0.94,-0.1),effect=c("Direct","Direct","Direct","Indirect"),predictor=c("MPD","Richness","Native MPD","Native MPD "))
sem.effect.data.trans = sem.effect.data.trans%>%mutate(predictor=factor(predictor,levels = c("MPD","Richness","Native MPD","Native MPD ")))
ggplot(data=sem.effect.data.trans,aes(x=predictor,y=coef))+
  geom_col(width=0.5,fill=NA,color="black",linewidth=0.05)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y="Standardized effects")+
  theme_test()+
  theme(plot.margin = unit(c(0.05,0.05,0.05,0.05),units="lines"))

model7 = bf(occurrence~mntd.scale+richness.scale+native.mntd.scale+(1|country/basin)+(1|exotic_species),family=bernoulli)
model8 = bf(mntd.scale~richness.scale+native.mntd.scale+(1|country/basin)+(1|exotic_species),family=gaussian)
brm.sem.translocation2 = brm(model7+model8+set_rescor(FALSE),data=data.distance.translocation.final.diversity,chains=4)
summary(brm.sem.translocation2)
bayes_R2(brm.sem.translocation2)#遇到内存不足，删除上面的brm.sem, brm.sem2释放内存，可以正常给结果！

#画出效应图，插入Omnigraffle组图 大小4.44*2.18
sem.effect.data.trans2 = data.frame(coef=c(-0.64,0.35,-1.01,0.04,-0.08),effect=c("Direct","Direct","Direct","Indirect","Indirect"),predictor=c("MNTD","Richness","Native MNTD","Richness ","Native MNTD "))
sem.effect.data.trans2 = sem.effect.data.trans2%>%mutate(predictor=factor(predictor,levels = c("MNTD","Richness","Native MNTD","Richness ","Native MNTD ")))
ggplot(data=sem.effect.data.trans2,aes(x=predictor,y=coef))+
  geom_col(width=0.5,fill=NA,color="black",linewidth=0.05)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(limits=c(-1.01,1.01),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y="Standardized effects")+
  theme_test()+
  theme(plot.margin = unit(c(0.05,0.05,0.05,0.05),units="lines"))

#检查共线性
glmer_mpd.diversity.translocation = glmer(occurrence~mpd.scale+richness.scale+native.mpd.scale+(1|country/basin)+(1|exotic_species),family=binomial,data=data.distance.translocation.final.diversity)
summary(glmer_mpd.diversity.translocation)
tab_model(glmer_mpd.diversity.translocation)
car::vif(glmer_mpd.diversity.translocation)

glmer_mntd.diversity.translocation = glmer(occurrence~mntd.scale+richness.scale+native.mntd.scale+(1|country/basin)+(1|exotic_species),family=binomial,data=data.distance.translocation.final.diversity)
summary(glmer_mntd.diversity.translocation)
tab_model(glmer_mntd.diversity.translocation)
car::vif(glmer_mntd.diversity.translocation)

#27.Create maps at the basin scale, showing the quantity, proportion, and phylogenetic distance of non-native species

library(sf)
basin.sf = read_sf("Basin042017_3119.shp")#Tedesco 2017,原始数据中3119个流域的地图
st_crs(basin.sf) = 4326

library(rnaturalearth)
W  =  ne_countries(scale = 50, returnclass = "sf")#世界地图

#-------------国外引入外来物种,采用原始数据data.used.final1

#计算国外引入外来种的数量和比例
data.basin.alien = 
  data.used.final1%>%
  group_by(X1.Basin.Name,X3.Native.Exotic.Status)%>%
  summarise(richness=n())%>%
  pivot_wider(names_from = X3.Native.Exotic.Status,values_from = richness)%>%
  replace_na(list(exotic=0,native=0))%>%
  rename(basin=X1.Basin.Name,alien.richness=exotic,native.richness=native)%>%
  mutate(ratio=alien.richness/(alien.richness+native.richness))

#将连续变量cut为分类，因为连续变量填充地图很花不好看
#data.basin.alien%>%ggplot()+geom_histogram(aes(x=alien.richness),color="white")
#data.basin.alien%>%ggplot()+geom_histogram(aes(x=ratio),color="white")
#data.basin.alien$alien.richness%>%range
#data.basin.alien$ratio%>%range

data.basin.alien = 
  data.basin.alien%>%
  mutate(alien.richness.cat=cut(alien.richness,breaks=c(0,1,5,20,36),labels=c("[0,1[","[1,5[","[5,20[","[20,35]"),right=F),
         ratio.cat=cut(ratio,breaks=c(0,0.01,0.05,0.20,1.1),labels=c("[0,0.01[","[0.01,0.05[","[0.05,0.20[","[0.20,1]"),right=F))

#将数量和比例数据整合进basin.sf
basin.alien.sf = 
  basin.sf%>%
  filter(BasinName%in%data.basin.alien$basin)%>%
  left_join(data.basin.alien,by=c("BasinName"="basin"))

#国外引入鱼类数量地图
map3 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.alien.sf,aes(fill=alien.richness.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+
  #scale_fill_brewer(name="",palette="RdYlBu",direction=-1)+#直接用RdYlBu配色没有上面挑的好
  #scale_fill_gradientn(colors=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+#做成连续色效果不好
  ggtitle("Exotic fish richness in global river basins (n = 3119 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))
#ggsave("map.png",dpi=300)

#国外引入鱼类比例地图
map4 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.alien.sf,aes(fill=ratio.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+
  ggtitle("Exotic fish percantage in global river basins (n = 3119 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#计算“出现国外引入外来鱼类”河流中alien-native平均和最近谱系距离。采用data.distance.final.richness即可，1580条河流
data.basin.alien1 = 
  data.distance.final.richness%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(mpd=mean(mpd),mntd=mean(mntd))

#将连续变量cut为分类
#data.basin.alien1%>%ggplot()+geom_histogram(aes(x=mpd),color="white")
#data.basin.alien1%>%ggplot()+geom_histogram(aes(x=mntd),color="white")
#data.basin.alien1$mpd%>%range
#data.basin.alien1$mntd%>%range
data.basin.alien1 = 
  data.basin.alien1%>%
  mutate(mpd.cat=cut(mpd,breaks=c(0,0.3,0.4,0.45,0.92),labels=c("[0.15,0.30[","[0.30,0.40[","[0.40,0.45[","[0.45,0.91]"),right=F),
         mntd.cat=cut(mntd,breaks=c(0,0.1,0.2,0.3,0.92),labels=c("[0,0.10[","[0.10,0.20[","[0.20,0.30[","[0.30,0.91]"),right=F))

#将mpd和mntd数据整合进basin.sf
basin.alien.sf1 = basin.sf%>%
  filter(BasinName%in%data.basin.alien1$basin)%>%
  left_join(data.basin.alien1,by=c("BasinName"="basin"))

#国外引入鱼类与本地种的mpd地图
map5 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.alien.sf1,aes(fill=mpd.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+
  ggtitle("Exotic-native MPD in global river basins (n = 1580 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#国外引入鱼类与本地种的mntd地图
map6 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.alien.sf1,aes(fill=mntd.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+
  ggtitle("Exotic-native MNTD in global river basins (n = 1580 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#------------流域间转移的外来物种，采用原始数据data.used.final.translocation1

#计算流域间转移外来种的数量和比例
data.basin.translocation = 
  data.used.final.translocation1%>%
  group_by(X1.Basin.Name,X3.Native.Exotic.Status)%>%
  summarise(richness=n())%>%
  pivot_wider(names_from = X3.Native.Exotic.Status,values_from = richness)%>%
  replace_na(list(exotic=0,native=0))%>%
  rename(basin=X1.Basin.Name,translocation.richness=exotic,native.richness=native)%>%
  mutate(ratio=translocation.richness/(translocation.richness+native.richness))

#将连续变量cut为分类，因为连续变量填充地图很花不好看
#data.basin.translocation%>%ggplot()+geom_histogram(aes(x=translocation.richness),color="white")
#data.basin.translocation%>%ggplot()+geom_histogram(aes(x=ratio),color="white")
#data.basin.translocation$translocation.richness%>%range
#data.basin.translocation$ratio%>%range
data.basin.translocation = 
  data.basin.translocation%>%
  mutate(translocation.richness.cat=cut(translocation.richness,breaks=c(0,1,5,20,75),labels=c("[0,1[","[1,5[","[5,20[","[20,73]"),right=F),
         ratio.cat=cut(ratio,breaks=c(0,0.01,0.05,0.20,1.1),labels=c("[0,0.01[","[0.01,0.05[","[0.05,0.20[","[0.20,1]"),right=F))

#将数量和比例数据整合进basin.sf
basin.translocation.sf = 
  basin.sf%>%
  filter(BasinName%in%data.basin.translocation$basin)%>%
  left_join(data.basin.translocation,by=c("BasinName"="basin"))

#流域间转移鱼类数量地图
map7 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.translocation.sf,aes(fill=translocation.richness.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=c("#2585a6","#90d7ec","#FEE090","#D73027"))+
  ggtitle("Translocated fish richness in global river basins (n = 3119 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#流域间转移鱼类比例地图
map8 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.translocation.sf,aes(fill=ratio.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=c("#2585a6","#90d7ec","#FEE090","#D73027"))+
  ggtitle("Translocated fish percantage in global river basins (n = 3119 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#计算“出现流域间转移鱼类”河流中translocation-native平均和最近谱系距离。采用data.distance.translocation.final.richness即可，602条河流
data.basin.translocation1 = 
  data.distance.translocation.final.richness%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(mpd=mean(mpd),mntd=mean(mntd))

#将连续变量cut为分类
#data.basin.translocation1%>%ggplot()+geom_histogram(aes(x=mpd),color="white")
#data.basin.translocation1%>%ggplot()+geom_histogram(aes(x=mntd),color="white")
#data.basin.translocation1$mpd%>%range
#data.basin.translocation1$mntd%>%range
data.basin.translocation1 = 
  data.basin.translocation1%>%
  mutate(mpd.cat=cut(mpd,breaks=c(0,0.3,0.4,0.45,0.85),labels=c("[0,0.30[","[0.30,0.40[","[0.40,0.45[","[0.45,0.83]"),right=F),
         mntd.cat=cut(mntd,breaks=c(0,0.1,0.2,0.3,0.5),labels=c("[0,0.10[","[0.10,0.20[","[0.20,0.30[","[0.30,0.48]"),right=F))

#将mpd和mntd数据整合进basin.sf
basin.translocation.sf1 = 
  basin.sf%>%
  filter(BasinName%in%data.basin.translocation1$basin)%>%
  left_join(data.basin.translocation1,by=c("BasinName"="basin"))

#流域间转移鱼类与本地种的mpd地图
map9 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.translocation.sf1,aes(fill=mpd.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=c("#2585a6","#90d7ec","#FEE090","#D73027"))+
  ggtitle("Translocation-native MPD in global river basins (n = 602 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#流域间转移鱼类与本地种的mntd地图
map10 = 
  ggplot()+
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=basin.translocation.sf1,aes(fill=mntd.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=c("#2585a6","#90d7ec","#FEE090","#D73027"))+
  ggtitle("Translocation-native MNTD in global river basins (n = 602 basins)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#合并以上3-10地图,大小16*14
library(ggpubr)
ggarrange(map3,map7,map4,map8,map5,map9,map6,map10,
          labels=c("a","e","b","f","c","g","d","h"),
          nrow=4,ncol=2,
          hjust=c(-2,-2.5,-2,-4.7,-2,-2,-2,-2),
          font.label=list(size=18),
          align="v")

#28.Create maps at the country scale, showing the quantity and proportion of non-native species

#-------------国外引入外来物种,采用原始数据data.used.final1

#计算国外引入外来种的数量和比例
data.map2 = 
  data.used.final1%>%left_join(drainage_basins)%>%
  distinct(X2.Country,X3.Native.Exotic.Status,valid_names,.keep_all = T)%>%select(c(1,3,4))%>%
  group_by(X2.Country,X3.Native.Exotic.Status)%>%summarise(n=n())%>%
  pivot_wider(names_from = X3.Native.Exotic.Status,values_from = n)%>%
  replace_na(list(exotic=0))%>%
  rename(alien.richness=exotic,native.richness=native)%>%
  mutate(ratio=alien.richness/(alien.richness+native.richness))#计算绝对数量及比例

#替换国家名
data.map2$X2.Country[!(data.map2$X2.Country %in% W$name_long)]
data.map2 = data.map2%>%mutate(X2.Country=if_else(X2.Country=="Laos","Lao PDR",
                                                 if_else(X2.Country=="Brunei","Brunei Darussalam",
                                                         if_else(X2.Country=="North Korea","Dem. Rep. Korea",
                                                                 if_else(X2.Country=="Palestina","Palestine",  
                                                                         if_else(X2.Country=="Republic of Congo","Democratic Republic of the Congo",
                                                                                 if_else(X2.Country=="Russia","Russian Federation",
                                                                                         if_else(X2.Country=="South Korea","Republic of Korea",X2.Country))))))))

#将连续变量cut为分类，因为连续变量填充地图很花不好看
#data.map2%>%ggplot()+geom_histogram(aes(x=alien.richness),color="white")
#data.map2%>%ggplot()+geom_histogram(aes(x=ratio),color="white")
#data.map2$alien.richness%>%range
#data.map2$ratio%>%range

data.map2 = 
  data.map2%>%
  mutate(alien.richness.cat=cut(alien.richness,breaks=c(0,5,15,25,85),labels=c("[0,5[","[5,15[","[15,25[","[25,84]"),right=F),
         ratio.cat=cut(ratio,breaks=c(0,0.01,0.1,0.20,0.65),labels=c("[0,0.01[","[0.01,0.10[","[0.10,0.20[","[0.20,0.65]"),right=F))

#将数量和比例数据整合进地图W
W2 = right_join(W,data.map2,by=c("name_long"="X2.Country"))                                     

#国外引入鱼类数量地图
map11 = 
  ggplot()+#大小9.05*6.28
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=W2,aes(fill=alien.richness.cat),linewidth=0.03,color="black")+#外来绝对数量
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+
  ggtitle("Exotic fish richness in the world (n = 143 countries)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#国外引入鱼类比例地图
map12 = 
  ggplot()+#大小9.05*6.28
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=W2,aes(fill=ratio.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=brewer.pal(10,"RdYlBu")[c(9,7,5,2)])+
  ggtitle("Exotic fish percentage in the world (n = 143 countries)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#------------流域间转移的外来物种，采用原始数据data.used.final.translocation1

#计算流域间转移外来种的数量和比例
data.map3 = 
  data.used.final.translocation1%>%left_join(drainage_basins)%>%
  distinct(X2.Country,X3.Native.Exotic.Status,valid_names,.keep_all = T)%>%select(c(1,3,4))%>%
  group_by(X2.Country,X3.Native.Exotic.Status)%>%summarise(n=n())%>%
  pivot_wider(names_from = X3.Native.Exotic.Status,values_from = n)%>%
  replace_na(list(exotic=0))%>%
  rename(translocation.richness=exotic,native.richness=native)%>%
  mutate(ratio=translocation.richness/native.richness)#计算绝对数量及比例，😄注意计算比例的方式分母不需要translocation.richness+native.richness,因为国内转移的物种本身一定是属于一个国家的本地种。

#替换国家名
data.map3$X2.Country[!(data.map3$X2.Country %in% W$name_long)]
data.map3 = data.map3%>%mutate(X2.Country=if_else(X2.Country=="Laos","Lao PDR",
                                                 if_else(X2.Country=="Brunei","Brunei Darussalam",
                                                         if_else(X2.Country=="North Korea","Dem. Rep. Korea",
                                                                 if_else(X2.Country=="Palestina","Palestine",  
                                                                         if_else(X2.Country=="Republic of Congo","Democratic Republic of the Congo",
                                                                                 if_else(X2.Country=="Russia","Russian Federation",
                                                                                         if_else(X2.Country=="South Korea","Republic of Korea",X2.Country))))))))

#将连续变量cut为分类，因为连续变量填充地图很花不好看
#data.map3%>%ggplot()+geom_histogram(aes(x=translocation.richness),color="white")
#data.map3%>%ggplot()+geom_histogram(aes(x=ratio),color="white")
#data.map3$translocation.richness%>%range
#data.map3$ratio%>%range

data.map3 = 
  data.map3%>%
  mutate(translocation.richness.cat=cut(translocation.richness,breaks=c(0,5,15,25,219),labels=c("[0,5[","[5,15[","[15,25[","[25,218]"),right=F),
         ratio.cat=cut(ratio,breaks=c(0,0.01,0.1,0.20,0.31),labels=c("[0,0.01[","[0.01,0.10[","[0.10,0.20[","[0.20,0.30]"),right=F))

#将数量和比例数据整合进地图W
W3 = right_join(W,data.map3,by=c("name_long"="X2.Country"))                                      

#流域间转移鱼类数量地图
map13 = 
  ggplot()+#大小9.05*6.28
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=W3,aes(fill=translocation.richness.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=c("#2585a6","#90d7ec","#FEE090","#D73027"))+
  ggtitle("Translocated fish richness in the world (n = 143 countries)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#流域间转移鱼类比例地图
map14 = 
  ggplot()+#大小9.05*6.28
  geom_sf(data=W,fill="grey",color="grey")+#和grey75一样
  geom_sf(data=W3,aes(fill=ratio.cat),linewidth=0.03,color="black")+
  coord_sf(ylim = c(-50, 80))+
  scale_fill_manual(name="",values=c("#2585a6","#90d7ec","#FEE090","#D73027"))+
  ggtitle("Translocated fish percentage in the world (n = 143 countries)")+
  theme_map(base_size = 16)+
  theme(legend.position=c(0.04,0),
        legend.background = element_blank(),
        plot.title=element_text(size=14,hjust=0.5,face="bold"))

#合并以上11-14地图,大小16*7
ggarrange(map11,map13,map12,map14,
          labels=c("a","b","c","d"),
          nrow=2,ncol=2,
          font.label=list(size=18),
          align="v")

#29.Considering the potential influence of basin area

#(1)对relatedness-occurrence关系的影响

#国外引入外来种
data.distance.final.area = data.distance.final%>%left_join(drainage_basins%>%select(c(1,9)),by=c("basin"="X1.Basin.Name"))%>%rename(area=X9.Surface.Area)

glmer_mpd.area = glmer(occurrence~mpd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd.area)

glmer_mntd.area = glmmTMB(occurrence~mntd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd.area)

#流域间转移外来种
data.distance.translocation.final.area = data.distance.translocation.final%>%left_join(drainage_basins%>%select(c(1,9)),by=c("basin"="X1.Basin.Name"))%>%rename(area=X9.Surface.Area)

glmer_mpd.translocation.area = glmer(occurrence~mpd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd.translocation.area)

glmer_mntd.translocation.area = glmer(occurrence~mntd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd.translocation.area)

#(2)对diversity-occurrence关系的影响

#国外引入外来种
data.distance.final.diversity.area = data.distance.final.diversity%>%left_join(drainage_basins%>%select(c(1,9)),by=c("basin"="X1.Basin.Name"))%>%rename(area=X9.Surface.Area)

glmer_richness.area = glmmTMB(occurrence~richness+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final.diversity.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_richness.area)

glmer_native.mpd.area = glmer(occurrence~native.mpd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final.diversity.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_native.mpd.area)

glmer_native.mntd.area = glmer(occurrence~native.mntd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final.diversity.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_native.mntd.area)

#流域间转移外来种
data.distance.translocation.final.diversity.area = data.distance.translocation.final.diversity%>%left_join(drainage_basins%>%select(c(1,9)),by=c("basin"="X1.Basin.Name"))%>%rename(area=X9.Surface.Area)

glmer_richness.translocation.area = glmmTMB(occurrence~richness+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final.diversity.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_richness.translocation.area)

glmer_native.mpd.translocation.area = glmer(occurrence~native.mpd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final.diversity.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_native.mpd.translocation.area)

glmer_native.mntd.translocation.area = glmer(occurrence~native.mntd+log(area)+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final.diversity.area)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_native.mntd.translocation.area)

#30.Re-analyze the data after defining alien and translocated species at the biogeographic scale

#--------------------分析国外引入的外来种

#(1).只利用国外引入的外来种数据开展分析，排除流域间转移的外来种。

#找出一个生态区中流域间转移的外来种，即在一个流域内被定义为exotic但在其他流域是本地种
data.used_ecoregion = data.used%>%
  left_join(drainage_basins)%>%
  select(1:5)%>%
  arrange(X3.Ecoregion,desc(X3.Native.Exotic.Status),X1.Basin.Name)

data.distinct_ecoregion = data.used_ecoregion%>%
  group_by(X3.Ecoregion,X3.Native.Exotic.Status)%>%
  distinct(valid_names,.keep_all = T)%>%ungroup()#一个生态区内外来和本地分别独一无二的种

data.dup_ecoregion = data.distinct_ecoregion%>%group_by(X3.Ecoregion)%>%filter(duplicated(valid_names))%>%ungroup()#一个流域内，外来和本地重复的种，这些种应该都是流域间转移种
data.dup_ecoregion = data.dup_ecoregion%>%unite(new.col,X3.Ecoregion,X3.Native.Exotic.Status,valid_names,remove=F)
data.dup_ecoregion%>%print(n=100)

#从data.used_ecoregion中将这些流域间转移种的属性exotic改为native,从而在国外引入计算exotic-native谱系距离时不考虑这些种为exotic
data.used.final_ecoregion = 
  data.used_ecoregion%>%unite(new.col,X3.Ecoregion,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  mutate(X3.Native.Exotic.Status=if_else(new.col %in% data.dup_ecoregion$new.col,"native",X3.Native.Exotic.Status))%>%
  mutate(new.col=NULL)

data.used.final_ecoregion%>%group_by(X3.Ecoregion,X3.Native.Exotic.Status)%>%
  distinct(valid_names,.keep_all = T)%>%ungroup()%>%
  group_by(X3.Ecoregion)%>%filter(duplicated(valid_names))#验证一下，确实生态区内没有外来和本地重复了，流域间转移的外来种已全部更改为native

data.used.final_ecoregion%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X3.Ecoregion)%>%unique()#验证下所有7个生态区都有exotic，否则没有办法算。

#总结这个最终使用的dataframe
data.used.final_ecoregion%>%summarise(n.ecoregion=n_distinct(X3.Ecoregion),n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))
data.used.final_ecoregion%>%group_by(X3.Native.Exotic.Status)%>%summarise(n.ecoregion=n_distinct(X3.Ecoregion),n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))

#(2).产生basin和species两两结合的数据框，通过tidyr::expend.
data.used.final_ecoregion_add.occur = data.used.final_ecoregion%>%
  mutate(occurrence=1,.after=valid_names)%>%
  mutate(X2.Country=NULL)%>%
  select(5,1:4)%>%
  arrange(X3.Ecoregion,X1.Basin.Name,X3.Native.Exotic.Status)#添加一列occurrence

data.expand_ecoregion = data.used.final_ecoregion%>%
  unite(new.col,X3.Native.Exotic.Status,valid_names,sep="/")%>%
  group_by(X3.Ecoregion)%>%
  tidyr::expand(X1.Basin.Name,new.col)%>%
  separate(new.col,into=c("X3.Native.Exotic.Status","valid_names"),sep="/")#产生basin和species两两结合的数据。expand()好啊，会把数据框两列中独一无二的数据两两结合

data.pair_ecoregion = left_join(data.expand_ecoregion,data.used.final_ecoregion_add.occur)%>%
  replace_na(list(occurrence=0))%>%
  arrange(X3.Ecoregion,X1.Basin.Name,X3.Native.Exotic.Status,valid_names)

data.pair_ecoregion%>%print(n=200)

#(3).产生统计模型需要的数据-计算每个生态区所有外来种与所有流域本地种的谱系距离
phylodist.function_ecoregion = function(data){
  
  data.out = NULL
  
  for(i in 1:length(unique(data$X1.Basin.Name))){
    
    basin = filter(data,X1.Basin.Name==unique(data$X1.Basin.Name)[i])
    
    exo = filter(basin,X3.Native.Exotic.Status=="exotic")%>%pull(valid_names)
    
    nat = filter(basin,X3.Native.Exotic.Status=="native"&occurrence!=0)%>%pull(valid_names)
    
    data.out.basin = list()
    
    for(j in 1:length(exo)){
      
      dist = distance_all[exo[j],nat]
      
      data.out.basin[[j]] = c(filter(basin,valid_names==exo[j])[1,1], filter(basin,valid_names==exo[j])[1,2], exo[j], filter(basin,valid_names==exo[j])[1,5],mean(dist,na.rm=T),min(dist,na.rm=T))#注意索引一定要是第一行[1,...]，因为用这个函数计算转移种-本地种谱系距离的时候，pair之后一个流域内可能有exotic和native 都是一个物种名。因为pair的数据是排过序的，exotic一定在第一行，所以索引第一行即可。也可以通过在filter里指定X3.Native.Exotic.Status=="exotic"&valid_names==exo[j]来实现。
    }
    
    data.out = append(data.out,data.out.basin) 
    
  }
  
  data.final = as.data.frame(do.call(rbind,data.out))
  names(data.final) = c("ecoregion","basin","exotic_species","occurrence","mpd","mntd")
  data.final
}

data.distance_real.ecoregion = phylodist.function_ecoregion(data.pair_ecoregion) I 

data.distance.final_real.ecoregion = data.distance_real.ecoregion%>%mutate(ecoregion=unlist(ecoregion),basin=unlist(basin),exotic_species=unlist(exotic_species),occurrence=unlist(occurrence),mpd=unlist(mpd),mntd=unlist(mntd))%>%tibble()#把每个成分的列表转为向量
data.distance.final_real.ecoregion = data.distance.final_real.ecoregion%>%filter(!is.na(mpd))#有36行mpd为NaN,mntd为Inf,排除它们。Palearctic区的Wadi.Nashu就3个国外引入种，没有本地种😌😌

#(4). 开展模型统计检验
glmer_mpd_ecoregion = glmer(occurrence~mpd+(1|ecoregion/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final_real.ecoregion)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd_ecoregion)
tab_model(glmer_mpd_ecoregion)

#手动画图
get.data_mpd_ecoregion = get_model_data(glmer_mpd_ecoregion,type="pred",terms="mpd [all]")
tibble(get.data_mpd_ecoregion)
(plot.mpd_ecoregion = 
    ggplot(data=get.data_mpd_ecoregion,aes(x=x,y=predicted))+
    geom_line(color="#00468BFF",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#00468BFF",alpha=0.2)+
    annotate(geom="text",x=c(0.527,0.506,0.486,0.548,0.565,0.505),y=c(0.025,0.022,0.019,0.016,0.013,0.01),label=c("italic(β)[mpd]==-5.34","italic(z)==-12.13","italic(P)<2e-16","italic(R^2)[marginal]==0.03","italic(R^2)[conditional]==0.72","italic(n)==126419"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.03,label="Exotic fish species",fontface="bold")+
    labs(x=NULL,y="Occurrence probability")+
    scale_y_continuous(labels=function(x) sprintf("%.2f",x),limits=c(0,0.032))+
    theme_classic()+
    theme(axis.title.y=element_text(margin = margin(0,0.4,0,0,'cm')))#控制轴标题与轴的距离
)

glmer_mntd_ecoregion = glmer(occurrence~mntd+(1|ecoregion/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final_real.ecoregion)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd_ecoregion)
tab_model(glmer_mntd_ecoregion)

#手动画图
get.data_mntd_ecoregion = get_model_data(glmer_mntd_ecoregion,type="pred",terms="mntd [all]")
tibble(get.data_mntd_ecoregion)
(plot.mntd_ecoregion = 
    ggplot(data=get.data_mntd_ecoregion,aes(x=x,y=predicted))+
    geom_line(color="#00468BFF",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#00468BFF",alpha=0.2)+
    annotate(geom="text",x=c(0.527,0.506,0.486,0.548,0.565,0.505),y=c(0.0125,0.011,0.0095,0.008,0.0065,0.005),label=c("italic(β)[mntd]==-4.49","italic(z)==-20.08","italic(P)<2e-16","italic(R^2)[marginal]==0.06","italic(R^2)[conditional]==0.72","italic(n)==126419"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.015,label="Exotic fish species",fontface="bold")+
    labs(x=NULL,y=NULL)+
    scale_y_continuous(labels=function(x) sprintf("%.3f",x),limits=c(0,0.016))+
    theme_classic()
)

#合并两个回归图,图片大小7.53*4.05
ggarrange(
  plot.mpd_ecoregion,plot.mntd_ecoregion,
  labels=c("(a)","(b)"),
  hjust=-1,vjust=1.2
)

#--------------------分析流域间转移的外来种

#(5).只利用流域间转移的外来种开展分析，排除国外引入的外来种数据。

#找出国外引入的外来种
data.foreign_ecoregion = data.used.final_ecoregion%>%filter(X3.Native.Exotic.Status=="exotic")%>%
  unite(new.col,X3.Ecoregion,X3.Native.Exotic.Status,valid_names,remove=F)

#从data.used_ecoregion中将这些国外引入的外来种的属性exotic改为native,从而在流域间尺度计算exotic-native谱系距离时不考虑这些种为exotic
data.used.final.translocation_ecoregion = 
  data.used_ecoregion%>%unite(new.col,X3.Ecoregion,X3.Native.Exotic.Status,valid_names,remove=F)%>%
  mutate(X3.Native.Exotic.Status=if_else(new.col %in% data.foreign_ecoregion$new.col,"native",X3.Native.Exotic.Status))%>%
  mutate(new.col=NULL)

ecoregion.used.final.translocation = data.used.final.translocation_ecoregion%>%filter(X3.Native.Exotic.Status=="exotic")%>%pull(X3.Ecoregion)%>%unique()#发现只有6个生态区都有转移的外来种，Oceania区没有无法算，需要去除。
data.used.final.translocation_ecoregion = data.used.final.translocation_ecoregion%>%filter(X3.Ecoregion%in%ecoregion.used.final.translocation)#过滤出这些有转移外来物种的生态区

#总结这个最终使用的dataframe
data.used.final.translocation_ecoregion%>%summarise(n.ecoregion=n_distinct(X3.Ecoregion),n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))
data.used.final.translocation_ecoregion%>%group_by(X3.Native.Exotic.Status)%>%summarise(n.ecoregion=n_distinct(X3.Ecoregion),n.country=n_distinct(X2.Country),n.basin=n_distinct(X1.Basin.Name),n.speces=n_distinct(valid_names))

#(6).产生basin和species两两结合的数据框，通过tidyr::expend.
data.used.final.translocation_ecoregion_add.occur = data.used.final.translocation_ecoregion%>%
  mutate(occurrence=1,.after=valid_names)%>%
  mutate(X2.Country=NULL)%>%
  select(5,1:4)%>%
  arrange(X3.Ecoregion,X1.Basin.Name,X3.Native.Exotic.Status)#添加一列occurrence

data.expand.translocation_ecoregion = data.used.final.translocation_ecoregion%>%
  unite(new.col,X3.Native.Exotic.Status,valid_names,sep="/")%>%
  group_by(X3.Ecoregion)%>%
  tidyr::expand(X1.Basin.Name,new.col)%>%
  separate(new.col,into=c("X3.Native.Exotic.Status","valid_names"),sep="/")#产生basin和species两两结合的数据。expand()好啊，会把数据框两列中独一无二的数据两两结合

data.pair.translocation_ecoregion = left_join(data.expand.translocation_ecoregion,data.used.final.translocation_ecoregion_add.occur)%>%
  replace_na(list(occurrence=0))%>%
  arrange(X3.Ecoregion,X1.Basin.Name,X3.Native.Exotic.Status,valid_names)

data.pair.translocation_ecoregion%>%print(n=200)

#(7).产生统计模型需要的数据-计算每个生态区所有转移外来种与本地种的谱系距离
data.distance.translocation_real.ecoregion = phylodist.function_ecoregion(data.pair.translocation_ecoregion)

data.distance.translocation.final_real.ecoregion = data.distance.translocation_real.ecoregion%>%mutate(ecoregion=unlist(ecoregion),basin=unlist(basin),exotic_species=unlist(exotic_species),occurrence=unlist(occurrence),mpd=unlist(mpd),mntd=unlist(mntd))%>%tibble()#把每个成分的列表转为向量
data.distance.translocation.final_real.ecoregion = data.distance.translocation.final_real.ecoregion%>%filter(!is.na(mpd))#有167行mpd为NaN,mntd为Inf,排除它们
data.distance.translocation.final_real.ecoregion = data.distance.translocation.final_real.ecoregion%>%filter(mntd!=0)#这个很重要，这是忽略每个国家内转移外来种出去的流域和外来种的谱系距离（数据两两配对后，流域内同种native发生，exotic没有发生的流域，即是转移外来种出去的流域。因为exotic-native同种，所以计算出来的外来种与这个流域的最近谱系距离为0，忽略mntd=0 这些行的数据即可），因为只是要算潜在接收外来种的流域与外来种的距离

#(8).开展流域间转移种的模型统计检验
glmer_mpd.translocation_ecoregion = glmer(occurrence~mpd+(1|ecoregion/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final_real.ecoregion)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd.translocation_ecoregion)
tab_model(glmer_mpd.translocation_ecoregion)

#手动画图
get.data_mpd.translocation_ecoregion = get_model_data(glmer_mpd.translocation_ecoregion,type="pred",terms="mpd [all]")
tibble(get.data_mpd.translocation_ecoregion)
(plot.mpd.translocation_ecoregion = 
    ggplot(data=get.data_mpd.translocation_ecoregion,aes(x=x,y=predicted))+
    geom_line(color="#f47920",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#f47920",alpha=0.2)+
    annotate(geom="text",x=c(0.531,0.506,0.486,0.548,0.55,0.505),y=c(0.0075,0.0065,0.0055,0.0045,0.0035,0.0025),label=c("italic(β)[mpd]==-4.19","italic(z)==-12.31","italic(P)<2e-16","italic(R^2)[marginal]==0.02","italic(R^2)[conditional]==0.80","italic(n)==281649"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.0092,label="Translocated fish species",fontface="bold")+
    labs(x="Nonnative-native MPD",y="Occurrence probability")+
    theme_classic()
)

glmer_mntd.translocation_ecoregion = glmer(occurrence~mntd+(1|ecoregion/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final_real.ecoregion)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd.translocation_ecoregion)
tab_model(glmer_mntd.translocation_ecoregion)

#手动画图
get.data_mntd.translocation_ecoregion = get_model_data(glmer_mntd.translocation_ecoregion,type="pred",terms="mntd [all]")
tibble(get.data_mntd.translocation_ecoregion)
(plot.mntd.translocation_ecoregion = 
    ggplot(data=get.data_mntd.translocation_ecoregion,aes(x=x,y=predicted))+
    geom_line(color="#f47920",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#f47920",alpha=0.2)+
    annotate(geom="text",x=c(0.531,0.506,0.486,0.548,0.55,0.505),y=c(0.0034,0.00295,0.0025,0.00205,0.0016,0.00115),label=c("italic(β)[mntd]==-4.16","italic(z)==-22.64","italic(P)<2e-16","italic(R^2)[marginal]==0.06","italic(R^2)[conditional]==0.80","italic(n)==281649"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.0042,label="Translocated fish species",fontface="bold")+
    labs(x="Nonnative-native MNTD",y=NULL)+
    ylim(0,0.0043)+
    theme_classic()
)

#合并两个回归图,图片大小7.53*4.05
ggarrange(
  plot.mpd.translocation_ecoregion,plot.mntd.translocation_ecoregion,
  labels=c("(c)","(d)"),
  hjust=-0.5,vjust=1.2
)

#合并国外引入和流域间转移的回归图.大小6.9*6.28.
ggarrange(
  plot.mpd_ecoregion,plot.mntd_ecoregion,
  plot.mpd.translocation_ecoregion,plot.mntd.translocation_ecoregion,
  labels=c("a","b","c","d"),
  ncol=2,nrow=2,
  font.label=list(size=16),
  hjust=-1,vjust=1.2
)

#31.Re-analyze the data after excluding two major families, Cyprinidae and Salmonidae. This exclusion is due to the potential distinct environmental preferences of these two families, which may lead to an overall pattern of pre-adaptation.

species_in_phylo = phylo_all_spp$Insertions_data%>%filter(insertions!="Not_inserted")%>%tibble

#--------------------分析国外引入的外来种
data.distance.final_exclude = data.distance.final%>%left_join(species_in_phylo[,1:2],by=c("exotic_species"="s"))
data.distance.final_exclude = data.distance.final_exclude%>%filter(f!="Cyprinidae"&f!="Salmonidae")

#开展模型统计检验
glmer_mpd_exclude = glmer(occurrence~mpd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final_exclude)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mpd_exclude)
tab_model(glmer_mpd_exclude)

#画图
get.data_mpd_exclude = get_model_data(glmer_mpd_exclude,type="pred",terms="mpd [all]")
tibble(get.data_mpd_exclude)
(plot.mpd_exclude = 
    ggplot(data=get.data_mpd_exclude,aes(x=x,y=predicted))+
    geom_line(color="#00468BFF",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#00468BFF",alpha=0.2)+
    annotate(geom="text",x=c(0.507,0.47,0.5,0.538,0.54,0.468),y=c(0.3,0.27,0.24,0.21,0.18,0.15),label=c("italic(β)[mpd]==-3.27","italic(z)==-7.15","italic(P)==8.86e-13","italic(R^2)[marginal]==0.007","italic(R^2)[conditional]==0.81","italic(n)==47537"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.36,label="Exotic fish species",fontface="bold")+
    labs(x=NULL,y="Occurrence probability")+
    scale_y_continuous(labels=function(x) sprintf("%.2f",x))+
    theme_classic()
)

glmer_mntd_exclude = glmer(occurrence~mntd+(1|country/basin)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.final_exclude)#有很多0，采用cloglog连接函数，不采用默认的logit
summary(glmer_mntd_exclude)
tab_model(glmer_mntd_exclude)

#画图
get.data_mntd_exclude = get_model_data(glmer_mntd_exclude,type="pred",terms="mntd [all]")
tibble(get.data_mntd_exclude)
(plot.mntd_exclude = 
    ggplot(data=get.data_mntd_exclude,aes(x=x,y=predicted))+
    geom_line(color="#00468BFF",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#00468BFF",alpha=0.2)+
    annotate(geom="text",x=c(0.5,0.48,0.462,0.52,0.523,0.467),y=c(0.235,0.21,0.185,0.16,0.135,0.11),label=c("italic(β)[mpd]==-3.95","italic(z)==-13.63","italic(P)<2e-16","italic(R^2)[marginal]==0.03","italic(R^2)[conditional]==0.80","italic(n)==47537"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.28,label="Exotic fish species",fontface="bold")+
    labs(x=NULL,y=NULL)+
    scale_y_continuous(labels=function(x) sprintf("%.2f",x),limits=c(0,0.28))+
    theme_classic()
)

#--------------------分析流域间转移的外来种
data.distance.translocation.final_exclude = data.distance.translocation.final%>%left_join(species_in_phylo[,1:2],by=c("exotic_species"="s"))
data.distance.translocation.final_exclude = data.distance.translocation.final_exclude%>%filter(f!="Cyprinidae"&f!="Salmonidae")

#开展流域间转移种的模型统计检验
glmer_mpd.translocation_exclude = glmer(occurrence~mpd+(1|country)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final_exclude)#很多的流域就只有一个数据了，basin再作为random effect时模型就不收敛了。
summary(glmer_mpd.translocation_exclude)
tab_model(glmer_mpd.translocation_exclude)

#画图
get.data_mpd.translocation_exclude = get_model_data(glmer_mpd.translocation_exclude,type="pred",terms="mpd [all]")
tibble(get.data_mpd.translocation_exclude)
(plot.mpd.translocation_exclude = 
    ggplot(data=get.data_mpd.translocation_exclude,aes(x=x,y=predicted))+
    geom_line(color="#f47920",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#f47920",alpha=0.2)+
    annotate(geom="text",x=c(0.5,0.465,0.49,0.52,0.535,0.46),y=c(0.3,0.27,0.24,0.21,0.18,0.15),label=c("italic(β)[mpd]==-2.46","italic(z)==-6.84","italic(P)==8.11e-12","italic(R^2)[marginal]==0.01","italic(R^2)[conditional]==0.71","italic(n)==57571"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.35,label="Translocated fish species",fontface="bold")+
    labs(x="Nonnative-native MPD",y="Occurrence probability")+
    scale_y_continuous(labels=function(x) sprintf("%.2f",x))+
    theme_classic()
)

glmer_mntd.translocation_exclude = glmer(occurrence~mntd+(1|country)+(1|exotic_species),family=binomial(link = "cloglog"),data=data.distance.translocation.final_exclude)#很多的流域就只有一个数据了，basin再作为random effect时模型就不收敛了。
summary(glmer_mntd.translocation_exclude)
tab_model(glmer_mntd.translocation_exclude)

#画图
get.data_mntd.translocation_exclude = get_model_data(glmer_mntd.translocation_exclude,type="pred",terms="mntd [all]")
tibble(get.data_mntd.translocation_exclude)
(plot.mntd.translocation_exclude = 
    ggplot(data=get.data_mntd.translocation_exclude,aes(x=x,y=predicted))+
    geom_line(color="#f47920",linewidth=1)+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#f47920",alpha=0.2)+
    annotate(geom="text",x=c(0.5,0.475,0.455,0.513,0.527,0.457),y=c(0.19,0.17,0.15,0.13,0.11,0.09),label=c("italic(β)[mntd]==-3.58","italic(z)==-18.62","italic(P)<2e-16","italic(R^2)[marginal]==0.06","italic(R^2)[conditional]==0.71","italic(n)==57571"),parse=T,size=3.5)+
    annotate(geom="text",x=0.5,y=0.22,label="Translocated fish species",fontface="bold")+
    labs(x="Nonnative-native MNTD",y=NULL)+
    theme_classic()
)

#合并国外引入和流域间转移的回归图.大小6.9*6.28.
ggarrange(
  plot.mpd_exclude,plot.mntd_exclude,
  plot.mpd.translocation_exclude,plot.mntd.translocation_exclude,
  labels=c("a","b","c","d"),
  ncol=2,nrow=2,
  font.label=list(size=16),
  hjust=-1,vjust=1.2
)

#32.Correlation relationships among multiple predictors

#--------------------分析国外引入的外来种
#本地谱系/物种多样性与MPD/MNTD的关系
(plot1 = 
    data.distance.final.diversity%>%
    group_by(basin)%>%
    summarise(richness=mean(richness),mpd=mean(mpd))%>%
    ggplot(aes(x=richness,y=mpd))+
    geom_point(shape=21,fill="grey",position=position_jitter(5,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native richness",y="MPD")+
    coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=70,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Exotic fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot2 = 
    data.distance.final.diversity%>%
    group_by(basin)%>%
    summarise(richness=mean(richness),mntd=mean(mntd))%>%
    ggplot(aes(x=richness,y=mntd))+
    geom_point(shape=21,fill="grey",position=position_jitter(10,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native richness",y="MNTD")+
    coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=100,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Exotic fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot3 = 
    data.distance.final.diversity%>%
    group_by(basin)%>%
    summarise(native.mpd=mean(native.mpd),mpd=mean(mpd))%>%
    ggplot(aes(x=native.mpd,y=mpd))+
    geom_point(shape=21,fill="grey",position=position_jitter(0.02,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native MPD",y="MPD")+
    coord_cartesian(ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.1,label.y=0.58)+
    theme_test()+
    facet_wrap(~"Exotic fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot4 = 
    data.distance.final.diversity%>%
    group_by(basin)%>%
    summarise(native.mntd=mean(native.mntd),mntd=mean(mntd))%>%
    ggplot(aes(x=native.mntd,y=mntd))+
    geom_point(shape=21,fill="grey",position=position_jitter(0.02,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native MNTD",y="MNTD")+
    coord_cartesian(ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.1,label.y=0.55)+
    theme_test()+
    facet_wrap(~"Exotic fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot5 = 
  data.distance.final.diversity%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),native.mpd=mean(native.mpd))%>%
  ggplot(aes(x=richness,y=native.mpd))+
  geom_point(shape=21,fill="grey",position=position_jitter(5,0.02),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Native richness",y="Native MPD")+
  coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=100,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Exotic fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot6 = 
  data.distance.final.diversity%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),native.mntd=mean(native.mntd))%>%
  ggplot(aes(x=richness,y=native.mntd))+
  geom_point(shape=21,fill="grey",position=position_jitter(5,0.02),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Native richness",y="Native MNTD")+
  coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=70,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Exotic fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

#外来物种多样性/比例与MPD/MNTD的关系
(plot7 = 
  data.distance.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),alien.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=alien.richness/(alien.richness+richness),.after=alien.richness)%>%
  ggplot(aes(x=alien.richness,y=mpd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.2),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Exotic fish richness",y="MPD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=6,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Exotic fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot8 = 
  data.distance.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),alien.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=alien.richness/(alien.richness+richness),.after=alien.richness)%>%
  ggplot(aes(x=alien.richness,y=mntd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.2),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Exotic fish richness",y="MNTD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=6,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Exotic fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot9 = 
  data.distance.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),alien.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=alien.richness/(alien.richness+richness),.after=alien.richness)%>%
  ggplot(aes(x=ratio,y=mpd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.02),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Exotic fish percentage",y="MPD")+
  coord_cartesian(ylim=c(0,0.75))+ 
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.1,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Exotic fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot10 = 
  data.distance.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),alien.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=alien.richness/(alien.richness+richness),.after=alien.richness)%>%
  ggplot(aes(x=ratio,y=mntd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.02),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Exotic fish percentage",y="MNTD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.1,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Exotic fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

#--------------------分析流域间转移的外来种
#本地谱系/物种多样性与MPD/MNTD的关系
(plot1_t = 
    data.distance.translocation.final.diversity%>%
    group_by(basin)%>%
    summarise(richness=mean(richness),mpd=mean(mpd))%>%
    ggplot(aes(x=richness,y=mpd))+
    geom_point(shape=21,fill="grey",position=position_jitter(5,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native richness",y="MPD")+
    coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=100,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Translocated fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot2_t = 
    data.distance.translocation.final.diversity%>%
    group_by(basin)%>%
    summarise(richness=mean(richness),mntd=mean(mntd))%>%
    ggplot(aes(x=richness,y=mntd))+
    geom_point(shape=21,fill="grey",position=position_jitter(10,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native richness",y="MNTD")+
    coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=70,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Translocated fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot3_t = 
    data.distance.translocation.final.diversity%>%
    group_by(basin)%>%
    summarise(native.mpd=mean(native.mpd),mpd=mean(mpd))%>%
    ggplot(aes(x=native.mpd,y=mpd))+
    geom_point(shape=21,fill="grey",position=position_jitter(0.02,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native MPD",y="MPD")+
    coord_cartesian(ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.1,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Translocated fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot4_t = 
    data.distance.translocation.final.diversity%>%
    group_by(basin)%>%
    summarise(native.mntd=mean(native.mntd),mntd=mean(mntd))%>%
    ggplot(aes(x=native.mntd,y=mntd))+
    geom_point(shape=21,fill="grey",position=position_jitter(0.02,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native MNTD",y="MNTD")+
    coord_cartesian(ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.1,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Translocated fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot5_t = 
    data.distance.translocation.final.diversity%>%
    group_by(basin)%>%
    summarise(richness=mean(richness),native.mpd=mean(native.mpd))%>%
    ggplot(aes(x=richness,y=native.mpd))+
    geom_point(shape=21,fill="grey",position=position_jitter(5,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native richness",y="Native MPD")+
    coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=100,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Translocated fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

(plot6_t = 
    data.distance.translocation.final.diversity%>%
    group_by(basin)%>%
    summarise(richness=mean(richness),native.mntd=mean(native.mntd))%>%
    ggplot(aes(x=richness,y=native.mntd))+
    geom_point(shape=21,fill="grey",position=position_jitter(5,0.02),size=1.5,alpha=0.2)+
    stat_smooth(method="lm",color="red")+
    labs(x="Native richness",y="Native MNTD")+
    coord_cartesian(xlim=c(0,500),ylim=c(0,0.75))+
    stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=70,label.y=0.6)+
    theme_test()+
    facet_wrap(~"Translocated fish species")+
    theme(strip.text=element_text(size=12,face="bold")))

#转移物种多样性/比例与MPD/MNTD的关系
(plot7_t = 
  data.distance.translocation.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),translocation.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=translocation.richness/(translocation.richness+richness),.after=translocation.richness)%>%
  ggplot(aes(x=translocation.richness,y=mpd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.2),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Translocated fish richness",y="MPD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=18,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Translocated fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot8_t = 
  data.distance.translocation.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),translocation.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=translocation.richness/(translocation.richness+richness),.after=translocation.richness)%>%
  ggplot(aes(x=translocation.richness,y=mntd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.2),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Translocated fish richness",y="MNTD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=18,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Translocated fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot9_t = 
  data.distance.translocation.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),translocation.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=translocation.richness/(translocation.richness+richness),.after=translocation.richness)%>%
  ggplot(aes(x=ratio,y=mpd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.02),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Translocated fish percentage",y="MPD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.2,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Translocated fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

(plot10_t = 
  data.distance.translocation.final.diversity%>%
  filter(occurrence==1)%>%
  group_by(basin)%>%
  summarise(richness=mean(richness),translocation.richness=n(),mpd=mean(mpd),mntd=mean(mntd),native.mpd=mean(native.mpd),native.mntd=mean(native.mntd))%>%
  mutate(ratio=translocation.richness/(translocation.richness+richness),.after=translocation.richness)%>%
  ggplot(aes(x=ratio,y=mntd))+
  geom_point(shape=21,fill="grey",position=position_jitter(0.02),size=1.5,alpha=0.2)+
  stat_smooth(method="lm",color="red")+
  labs(x="Translocated fish percentage",y="MNTD")+
  coord_cartesian(ylim=c(0,0.75))+
  stat_cor(aes(label=paste(after_stat(r.label),after_stat(p.label),sep="~','~")),label.x=0.2,label.y=0.6)+
  theme_test()+
  facet_wrap(~"Translocated fish species")+
  theme(strip.text=element_text(size=12,face="bold")))

#合并以上相关图，大小14*10
ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,
          plot1_t,plot2_t,plot3_t,plot4_t,plot5_t,plot6_t,plot7_t,plot8_t,plot9_t,plot10_t,
          labels=letters[1:20],
          label.x=c(0,0,0,0,0,0.05,0,0,0.05,0.05,0,0.05,-0.05,0,0,0,0,0.05,0,0.05),
          ncol=5,nrow=4,
          align="hv",
          hjust=-2,vjust=1.5)
