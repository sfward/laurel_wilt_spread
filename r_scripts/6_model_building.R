# Load packages  --------------------------------------------
#
#
#---
source("r_scripts/packages.R")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# CPH models  --------------------------------------------
#
#
#---
CPH_df <- fread("data/CPH_data/CPH_LWD_COMBINED.csv", na=".")

# fix the FIPS (change from 4 characters to 5)
CPH_df$FIPS <- ifelse(nchar(CPH_df$FIPS)==4, paste(0,CPH_df$FIPS,sep=""), paste(CPH_df$FIPS))
table(nchar(CPH_df$FIPS))
summary(CPH_df)

CPH_df_raw <- CPH_df
table((CPH_df_raw$FIPS))

# start in 2005
CPH_df <- CPH_df[which(CPH_df$time1 >= 2006),]


# transform relevant data
CPH_df$ln_prop_press_dist2 <- log(CPH_df$prop_press_dist2)
hist(CPH_df$prop_press_dist2)
hist(CPH_df$ln_prop_press_dist2)

#
CPH_df$ln_totpop <- log(CPH_df$POPESTIMATE2019)
hist(CPH_df$ln_totpop)
#
CPH_df$ln_camp <- log(CPH_df$n_campgrounds+1)
hist(CPH_df$ln_camp)
#
CPH_df$ln_med_income <- log(CPH_df$med_income)
hist(CPH_df$ln_med_income)
#
CPH_df$ln_host <- log(CPH_df$host_biomass+1)
hist(CPH_df$ln_host)
#
CPH_df$ln_host_red <- log(CPH_df$redbay_biomass+1)
hist(CPH_df$ln_host_red)
#
#
CPH_df$ln_host_sas <- log(CPH_df$sassafras_biomass+1)
hist(CPH_df$ln_host_sas)
#
CPH_df$ln_non_host<- log(CPH_df$non_host_biomass+1)
hist(CPH_df$ln_non_host)

names(CPH_df)
CPH_df <- CPH_df[, c("FIPS","State", "NAME", "Invaded","YrInv", "time0", "time1", "invaded_time1",
                     "ln_prop_press_dist2", "ln_totpop", "ln_camp","ln_med_income",
                     "ln_host", "ln_host_red", "ln_host_sas", "ln_non_host","MINMOTEMP","MAP")]
CPH_df_scaled <- CPH_df[,c("ln_prop_press_dist2", "ln_totpop", "ln_camp","ln_med_income",
                     "ln_host", "ln_host_red", "ln_host_sas", "ln_non_host","MINMOTEMP","MAP")]
CPH_df_scaled <- scale(CPH_df_scaled)
colnames(CPH_df_scaled) <- c("Contagion","Humans","Campgrounds","Income",
                             "Hosts", "Redbay", "Sassafras", "Nonhosts","MinTemp","Precip")

CPH_df <- cbind.data.frame(CPH_df, CPH_df_scaled)




# remove the 2021-2022 interval, meaning observations for counties that were invaded in 2021
# note these counties still appear in the dataset as uninvaded (as of 2020)
# View(head(CPH_df,500))
# View(CPH_df[which(CPH_df$YrInv %in% 2021),])
nrow(CPH_df) 

CPH_df.first <- CPH_df[match(unique(CPH_df$FIPS), CPH_df$FIPS),]
training <- data.frame(table(CPH_df.first$YrInv))
colnames(training) <- c("year", "n_counties")
training$n_cumul <- cumsum(training$n_counties)
training$percent_inv <- training$n_cumul/(275-3)
(30+33+16)/275

#
CPH_df_train <- CPH_df[which(CPH_df$time1 <= 2019),]
CPH_df_train[which(CPH_df_train$YrInv == 2017),]
CPH_df_train[which(CPH_df_train$YrInv == 2018),]
CPH_df_train[which(CPH_df_train$YrInv == 2019),] # check that counties invaded in 2019 are not marked invaded in the dataset
table(CPH_df_train$time1)
nrow(CPH_df_train) # original dataset ending in 2021 (invasion data 2004-2020) was 51231 rows
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




# Create three data frames --------------------------------------------
#
#
#---
# Model 1
CPH_df_hosts <- CPH_df_train[which(CPH_df_train$ln_host > 0),] # subset counties to those with hosts
length(unique(CPH_df_hosts$FIPS))

# Model 2
# make subset to buffered area (counties_buff_500km), made in "estimating spread" R script
counties_within_buffer = st_read("shapefiles/invasion_data/buffered_invasion.shp")
CPH_df_buffer <- CPH_df_train[which(CPH_df_train$FIPS %in% counties_within_buffer$FIPS),]
nrow(counties_within_buffer)
length(unique(CPH_df_buffer$FIPS))


# Model 3
# uses all of the data, so no subset
length(unique(CPH_df_train$FIPS))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# CPH models COUNTIES WITH HOSTS --------------------------------------------
#
#
#---
matrix_cor_hosts <-cor(CPH_df_hosts[, c("Contagion","Humans","Campgrounds","Income",
                                        "Hosts","Nonhosts","MinTemp","Precip")])
matrix_cor_hosts[matrix_cor_hosts>0.7]

summary(CPH_df_hosts)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# full model
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
fit_lwd_full_host <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                           Hosts + Nonhosts + MinTemp + Precip,
                      data=CPH_df_hosts)
summary(fit_lwd_full_host)

options(na.action = "na.fail")
dd_host <- dredge(fit_lwd_full_host, rank=AIC)
subset(dd_host, delta < 10)
write.excel(subset(dd_host, delta < 10))

model_averaged_host <- model.avg(dd_host, subset = delta < 2)
write.excel(cbind(row.names(confint(model_averaged_host)), confint(model_averaged_host)))


# reduced 
fit_lwd_red_host <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Nonhosts + MinTemp,
                           data=CPH_df_hosts)
summary(fit_lwd_red_host)
AIC(fit_lwd_red_host)
round(summary(fit_lwd_red_host)$coef,4)
round(summary(fit_lwd_red_host)$coef,2)
write.excel(cbind(rownames(summary(fit_lwd_red_host)$coef),round(summary(fit_lwd_red_host)$coef,4)))
write.excel(cbind(rownames(summary(fit_lwd_red_host)$coef),round(summary(fit_lwd_red_host)$coef,2)))

# https://stackoverflow.com/questions/54962119/how-to-plot-from-mumin-model-avg-summary
options(na.action = "na.fail") # needed for dredge to work
mA_host <-summary(model_averaged_host) #pulling out model averages
df1_host<-as.data.frame(mA_host$coefmat.full) #selecting full model coefficient averages
CI_host <- as.data.frame(confint(model_averaged_host, full=T)) # get confidence intervals for full model
df1_host$CI.min <-CI_host$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1_host$CI.max <-CI_host$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1_host, keep.rownames = "coefficient") #put rownames into column
names(df1_host) <- gsub(" ", "", names(df1_host)) # remove spaces from column headers
#
CI_host_red <- cbind(data.frame(row.names(confint(fit_lwd_red_host))),
                     data.frame(coef(fit_lwd_red_host)),
                     data.frame(confint(fit_lwd_red_host))) # get confidence intervals for reduced model
colnames(CI_host_red) <- c("coefficient", "Estimate", "CI.min" ,"CI.max")

every_variable <- row.names(confint(fit_lwd_full_host)) 
missing_to_add <- every_variable[which(every_variable %!in% df1_host$coefficient )]
if(length(missing_to_add) > 0){
missing_to_add <- data.frame(coefficient=missing_to_add, Estimate=NA, CI.min=NA, CI.max=NA)
CI_host_red <- rbind.data.frame(CI_host_red,missing_to_add)} else{
  CI_host_red <- CI_host_red}


df1_host$type <- "Model-averaged"
CI_host_red$type <- "Lowest AIC"
df_fin_host <- rbind(df1_host[,c("coefficient", "Estimate", "CI.min" ,"CI.max", "type")],
                     CI_host_red)
#
df_fin_host$coefficient <- factor(df_fin_host$coefficient, levels=rev(c("Contagion", "Humans", "Campgrounds","Income",
                                                                            "Hosts", "Nonhosts", "MinTemp", "Precip")))

df_fin_host$type <- factor(df_fin_host$type, levels = c("Model-averaged", "Lowest AIC"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''











# CPH models COUNTIES WITH BUFFER --------------------------------------------
#
#
#---
matrix_cor_buffer <-cor(CPH_df_buffer[, c("Contagion","Humans","Campgrounds","Income",
                                        "Hosts","Nonhosts","MinTemp","Precip")])
matrix_cor_buffer[matrix_cor_buffer>0.7]

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# full model
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
fit_lwd_full_buffer <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                               Hosts + Nonhosts + MinTemp + Precip,
                      data=CPH_df_buffer)
summary(fit_lwd_full_buffer)
AIC(fit_lwd_full_buffer)

dd_buffer <- dredge(fit_lwd_full_buffer)
subset(dd_buffer, delta < 2)
write.excel(subset(dd_buffer, delta < 2))

model_averaged_buffer <- model.avg(dd_buffer, subset = delta < 2)# AQUI 
model_averaged_buffer <- model.avg(rbind(dd_buffer[1,],dd_buffer[1,])) # AQUI 
write.excel(cbind(row.names(confint(model_averaged_buffer)), confint(model_averaged_buffer)))

# reduced 
fit_lwd_red_buffer <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                              Hosts + MinTemp + Precip,
                          data=CPH_df_buffer)
summary(fit_lwd_red_buffer)
AIC(fit_lwd_red_buffer)

write.excel(cbind(rownames(summary(fit_lwd_red_buffer)$coef),round(summary(fit_lwd_red_buffer)$coef,4)))
write.excel(cbind(rownames(summary(fit_lwd_red_buffer)$coef),round(summary(fit_lwd_red_buffer)$coef,2)))

# https://stackoverflow.com/questions/54962119/how-to-plot-from-mumin-model-avg-summary
options(na.action = "na.fail") # needed for dredge to work
mA_buffer <-summary(model_averaged_buffer) #pulling out model averages
df1_buffer<-as.data.frame(mA_buffer$coefmat.full) #selecting full model coefficient averages
CI_buffer <- as.data.frame(confint(model_averaged_buffer, full=T)) # get confidence intervals for full model
df1_buffer$CI.min <-CI_buffer$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1_buffer$CI.max <-CI_buffer$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1_buffer, keep.rownames = "coefficient") #put rownames into column
names(df1_buffer) <- gsub(" ", "", names(df1_buffer)) # remove spaces from column headers
#
CI_buffer_red <- cbind(data.frame(row.names(confint(fit_lwd_red_buffer))),
                     data.frame(coef(fit_lwd_red_buffer)),
                     data.frame(confint(fit_lwd_red_buffer))) # get confidence intervals for reduced model
colnames(CI_buffer_red) <- c("coefficient", "Estimate", "CI.min" ,"CI.max")

every_variable <- row.names(confint(fit_lwd_full_buffer)) 
missing_to_add <- every_variable[which(every_variable %!in% df1_buffer$coefficient )]
missing_to_add <- data.frame(coefficient=missing_to_add, Estimate=NA, CI.min=NA, CI.max=NA)
CI_buffer_red <- rbind.data.frame(CI_buffer_red,missing_to_add)


df1_buffer$type <- "Model-averaged"
CI_buffer_red$type <- "Lowest AIC"
df_fin_buffer <- rbind(df1_buffer[,c("coefficient", "Estimate", "CI.min" ,"CI.max", "type")],
                     CI_buffer_red)
#
df_fin_buffer$coefficient <- factor(df_fin_buffer$coefficient, levels=rev(c("Contagion", "Humans", "Campgrounds","Income",
                                                                      "Hosts", "Nonhosts", "MinTemp", "Precip")))

df_fin_buffer$type <- factor(df_fin_buffer$type, levels = c("Model-averaged", "Lowest AIC"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''














# CPH models CONTIGUOUS USA --------------------------------------------
#
#
#---
matrix_cor <-cor(CPH_df[, c("Contagion","Humans","Campgrounds","Income",
                                          "Hosts","Nonhosts","MinTemp","Precip")])
matrix_cor[matrix_cor>0.7]



#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# full model
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
fit_lwd_full_all <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                            Hosts + Nonhosts + MinTemp + Precip,
                      data=CPH_df_train)
summary(fit_lwd_full_all)
AIC(fit_lwd_full_all)

dd_all <- dredge(fit_lwd_full_all)
subset(dd_all, delta < 2)
write.excel(subset(dd_all, delta < 2))

model_averaged_all<- model.avg(dd_all, subset = delta < 2) # AQUI 
model_averaged_all <- model.avg(rbind(dd_all[1,],dd_all[1,])) # AQUI 
write.excel(cbind(row.names(confint(model_averaged_all)), confint(model_averaged_all)))


# reduced 
fit_lwd_red_all <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                           Hosts + MinTemp + Precip,
                            data=CPH_df_train)
summary(fit_lwd_red_all)
AIC(fit_lwd_red_all)

write.excel(cbind(rownames(summary(fit_lwd_red_all)$coef),round(summary(fit_lwd_red_all)$coef,4)))
write.excel(cbind(rownames(summary(fit_lwd_red_all)$coef),round(summary(fit_lwd_red_all)$coef,2)))


# https://stackoverflow.com/questions/54962119/how-to-plot-from-mumin-model-avg-summary
options(na.action = "na.fail") # needed for dredge to work
mA_all <-summary(model_averaged_all) #pulling out model averages
df1_all<-as.data.frame(mA_all$coefmat.full) #selecting full model coefficient averages
CI_all <- as.data.frame(confint(model_averaged_all, full=T)) # get confidence intervals for full model
df1_all$CI.min <-CI_all$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1_all$CI.max <-CI_all$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1_all, keep.rownames = "coefficient") #put rownames into column
names(df1_all) <- gsub(" ", "", names(df1_all)) # remove spaces from column headers
#
CI_all_red <- cbind(data.frame(row.names(confint(fit_lwd_red_all))),
                       data.frame(coef(fit_lwd_red_all)),
                       data.frame(confint(fit_lwd_red_all))) # get confidence intervals for reduced model
colnames(CI_all_red) <- c("coefficient", "Estimate", "CI.min" ,"CI.max")

every_variable <- row.names(confint(fit_lwd_full_all)) 
missing_to_add <- every_variable[which(every_variable %!in% df1_all$coefficient )]
missing_to_add <- data.frame(coefficient=missing_to_add, Estimate=NA, CI.min=NA, CI.max=NA)
CI_all_red <- rbind.data.frame(CI_all_red,missing_to_add)


df1_all$type <- "Model-averaged"
CI_all_red$type <- "Lowest AIC"
df_fin_all <- rbind(df1_all[,c("coefficient", "Estimate", "CI.min" ,"CI.max", "type")],
                       CI_all_red)
#
df_fin_all$coefficient <- factor(df_fin_all$coefficient, levels=rev(c("Contagion", "Humans", "Campgrounds","Income",
                                                                        "Hosts", "Nonhosts", "MinTemp", "Precip")))

df_fin_all$type <- factor(df_fin_all$type, levels = c("Model-averaged", "Lowest AIC"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''








# Figure: coefficient graph --------------------------------------------
#
#
#---
dodge <- position_dodge(width=0.6)  
base_size_text <- 12
point_size = 1.9


host_coef <- ggplot(data=df_fin_host, aes(x=coefficient, y=Estimate, color=type))+ #again, excluding intercept because estimates so much larger
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1)+ #add dashed line at zero
  #
  geom_point(position=dodge, size=point_size) +
  geom_errorbar(aes(ymax=CI.max,ymin=CI.min),position = dodge, width=0, size=0) + 
  scale_colour_manual(values = c(viridis(6)[3],viridis(6)[1])) +
  #
  coord_flip()+ # flipping x and y axes
  theme_classic(base_size = base_size_text)+ ylab("")+ xlab("")+
  theme(legend.position =  "none", legend.title = element_blank())+
  scale_y_continuous(limits=c(-1, 6))

buffer_coef <- ggplot(data=df_fin_buffer, aes(x=coefficient, y=Estimate, color=type))+ #again, excluding intercept because estimates so much larger
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1)+ #add dashed line at zero
  #
  geom_point(position=dodge, size=point_size) +
  geom_errorbar(aes(ymax=CI.max,ymin=CI.min),position = dodge, width=0, size=0) + 
  scale_colour_manual(values = c(viridis(6)[3],viridis(6)[1])) +
  #
  coord_flip()+ # flipping x and y axes
  theme_classic(base_size = base_size_text)+ ylab("")+ xlab("")+
  theme(legend.position = "none", legend.title = element_blank())+
  scale_y_continuous(limits=c(-1, 6))

all_coef <- ggplot(data=df_fin_all, aes(x=coefficient, y=Estimate, color=type))+ #again, excluding intercept because estimates so much larger
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1)+ #add dashed line at zero
  #
  geom_point(position=dodge, size=point_size) +
  geom_errorbar(aes(ymax=CI.max,ymin=CI.min),position = dodge, width=0, size=0) + 
  scale_colour_manual(values = c(viridis(6)[3],viridis(6)[1])) +
  guides(color = guide_legend(reverse = T))+
  #
  coord_flip()+ # flipping x and y axes
  theme_classic(base_size = base_size_text)+ ylab("Slope coefficient ± 95% CL")+ xlab("")+
  theme(legend.position = c(0.8, 0.3), legend.title = element_blank())+
  scale_y_continuous(limits=c(-1, 6))

resize.win(6.85,8)
ggarrange(host_coef, buffer_coef, all_coef, ncol = 1, nrow = 3,  align = "hv",
          labels="auto", label.x =0.03, label.y =1)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




beep(1)



# Forecasts --------------------------------------------
#
#
#---
#'''''''''''''''''''''''''
# Load in shapefile for mapping
#'''''''''''''''''''''''''
spread_ALBERS = st_read("shapefiles/invasion_data/uscounties_ALBERS.shp")
spread_ALBERS_hosts_crop <- spread_ALBERS[which(spread_ALBERS$FIPS %in% CPH_df_buffer$FIPS),]
spread_ALBERS_hosts <- st_crop(spread_ALBERS, spread_ALBERS_hosts_crop)
#
statesUS = st_read("shapefiles/states/cb_2018_us_state_20m.shp")
statesUS_albers <- st_transform(statesUS, st_crs(spread_ALBERS_hosts_crop))
#
statesUS_albers_hosts <- st_crop(statesUS_albers, spread_ALBERS_hosts_crop)
nrow(spread_ALBERS_hosts)


#'''''''''''''''''''''''''
# Hosts-only model
#'''''''''''''''''''''''''
# CPH_df_hosts_predict <- CPH_df[which(CPH_df$time1 %in% 2019 & CPH_df$YrInv %in% c(NA,2021)),]
# CPH_df_hosts_predict$predicted_vals_hosts <- 1-exp(-predict(fit_lwd_red_host, newdata=CPH_df_hosts_predict,type=c("expected")))
# nrow(CPH_df_hosts_predict)
# summary(CPH_df_hosts_predict)

CPH_df_hosts_predict_2019 <- CPH_df[which(CPH_df$time1 %in% 2020 & CPH_df$YrInv %in% c(NA,2019,2020,2021)),]
CPH_df_hosts_predict_2019$time0 <- 2018
CPH_df_hosts_predict_2019$time1 <- 2019
CPH_df_hosts_predict_2019$year_predict <- 2019
CPH_df_hosts_predict_2019$predicted_vals_hosts <- 1-exp(-predict(fit_lwd_red_host, newdata=CPH_df_hosts_predict_2019,type=c("expected")))
summary(CPH_df_hosts_predict_2019)
#

CPH_df_hosts_predict_2020 <- CPH_df[which(CPH_df$time1 %in% 2021 & CPH_df$YrInv %in% c(NA,2020,2021)),]
CPH_df_hosts_predict_2020$time0 <- 2018
CPH_df_hosts_predict_2020$time1 <- 2019
CPH_df_hosts_predict_2020$year_predict <- 2020
CPH_df_hosts_predict_2020$predicted_vals_hosts <- 1-exp(-predict(fit_lwd_red_host, newdata=CPH_df_hosts_predict_2020,type=c("expected")))
#
CPH_df_hosts_predict_2021 <- CPH_df[which(CPH_df$time1 %in% 2022 & CPH_df$YrInv %in% c(NA,2021)),]
CPH_df_hosts_predict_2021$time0 <- 2018
CPH_df_hosts_predict_2021$time1 <- 2019
CPH_df_hosts_predict_2021$year_predict <- 2021
CPH_df_hosts_predict_2021$predicted_vals_hosts <- 1-exp(-predict(fit_lwd_red_host, newdata=CPH_df_hosts_predict_2021,type=c("expected")))


CPH_df_hosts_predict <- rbind(CPH_df_hosts_predict_2019,CPH_df_hosts_predict_2020,CPH_df_hosts_predict_2021)


#'''''''''''''''''''''''''
# Buffer model
#'''''''''''''''''''''''''
CPH_df_buffer_predict_2019 <- CPH_df[which(CPH_df$time1 %in% 2020 & CPH_df$YrInv %in% c(NA,2019,2020,2021)),]
CPH_df_buffer_predict_2019$time0 <- 2018
CPH_df_buffer_predict_2019$time1 <- 2019
summary(CPH_df_buffer_predict_2019)
CPH_df_buffer_predict_2019$predicted_vals_buffer <- 1-exp(-predict(fit_lwd_red_buffer, newdata=CPH_df_buffer_predict_2019,type=c("expected")))
#

CPH_df_buffer_predict_2020 <- CPH_df[which(CPH_df$time1 %in% 2021 & CPH_df$YrInv %in% c(NA,2020,2021)),]
CPH_df_buffer_predict_2020$time0 <- 2018
CPH_df_buffer_predict_2020$time1 <- 2019
CPH_df_buffer_predict_2020$predicted_vals_buffer <- 1-exp(-predict(fit_lwd_red_buffer, newdata=CPH_df_buffer_predict_2020,type=c("expected")))

#
CPH_df_buffer_predict_2021 <- CPH_df[which(CPH_df$time1 %in% 2022 & CPH_df$YrInv %in% c(NA,2021)),]
CPH_df_buffer_predict_2021$time0 <- 2018
CPH_df_buffer_predict_2021$time1 <- 2019
CPH_df_buffer_predict_2021$predicted_vals_buffer <- 1-exp(-predict(fit_lwd_red_buffer, newdata=CPH_df_buffer_predict_2021,type=c("expected")))

CPH_df_buffer_predict <- rbind(CPH_df_buffer_predict_2019,CPH_df_buffer_predict_2020,CPH_df_buffer_predict_2021)

#'''''''''''''''''''''''''
# All data model
#'''''''''''''''''''''''''
CPH_df_predict_2019 <- CPH_df[which(CPH_df$time1 %in% 2020 & CPH_df$YrInv %in% c(NA,2019,2020,2021)),]
CPH_df_predict_2019$time0 <- 2018
CPH_df_predict_2019$time1 <- 2019
CPH_df_predict_2019$predicted_vals_all <- 1-exp(-predict(fit_lwd_red_all, newdata=CPH_df_predict_2019,type=c("expected")))
#

CPH_df_predict_2020 <- CPH_df[which(CPH_df$time1 %in% 2021 & CPH_df$YrInv %in% c(NA,2020,2021)),]
CPH_df_predict_2020$time0 <- 2018
CPH_df_predict_2020$time1 <- 2019
CPH_df_predict_2020$predicted_vals_all <- 1-exp(-predict(fit_lwd_red_all, newdata=CPH_df_predict_2020,type=c("expected")))
#
CPH_df_predict_2021 <- CPH_df[which(CPH_df$time1 %in% 2022 & CPH_df$YrInv %in% c(NA,2021)),]
CPH_df_predict_2021$time0 <- 2018
CPH_df_predict_2021$time1 <- 2019
CPH_df_predict_2021$predicted_vals_all <- 1-exp(-predict(fit_lwd_red_all, newdata=CPH_df_predict_2021,type=c("expected")))

CPH_df_predict <- rbind(CPH_df_predict_2019,CPH_df_predict_2020,CPH_df_predict_2021)

1-exp(-predict(fit_lwd_red_all, newdata=CPH_df_predict_2021,type=c("expected")))==
1-predict(fit_lwd_red_all, newdata=CPH_df_predict_2021,type=c("survival"))




#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

nrow(CPH_df_hosts_predict)
nrow(CPH_df_buffer_predict)
nrow(CPH_df_predict)
CPH_df_buffer_predict$FIPS_BUFFER <- CPH_df_buffer_predict$FIPS
CPH_df_predict$FIPS_ALL <- CPH_df_predict$FIPS
all_hosts_buffer_forecast <- cbind(CPH_df_hosts_predict[,c("FIPS", "predicted_vals_hosts", "YrInv", "year_predict")],
                               CPH_df_buffer_predict[,c("FIPS_BUFFER", "predicted_vals_buffer")], CPH_df_predict[,c("FIPS_ALL", "predicted_vals_all")])
nrow(all_hosts_buffer_forecast)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




# Figure: Forecast comparison using each model --------------------------------------------
# 
#
#---
theme_comps <-   theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))

comp_p1 <- ggplot(all_hosts_buffer_forecast, aes(x=predicted_vals_hosts, y=predicted_vals_buffer)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  xlab("Hosts")+ylab("Buffer")+
  theme_bw() + geom_abline(intercept = 0, slope = 1, linetype="dashed")+theme_comps
comp_p2 <- ggplot(all_hosts_buffer_forecast, aes(x=predicted_vals_hosts, y=predicted_vals_all)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  xlab("Hosts")+ylab("All")+
  theme_bw() + geom_abline(intercept = 0, slope = 1, linetype="dashed")+theme_comps
comp_p3 <- ggplot(all_hosts_buffer_forecast, aes(x=predicted_vals_buffer, y=predicted_vals_all)) +
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  xlab("Buffer")+ylab("All")+
  theme_bw() + geom_abline(intercept = 0, slope = 1, linetype="dashed")+theme_comps


resize.win(9.85,3.2)
ggarrange(comp_p1, comp_p2, comp_p3, ncol = 3, nrow = 1,  align = "hv",
          labels="auto", label.x =0.18, label.y =1)
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# Check estimated vs. observed --------------------------------------------
#
#
#---
# from Lynn: If the lw_detection_year is blank, it hasn't been found there yet. 
LWD_cty <- read_excel(path="data/LWD_county_level/Laurel_Wilt_Counties_2_9_2022.xlsx",
                      sheet = "Laurel_Wilt_Counties_2_9_2022",
                      range = "A1:H1264", col_names=T)
table(nchar(LWD_cty$FIPS))
LWD_cty$FIPS <- ifelse(nchar(LWD_cty$FIPS) == 4, paste("0", LWD_cty$FIPS, sep=""), paste(LWD_cty$FIPS))
table(nchar(LWD_cty$FIPS))


# get 2019 data
LWD_cty_2019 <- LWD_cty[which(LWD_cty$lw_detection_year == 2019),]
nrow(LWD_cty_2019)
projected_2019 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$FIPS %in% LWD_cty_2019$FIPS),]
projected_2019_no_inv <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$FIPS %!in% LWD_cty_2019$FIPS & is.na(all_hosts_buffer_forecast$YrInv)),]

all_hosts_buffer_forecast$Invaded_2019 <- ifelse(all_hosts_buffer_forecast$FIPS %in% LWD_cty_2019$FIPS, "Yes", "No")
table(all_hosts_buffer_forecast$Invaded_2019)

# get 2020 data
LWD_cty_2020 <- LWD_cty[which(LWD_cty$lw_detection_year == 2020),]
nrow(LWD_cty_2020)
projected_2020 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$FIPS %in% LWD_cty_2020$FIPS),]
projected_2020_no_inv <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$FIPS %!in% LWD_cty_2020$FIPS & is.na(all_hosts_buffer_forecast$YrInv)),]

all_hosts_buffer_forecast$Invaded_2020 <- ifelse(all_hosts_buffer_forecast$FIPS %in% LWD_cty_2020$FIPS, "Yes", "No")
table(all_hosts_buffer_forecast$Invaded_2020)

# get 2021 data
LWD_cty_2021 <- LWD_cty[which(LWD_cty$lw_detection_year == 2021),]
nrow(LWD_cty_2021)
projected_2021 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$FIPS %in% LWD_cty_2021$FIPS),]
projected_2021_no_inv <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$FIPS %!in% LWD_cty_2021$FIPS & is.na(all_hosts_buffer_forecast$YrInv)),]

all_hosts_buffer_forecast$Invaded_2021 <- ifelse(all_hosts_buffer_forecast$FIPS %in% LWD_cty_2021$FIPS, "Yes", "No")
table(all_hosts_buffer_forecast$Invaded_2021)
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''






# Density plots of forecasts --------------------------------------------
# Prob of invasion for counties invaded in 2021 vs. not invaded in 2021
# http://www.sthda.com/english/wiki/ggplot2-density-plot-quick-start-guide-r-software-and-data-visualization
#---
Invaded_2019 <- data.frame(Probability = projected_2019$predicted_vals_all, type="Invaded")
Non_invaded_2019 <- data.frame(Probability = projected_2019_no_inv$predicted_vals_all, type="Non-invaded")
density_df_2019 <- rbind.data.frame(Invaded_2019,Non_invaded_2019)
#
Invaded_2020 <- data.frame(Probability = projected_2020$predicted_vals_all, type="Invaded")
Non_invaded_2020 <- data.frame(Probability = projected_2020_no_inv$predicted_vals_all, type="Non-invaded")
density_df_2020 <- rbind.data.frame(Invaded_2020,Non_invaded_2020)
#
Invaded_2021 <- data.frame(Probability = projected_2021$predicted_vals_all, type="Invaded")
Non_invaded_2021 <- data.frame(Probability = projected_2021_no_inv$predicted_vals_all, type="Non-invaded")
density_df_2021 <- rbind.data.frame(Invaded_2021,Non_invaded_2021)

density_df <- rbind.data.frame(Invaded_2019,Non_invaded_2019, Invaded_2020,Non_invaded_2020, Invaded_2021,Non_invaded_2021)



resize.win(6.85,4)
forecast_dens_cols <- viridis(12)[c(5,10)]
p_est_obs <-ggplot(density_df, aes(x=Probability, y=..scaled.., fill=type)) +
  geom_density(alpha=0.4) +theme_bw()+theme_comps+
  scale_fill_manual(values=c(viridis(6)[1], "dark gray"))+
  ylab("Density")+
  xlab("P(Invasion)")+
  scale_x_continuous(limits = c(0, 0.8))+
  theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8), legend.title = element_blank(),
        legend.position=c(0.6,0.85))
p_est_obs # error because you truncate x-axis, thus removing points

# final_cast$vals_in_0.3to0.4 <- ifelse(final_cast$predicted_ensemble > 0.3 & final_cast$predicted_ensemble < 0.4, "Yes", "No")
# table(final_cast$vals_in_0.3to0.4)
# bimodal_why_map <- ggplot() +
#   geom_sf(data = final_cast, aes(fill = vals_in_0.3to0.4)) + #polygons filled based on the density value
#   theme_bw()+theme_void()+
#   geom_sf(data=statesUS_albers_crop, fill="transparent", color="black", lwd=1)+
#   scale_fill_manual(values= c(adj_col,iso_col),
#                     name="", na.translate=FALSE)+
#   theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8), legend.justification = "left",
#         legend.position=c(0.3,0.25))

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# ROC plot construction --------------------------------------------
#
#
#---
library(pROC)
forcast_2019 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$year_predict == 2019),c("predicted_vals_all","Invaded_2019")] 
colnames(forcast_2019) <- c("predictions","status")
table(forcast_2019$status)
forcast_2020 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$year_predict == 2020),c("predicted_vals_all","Invaded_2020")] 
colnames(forcast_2020) <- c("predictions","status")
table(forcast_2020$status)
forcast_2021 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$year_predict == 2021),c("predicted_vals_all","Invaded_2021")] 
colnames(forcast_2021) <- c("predictions","status")
table(forcast_2021$status)

# combine datasets, set them up for ROC plot
forcast_fin <- rbind(forcast_2019,forcast_2019,forcast_2019)
rocobj <- roc(forcast_fin$status, forcast_fin$predictions)
auc <- round(auc(forcast_fin$status, forcast_fin$predictions),2)

#create ROC plot
roc_plot <- ggroc(rocobj, colour = viridis(6)[1], size = 2) +
  xlab("False-positive rate")+ylab("True-positive rate")+
  theme_bw() +theme_comps+
  annotate("text", x=0.5,y= 0.9, label= paste0('AUC = ', auc),
           size=4, hjust=0, vjust=1, fontface =1)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




# Figure: ROC and Density plot --------------------------------------------
#
#
#---
resize.win(6.85/2,5)
ggarrange(roc_plot,p_est_obs, ncol = 1, nrow = 2,  align = "hv",
          labels="auto", label.x =0.01, label.y =1)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# Forecast maps (prediction data sets) --------------------------------------------
#
#
#---
#'''''''''''''''''''''''''
# 2019
#'''''''''''''''''''''''''
prediction2019 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$year_predict == 2019),c("FIPS", "predicted_vals_all")]
nrow(prediction2019)
final_cast2019 <- merge(spread_ALBERS_hosts,prediction2019,by="FIPS", all.x=T)
summary(final_cast2019)
nrow(final_cast2019)


#'''''''''''''''''''''''''
# 2020
#'''''''''''''''''''''''''
prediction2020 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$year_predict == 2020),c("FIPS", "predicted_vals_all")]
nrow(prediction2020)
final_cast2020 <- merge(spread_ALBERS_hosts,prediction2020,by="FIPS", all.x=T)
summary(final_cast2020)
nrow(final_cast2020)


#'''''''''''''''''''''''''
# 2021
#'''''''''''''''''''''''''
prediction2021 <- all_hosts_buffer_forecast[which(all_hosts_buffer_forecast$year_predict == 2021),c("FIPS", "predicted_vals_all")]
nrow(prediction2021)
final_cast2021 <- merge(spread_ALBERS_hosts,prediction2021,by="FIPS", all.x=T)
summary(final_cast2021)
nrow(final_cast2021)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




# truncating the forecasts (removing really low estimates) --------------------------------------------
#
#
#---
final_cast2019$predicted_trunc <- ifelse(final_cast2019$predicted_vals_all < 0.01, NA, final_cast2019$predicted_vals_all)
table(is.na(final_cast2019$predicted_trunc))
#
final_cast2020$predicted_trunc <- ifelse(final_cast2020$predicted_vals_all < 0.01, NA, final_cast2020$predicted_vals_all)
table(is.na(final_cast2020$predicted_trunc))
#
final_cast2021$predicted_trunc <- ifelse(final_cast2021$predicted_vals_all < 0.01, NA, final_cast2021$predicted_vals_all)
table(is.na(final_cast2021$predicted_trunc))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# Forecast maps (summary statistics) --------------------------------------------
#
#
#---
# aqui
#'''''''''''''''''''''''''
# 2019
#'''''''''''''''''''''''''
final_cast_summary_I2019 <- final_cast2019[which(final_cast2019$YrInv %in% 2019),]
median(final_cast_summary_I2019$predicted_vals_all)
#
final_cast_summary_NI2019 <- final_cast2019[is.na(final_cast2019$YrInv),]
median(final_cast_summary_NI2019$predicted_vals_all)
counties_for_forecast_2019 <- rbind(final_cast_summary_NI2019,final_cast_summary_I2019)
#
counties_for_forecast_0.01_2019 <-  final_cast2019[which(final_cast2019$predicted_vals_all >= 0.01),]
nrow(counties_for_forecast_0.01_2019)/nrow(spread_ALBERS)
(nrow(spread_ALBERS)- nrow(counties_for_forecast_0.01_2019))
(nrow(spread_ALBERS)- nrow(counties_for_forecast_0.01_2019))/nrow(spread_ALBERS)



#'''''''''''''''''''''''''
# 2020
#'''''''''''''''''''''''''
final_cast_summary_I2020 <- final_cast2020[which(final_cast2020$YrInv %in% 2020),]
median(final_cast_summary_I2020$predicted_vals_all)
#
final_cast_summary_NI2020 <- final_cast2020[is.na(final_cast2020$YrInv),]
median(final_cast_summary_NI2020$predicted_vals_all)
counties_for_forecast_2020 <- rbind(final_cast_summary_NI2020,final_cast_summary_I2020)
#
counties_for_forecast_0.01_2020 <-  final_cast2020[which(final_cast2020$predicted_vals_all >=0.01),]
nrow(counties_for_forecast_0.01_2020)/nrow(spread_ALBERS)
(nrow(spread_ALBERS)- nrow(counties_for_forecast_0.01_2020))
(nrow(spread_ALBERS)- nrow(counties_for_forecast_0.01_2020))/nrow(spread_ALBERS)

#'''''''''''''''''''''''''
# 2021
#'''''''''''''''''''''''''
final_cast_summary_I2021 <- final_cast2021[which(final_cast2021$YrInv %in% 2021),]
median(final_cast_summary_I2021$predicted_vals_all)
#
final_cast_summary_NI2021 <- final_cast2021[is.na(final_cast2021$YrInv),]
median(final_cast_summary_NI2021$predicted_vals_all)
counties_for_forecast_2021 <- rbind(final_cast_summary_NI2021,final_cast_summary_I2021)
#
counties_for_forecast_0.01_2021 <-  final_cast2021[which(final_cast2021$predicted_vals_all >=0.01),]
nrow(counties_for_forecast_0.01_2021)/nrow(spread_ALBERS)
(nrow(spread_ALBERS)- nrow(counties_for_forecast_0.01_2021))
(nrow(spread_ALBERS)- nrow(counties_for_forecast_0.01_2021))/nrow(spread_ALBERS)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# Figure: Forecast maps --------------------------------------------
#
#
#---
my_breaks_forecast <- c(0.01, 0.2, 0.4,0.6, 0.8)
summary(final_cast2019$predicted_vals_all)
summary(final_cast2020$predicted_vals_all)
summary(final_cast2021$predicted_vals_all)

cast_p_2019 <- ggplot() +
  geom_sf(data = final_cast2019, aes(fill = predicted_trunc), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=final_cast2019[which(final_cast2019$YrInv < 2019), ], fill=alpha(viridis(24)[24],0.2), color="transparent", lwd=1)+
  scale_fill_continuous(low="gray", high=viridis(12)[1], 
                        na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                        breaks=my_breaks_forecast, labels=my_breaks_forecast,limits=c(0, 0.8))+
  geom_sf(data=statesUS_albers_hosts, fill="transparent", color="light gray", lwd=1)+
  theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left", 
        legend.title = element_blank(), legend.position=c(0.86,0.38))+
  geom_sf(data=final_cast2019[which(final_cast2019$YrInv==2019),], fill="transparent", color="black", lwd=0.01)+
  coord_sf(datum=st_crs(spread_ALBERS_hosts))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

cast_p_2020 <- ggplot() +
  geom_sf(data = final_cast2020, aes(fill = predicted_trunc), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=final_cast2020[which(final_cast2020$YrInv < 2020), ], fill=alpha(viridis(24)[24],0.2), color="transparent", lwd=1)+
  scale_fill_continuous(low="gray", high=viridis(12)[1], 
                        na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                        breaks=my_breaks_forecast, labels=my_breaks_forecast, limits=c(0, 0.8))+
  geom_sf(data=statesUS_albers_hosts, fill="transparent", color="light gray", lwd=1)+
  theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left", 
        legend.title = element_blank(), legend.position=c(0.86,0.38))+
  geom_sf(data=final_cast2020[which(final_cast2020$YrInv==2020),], fill="transparent", color="black", lwd=0.01)+
  coord_sf(datum=st_crs(spread_ALBERS_hosts))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

cast_p_2021 <- ggplot() +
  geom_sf(data = final_cast2021, aes(fill = predicted_trunc), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=final_cast2021[which(final_cast2021$YrInv < 2021), ], fill=alpha(viridis(24)[24],0.2), color="transparent", lwd=1)+
  scale_fill_continuous(low="gray", high=viridis(12)[1], 
                        na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                        breaks=my_breaks_forecast, labels=my_breaks_forecast, limits=c(0, 0.8))+
  geom_sf(data=statesUS_albers_hosts, fill="transparent", color="light gray", lwd=1)+
  theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left", 
        legend.title = element_blank(), legend.position=c(0.86,0.38))+
  geom_sf(data=final_cast2021[which(final_cast2021$YrInv==2021),], fill="transparent", color="black", lwd=0.01)+
  coord_sf(datum=st_crs(spread_ALBERS_hosts))+
  theme(plot.margin=grid::unit(c(-5,0,0,0), "mm"))



resize.win(6.85,9)
fig7 <- ggarrange(cast_p_2019, cast_p_2020, cast_p_2021, ncol = 1, nrow = 3,  align = "hv",
            label.x =0.24, label.y =0.93, labels=c(2019,2020,2021));fig7
ggsave("manuscript_spatial_drivers/figures/Fig7.tiff", plot=fig7, width = 6.85,  height = 9,  units = "in")
new7 <- image_read("manuscript_spatial_drivers/figures/Fig7.tiff")
new7_1 <- image_crop(new7, "1060x0+470+0")
image_write(new7_1, "manuscript_spatial_drivers/figures/Fig7_1.tiff")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

































# ANISOTROPY - add wedge predictor --------------------------------------------
#
#
#---
wedge_membership <- st_read("shapefiles/wedges/wedge_ID_sp.shp")
CPH_df$Wedge <- NA
cty_fips <- unique(CPH_df$FIPS)
i <- cty_fips[1]

for(i in cty_fips){
  curr_wedge <- wedge_membership[which(wedge_membership$FIPS %in% i),]
  CPH_df[which(CPH_df$FIPS %in% i), "Wedge"] <- curr_wedge$FID
}
CPH_df$Bearing <- as.factor(CPH_df$Wedge)
summary(CPH_df)
table(CPH_df$Bearing)


# ANISOTROPY - Full model no wedge --------------------------------------------
#
#
#---
fit_lwd_full_all_ANISO <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                              Redbay + Sassafras + Nonhosts + MinTemp + Precip,
                            data=CPH_df)
summary(fit_lwd_full_all_ANISO)
library(car)
Anova(fit_lwd_full_all_ANISO, type="III")
AIC(fit_lwd_full_all_ANISO)

options(na.action = "na.fail")
dd_all_ANISO <- dredge(fit_lwd_full_all_ANISO)
subset(dd_all_ANISO, delta < 2)
write.excel(subset(dd_all_ANISO, delta < 2))

model_averaged_all_ANISO<- model.avg(dd_all_ANISO)
write.excel(cbind(row.names(confint(model_averaged_all_ANISO)), confint(model_averaged_all_ANISO)))


# reduced 
fit_lwd_red_all_ANISO <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds +
                                 Redbay + Sassafras + MinTemp,
                               data=CPH_df_train)
summary(fit_lwd_red_all_ANISO)
AIC(fit_lwd_red_all_ANISO)

write.excel(cbind(rownames(summary(fit_lwd_red_all_ANISO)$coef),round(summary(fit_lwd_red_all_ANISO)$coef,4)))
write.excel(cbind(rownames(summary(fit_lwd_red_all_ANISO)$coef),round(summary(fit_lwd_red_all_ANISO)$coef,2)))


# https://stackoverflow.com/questions/54962119/how-to-plot-from-mumin-model-avg-summary
options(na.action = "na.fail") # needed for dredge to work
mA_all_ANISO <-summary(model_averaged_all_ANISO) #pulling out model averages
df1_all_ANISO<-as.data.frame(mA_all_ANISO$coefmat.full) #selecting full model coefficient averages
CI_all_ANISO <- as.data.frame(confint(model_averaged_all_ANISO, full=T)) # get confidence intervals for full model
df1_all_ANISO$CI.min <-CI_all_ANISO$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1_all_ANISO$CI.max <-CI_all_ANISO$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1_all_ANISO, keep.rownames = "coefficient") #put rownames into column
names(df1_all_ANISO) <- gsub(" ", "", names(df1_all_ANISO)) # remove spaces from column headers
#
CI_all_ANISO_red <- cbind(data.frame(row.names(confint(fit_lwd_red_all_ANISO))),
                          data.frame(coef(fit_lwd_red_all_ANISO)),
                          data.frame(confint(fit_lwd_red_all_ANISO))) # get confidence intervals for reduced model
colnames(CI_all_ANISO_red) <- c("coefficient", "Estimate", "CI.min" ,"CI.max")

every_variable <- row.names(confint(fit_lwd_full_all_ANISO)) 
missing_to_add <- every_variable[which(every_variable %!in% CI_all_ANISO_red$coefficient )]
missing_to_add <- data.frame(coefficient=missing_to_add, Estimate=NA, CI.min=NA, CI.max=NA)
CI_all_ANISO_red <- rbind.data.frame(CI_all_ANISO_red,missing_to_add)


df1_all_ANISO$type <- "Model-averaged"
CI_all_ANISO_red$type <- "Lowest AIC"
df_fin_all_ANISO <- rbind(df1_all_ANISO[,c("coefficient", "Estimate", "CI.min" ,"CI.max", "type")],
                          CI_all_ANISO_red)
#
df_fin_all_ANISO$coefficient <- factor(df_fin_all_ANISO$coefficient, levels=rev(c("Contagion", "Humans", "Campgrounds","Income",
                                                                                  "Redbay", "Sassafras", "Nonhosts", "MinTemp", "Precip")))

df_fin_all_ANISO$type <- factor(df_fin_all_ANISO$type, levels = c("Model-averaged", "Lowest AIC"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# ANISOTROPY Figure: coefficient graph NO wedge --------------------------------------------
#
#
#---
dodge <- position_dodge(width=0.6)  
base_size_text <- 12
point_size = 1.9

resize.win(6.85,4)
all_coef_no_wedge <- ggplot(data=df_fin_all_ANISO, aes(x=coefficient, y=Estimate, color=type))+ #again, excluding intercept because estimates so much larger
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1)+ #add dashed line at zero
  #
  geom_point(position=dodge, size=point_size) +
  geom_errorbar(aes(ymax=CI.max,ymin=CI.min),position = dodge, width=0, size=0) + 
  scale_colour_manual(values = c("gray", "black")) +
  guides(color = guide_legend(reverse = T))+
  #
  coord_flip()+ # flipping x and y axes
  theme_classic(base_size = base_size_text)+ xlab("")+ylab("")+
  theme(legend.position = c(0.8, 0.3), legend.title = element_blank())+
  scale_y_continuous(limits=c(-1, 6))
all_coef_no_wedge
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''








# ANISOTROPY - Full model WITH wedge --------------------------------------------
#
#
#---
fit_lwd_full_all_ANISO_wedge <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                                  Redbay + Sassafras + Nonhosts + MinTemp + Precip + Bearing,
                                data=CPH_df)
summary(fit_lwd_full_all_ANISO_wedge)
library(car)
Anova(fit_lwd_full_all_ANISO_wedge, type="III")
AIC(fit_lwd_full_all_ANISO_wedge)

options(na.action = "na.fail")
dd_all_ANISO_wedge <- dredge(fit_lwd_full_all_ANISO_wedge)
subset(dd_all_ANISO_wedge, delta < 2)
write.excel(subset(dd_all_ANISO_wedge, delta < 2))

model_averaged_all_ANISO_wedge <- model.avg(dd_all_ANISO_wedge, delta < 2) # AQUI 
write.excel(cbind(row.names(confint(model_averaged_all_ANISO_wedge)), confint(model_averaged_all_ANISO_wedge)))


# reduced 
fit_lwd_red_all_ANISO_wedge <- coxph(Surv(time0, time1, invaded_time1) ~ Contagion + Humans + Campgrounds + Income +
                                 Redbay + Sassafras + Bearing, data=CPH_df)
summary(fit_lwd_red_all_ANISO_wedge)
AIC(fit_lwd_red_all_ANISO_wedge)

write.excel(cbind(rownames(summary(fit_lwd_red_all_ANISO_wedge)$coef),round(summary(fit_lwd_red_all_ANISO_wedge)$coef,4)))
write.excel(cbind(rownames(summary(fit_lwd_red_all_ANISO_wedge)$coef),round(summary(fit_lwd_red_all_ANISO_wedge)$coef,2)))


# https://stackoverflow.com/questions/54962119/how-to-plot-from-mumin-model-avg-summary
options(na.action = "na.fail") # needed for dredge to work
mA_all_ANISO_wedge <-summary(model_averaged_all_ANISO_wedge) #pulling out model averages
df1_all_ANISO_wedge<-as.data.frame(mA_all_ANISO_wedge$coefmat.full) #selecting full model coefficient averages
CI_all_ANISO_wedge <- as.data.frame(confint(model_averaged_all_ANISO_wedge, full=T)) # get confidence intervals for full model
df1_all_ANISO_wedge$CI.min <-CI_all_ANISO_wedge$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1_all_ANISO_wedge$CI.max <-CI_all_ANISO_wedge$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1_all_ANISO_wedge, keep.rownames = "coefficient") #put rownames into column
names(df1_all_ANISO_wedge) <- gsub(" ", "", names(df1_all_ANISO_wedge)) # remove spaces from column headers
#
CI_all_ANISO_wedge_red <- cbind(data.frame(row.names(confint(fit_lwd_red_all_ANISO_wedge))),
                          data.frame(coef(fit_lwd_red_all_ANISO_wedge)),
                          data.frame(confint(fit_lwd_red_all_ANISO_wedge))) # get confidence intervals for reduced model
colnames(CI_all_ANISO_wedge_red) <- c("coefficient", "Estimate", "CI.min" ,"CI.max")

every_variable <- row.names(confint(fit_lwd_full_all_ANISO_wedge)) 
missing_to_add <- every_variable[which(every_variable %!in% CI_all_ANISO_wedge_red$coefficient)]
missing_to_add <- data.frame(coefficient=missing_to_add, Estimate=NA, CI.min=NA, CI.max=NA)
CI_all_ANISO_wedge_red <- rbind.data.frame(CI_all_ANISO_wedge_red,missing_to_add)


df1_all_ANISO_wedge$type <- "Model-averaged"
CI_all_ANISO_wedge_red$type <- "Lowest AIC"
df_fin_all_ANISO_wedge <- rbind(df1_all_ANISO_wedge[,c("coefficient", "Estimate", "CI.min" ,"CI.max", "type")],
                          CI_all_ANISO_wedge_red)
#
bearings <- c("22.5-45°","157.5-180°","180-202.5°","202.5-225°","225-247.5°","247.5-270°","270-292.5°","292.5-315°","315-337.5°","337.5-0°")
df_fin_all_ANISO_wedge[which(df_fin_all_ANISO_wedge$coefficient %in%  paste("Bearing",c(2,8:16), sep="")),"coefficient"] <- rep(bearings,2)


df_fin_all_ANISO_wedge$coefficient <- factor(df_fin_all_ANISO_wedge$coefficient, levels=rev(c("Contagion", "Humans", "Campgrounds","Income",
                                                                                  "Redbay", "Sassafras", "Nonhosts", "MinTemp", "Precip",
                                                                                  paste(bearings))))

df_fin_all_ANISO_wedge$type <- factor(df_fin_all_ANISO_wedge$type, levels = c("Model-averaged", "Lowest AIC"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# ANISOTROPY Figure: coefficient graph wedge--------------------------------------------
#
#
#---
dodge <- position_dodge(width=0.6)  
base_size_text <- 12
point_size = 1.9

resize.win(6.85,4)
all_coef_wedge <- ggplot(data=df_fin_all_ANISO_wedge, aes(x=coefficient, y=Estimate, color=type))+ #again, excluding intercept because estimates so much larger
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1)+ #add dashed line at zero
  #
  geom_point(position=dodge, size=point_size) +
  geom_errorbar(aes(ymax=CI.max,ymin=CI.min),position = dodge, width=0, size=0) + 
  scale_colour_manual(values = c("gray", "black")) +
  guides(color = guide_legend(reverse = T))+
  #
  coord_flip()+ # flipping x and y axes
  theme_classic(base_size = base_size_text)+ ylab("Slope coefficient ± 95% CL")+ xlab("")+
  theme(legend.position = c(0.8, 0.15), legend.title = element_blank())+
  scale_y_continuous(limits=c(-1, 6))
all_coef_wedge
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# ANISOTROPY combined figure --------------------------------------------
#
#
#---
resize.win(6.85,7)
ggarrange(all_coef_no_wedge,all_coef_wedge, ncol = 1, nrow = 2,  align = "hv",
          labels="auto", label.x =0.03, label.y =1, heights=c(1,2))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


