#1-day risk analysis for progression to severe disease or death in patients hospitalized with COVID-19  
#Hospital-Site Cross-Validation (with the main hospital site 1 always included)
#Hospital sites 2 - 5 are left out one at a time, sequentially 

library(tidyverse)
library(magrittr)
library(dplyr)
library(pROC)
library(rfSLAM)
library(splines)
library(grid)
library(gridExtra)

# rf.cpiu.1 is the dataframe with data in the 6-hour interval format

############ hospital 2 ###############

i = 2

  testIndexes <- which(rf.cpiu.1$admit_hospital_num==i,arr.ind=TRUE)
  testData <- rf.cpiu.1[rf.cpiu.1$osler_id %in% ids[testIndexes], ]
  trainData <- rf.cpiu.1[rf.cpiu.1$osler_id %nin% ids[testIndexes], ]
  
  id.1 <- unique(trainData$pid)
  rfslam.samp.1 <- boot.samp(ntree = 100, id = id.1, boot.df.rep = trainData) 
  rf.df.1 <- trainData[,sel_var]
  
  rf.df.1[, sapply(rf.df.1, class) == 'character'] <- lapply(rf.df.1[, sapply(rf.df.1, class) == 'character'], factor)
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  rf.df.1$i.sev_died.1 <- trainData$i.sev_died.1
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("sev_died_within_1_days")]
  
  set.seed(321)
  node.size <- round((dim(trainData)[1])*0.1)
  mymodel.1 <- rfSLAM::rfsrc(i.sev_died.1 ~., data = rf.df.1, nodesize = node.size, ntree = 100, na.action = "na.impute", splitrule="poisson.split1", risk.time = trainData$rt_1p, stratifier=rf.df.1$q6, membership = TRUE, bootstrap = "by.user", samp = rfslam.samp.1)
  
  ## obtain predicted event rates 
  
  p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, alpha.tm = 0)
  
  ## Calibrate predictions
  
  rf.df.1$p.hat <- 1-exp(-p.cpiu.be.ppl) # predicted event probabilities 
  rf.df.1$ni.sca <- as.numeric(rf.df.1$i.sev_died.1)-1
  
  lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),  
                 data = rf.df.1, family = binomial)
  
  ## new predictions
  
  test <- testData[,sel_var]
  
  test[, sapply(test, class) == 'character'] <- lapply(test[, sapply(test, class) == 'character'], factor)
  test <- test[,!names(test) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  test <- test[,!names(test) %in% c("sev_died_within_1_days")]
  
  new.pred <- predict(mymodel.1, newdata = test, na.action = "na.impute", membership = TRUE)
  
  new.p <- risk.adjust.new.be(rf=mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, pre.mem = new.pred, k = 2, alpha.tm = 0)
  
  test.cpiu <- test
  test.cpiu$p.hat <- 1-exp(-new.p)
  
  p.hat.df <- data.frame(p.hat = test.cpiu$p.hat, q6 = test.cpiu$q6)
  
  lr_probs = predict(lr_model,  
                     newdata = p.hat.df, 
                     type = "response")
  
  testData$p.hat <- lr_probs
  
  
  phat.list[[i]] <- testData

  rm(mymodel.1)

############ hospital 3 ###############

i = 3
  testIndexes <- which(rf.cpiu.1$admit_hospital_num==i,arr.ind=TRUE)
  testData <- rf.cpiu.1[rf.cpiu.1$osler_id %in% ids[testIndexes], ]
  trainData <- rf.cpiu.1[rf.cpiu.1$osler_id %nin% ids[testIndexes], ]
  
  id.1 <- unique(trainData$pid)
  rfslam.samp.1 <- boot.samp(ntree = 100, id = id.1, boot.df.rep = trainData) 
  rf.df.1 <- trainData[,sel_var]
  
  rf.df.1[, sapply(rf.df.1, class) == 'character'] <- lapply(rf.df.1[, sapply(rf.df.1, class) == 'character'], factor)
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  rf.df.1$i.sev_died.1 <- trainData$i.sev_died.1
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("sev_died_within_1_days")]
  
  set.seed(321)
  node.size <- round((dim(trainData)[1])*0.1)
  mymodel.1 <- rfSLAM::rfsrc(i.sev_died.1 ~., data = rf.df.1, nodesize = node.size, ntree = 100, na.action = "na.impute", splitrule="poisson.split1", risk.time = trainData$rt_1p, stratifier=rf.df.1$q6, membership = TRUE, bootstrap = "by.user", samp = rfslam.samp.1)
  
  ## obtain predicted event rates 
  
  p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, alpha.tm = 0)
  
  ## Calibrate predictions
  
  rf.df.1$p.hat <- 1-exp(-p.cpiu.be.ppl) # predicted event probabilities 
  rf.df.1$ni.sca <- as.numeric(rf.df.1$i.sev_died.1)-1
  
  lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),  
                 data = rf.df.1, family = binomial)
  
  ## new predictions
  
  test <- testData[,sel_var]
  
  test[, sapply(test, class) == 'character'] <- lapply(test[, sapply(test, class) == 'character'], factor)
  test <- test[,!names(test) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  test <- test[,!names(test) %in% c("sev_died_within_1_days")]
  
  new.pred <- predict(mymodel.1, newdata = test, na.action = "na.impute", membership = TRUE)
  
  new.p <- risk.adjust.new.be(rf=mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, pre.mem = new.pred, k = 2, alpha.tm = 0)
  
  test.cpiu <- test
  test.cpiu$p.hat <- 1-exp(-new.p)
  
  p.hat.df <- data.frame(p.hat = test.cpiu$p.hat, q6 = test.cpiu$q6)
  
  lr_probs = predict(lr_model,  
                     newdata = p.hat.df, 
                     type = "response")
  
  testData$p.hat <- lr_probs
  
  
  phat.list[[i]] <- testData
  
  rm(mymodel.1)


############ hospital 4 ###############

i = 4
  
  testIndexes <- which(rf.cpiu.1$admit_hospital_num==i,arr.ind=TRUE)
  testData <- rf.cpiu.1[rf.cpiu.1$osler_id %in% ids[testIndexes], ]
  trainData <- rf.cpiu.1[rf.cpiu.1$osler_id %nin% ids[testIndexes], ]
  
  id.1 <- unique(trainData$pid)
  rfslam.samp.1 <- boot.samp(ntree = 100, id = id.1, boot.df.rep = trainData) 
  rf.df.1 <- trainData[,sel_var]
  
  rf.df.1[, sapply(rf.df.1, class) == 'character'] <- lapply(rf.df.1[, sapply(rf.df.1, class) == 'character'], factor)
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  rf.df.1$i.sev_died.1 <- trainData$i.sev_died.1
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("sev_died_within_1_days")]
  
  set.seed(321)
  node.size <- round((dim(trainData)[1])*0.1)
  mymodel.1 <- rfSLAM::rfsrc(i.sev_died.1 ~., data = rf.df.1, nodesize = node.size, ntree = 100, na.action = "na.impute", splitrule="poisson.split1", risk.time = trainData$rt_1p, stratifier=rf.df.1$q6, membership = TRUE, bootstrap = "by.user", samp = rfslam.samp.1)
  
  ## obtain predicted event rates 
  
  p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, alpha.tm = 0)
  
  ## calibrate predictions
  
  rf.df.1$p.hat <- 1-exp(-p.cpiu.be.ppl) # predicted event probabilities 
  rf.df.1$ni.sca <- as.numeric(rf.df.1$i.sev_died.1)-1
  
  lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),  
                 data = rf.df.1, family = binomial)
  
  ## new predictions
  
  test <- testData[,sel_var]
  
  test[, sapply(test, class) == 'character'] <- lapply(test[, sapply(test, class) == 'character'], factor)
  test <- test[,!names(test) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  test <- test[,!names(test) %in% c("sev_died_within_1_days")]
  
  new.pred <- predict(mymodel.1, newdata = test, na.action = "na.impute", membership = TRUE)
  
  new.p <- risk.adjust.new.be(rf=mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, pre.mem = new.pred, k = 2, alpha.tm = 0)
  
  test.cpiu <- test
  test.cpiu$p.hat <- 1-exp(-new.p)
  
  p.hat.df <- data.frame(p.hat = test.cpiu$p.hat, q6 = test.cpiu$q6)
  
  lr_probs = predict(lr_model,  
                     newdata = p.hat.df, 
                     type = "response")
  
  testData$p.hat <- lr_probs
  
  
  phat.list[[i]] <- testData
  
  rm(mymodel.1)


############ hospital 5 ###############

i = 5
  
  testIndexes <- which(rf.cpiu.1$admit_hospital_num==i,arr.ind=TRUE)
  testData <- rf.cpiu.1[rf.cpiu.1$osler_id %in% ids[testIndexes], ]
  trainData <- rf.cpiu.1[rf.cpiu.1$osler_id %nin% ids[testIndexes], ]
  
  id.1 <- unique(trainData$pid)
  rfslam.samp.1 <- boot.samp(ntree = 100, id = id.1, boot.df.rep = trainData) 
  rf.df.1 <- trainData[,sel_var]
  
  rf.df.1[, sapply(rf.df.1, class) == 'character'] <- lapply(rf.df.1[, sapply(rf.df.1, class) == 'character'], factor)
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  rf.df.1$i.sev_died.1 <- trainData$i.sev_died.1
  rf.df.1 <- rf.df.1[,!names(rf.df.1) %in% c("sev_died_within_1_days")]
  
  set.seed(321)
  node.size <- round((dim(trainData)[1])*0.1)
  mymodel.1 <- rfSLAM::rfsrc(i.sev_died.1 ~., data = rf.df.1, nodesize = node.size, ntree = 100, na.action = "na.impute", splitrule="poisson.split1", risk.time = trainData$rt_1p, stratifier=rf.df.1$q6, membership = TRUE, bootstrap = "by.user", samp = rfslam.samp.1)
  
  ## obtain predicted event rates 
  
  p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, alpha.tm = 0)
  
  ## calibrate predictions
  
  rf.df.1$p.hat <- 1-exp(-p.cpiu.be.ppl) # predicted event probabilities 
  rf.df.1$ni.sca <- as.numeric(rf.df.1$i.sev_died.1)-1
  
  lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),  
                 data = rf.df.1, family = binomial)
  
  ## New predictions
  
  test <- testData[,sel_var]
  
  test[, sapply(test, class) == 'character'] <- lapply(test[, sapply(test, class) == 'character'], factor)
  test <- test[,!names(test) %in% c("osler_id", "time_to_death", "time_to_severe_or_death", "time_to_last_known_status")]
  
  test <- test[,!names(test) %in% c("sev_died_within_1_days")]
  
  new.pred <- predict(mymodel.1, newdata = test, na.action = "na.impute", membership = TRUE)
  
  new.p <- risk.adjust.new.be(rf=mymodel.1, status = rf.df.1$i.sev_died.1, rt = trainData$rt_1p, pre.mem = new.pred, k = 2, alpha.tm = 0)
  
  test.cpiu <- test
  test.cpiu$p.hat <- 1-exp(-new.p)
  
  p.hat.df <- data.frame(p.hat = test.cpiu$p.hat, q6 = test.cpiu$q6)
  
  lr_probs = predict(lr_model,  
                     newdata = p.hat.df, 
                     type = "response")
  
  testData$p.hat <- lr_probs
  
  
  phat.list[[i]] <- testData
  
  rm(mymodel.1)

##### performance metrics

auc.df <- data.frame(matrix(NA, nrow = 8, ncol = 5))

names(auc.df) <- c("time", "lci", "auc", "uci", "hosp")

auc.df$time <- c(1, 2)
auc.df$hosp <- c(2, 2, 3, 3, 4, 4, 5, 5)

db2 <- phat.list[[2]]

# week 1
t0 <- db2 %>% filter(q6 <= 27)
auc.df[1,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


# week 2
t0 <- db2 %>% filter(q6 > 27)
auc.df[2,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


## calibration plot
library(viridis)
library(gridExtra)

g1 <- mutate(db2, bin = ntile(p.hat, 10)) %>% 
  # Bin prediction into 10ths
  group_by(bin) %>%
  mutate(n = n(), # Get ests and CIs
         bin_pred = mean(p.hat), 
         bin_prob = mean(as.numeric(i.sev_died.1) - 1), 
         se = sqrt((bin_prob * (1 - bin_prob)) / n), 
         ul = bin_prob + 1.96 * se, 
         ll = bin_prob - 1.96 * se) %>%
  ungroup() %>%
  ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  geom_pointrange(size = 0.5, color = "black") +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
  geom_abline() + geom_smooth(method = "loess", se = FALSE, linetype = "dashed", 
                                color = "black") +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggtitle("Calibration")

g2 <- ggplot(db2, aes(x = p.hat)) +
  geom_histogram(fill = "black", bins = 200) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous(limits = c(0,3000),breaks = seq(0, 3000, by = 1000)) +
  theme(panel.grid.minor = element_blank())

g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
grid.newpage()
grid.draw(g)



db2 <- phat.list[[3]]

# week 1
t0 <- db2 %>% filter(q6 <= 27)
auc.df[3,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


# week 2
t0 <- db2 %>% filter(q6 > 27)
auc.df[4,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


## calibration plot

g1 <- mutate(db2, bin = ntile(p.hat, 10)) %>% 
  # Bin prediction into 10ths
  group_by(bin) %>%
  mutate(n = n(), # Get ests and CIs
         bin_pred = mean(p.hat), 
         bin_prob = mean(as.numeric(i.sev_died.1) - 1), 
         se = sqrt((bin_prob * (1 - bin_prob)) / n), 
         ul = bin_prob + 1.96 * se, 
         ll = bin_prob - 1.96 * se) %>%
  ungroup() %>%
  ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  geom_pointrange(size = 0.5, color = "black") +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
  geom_abline() + geom_smooth(method = "loess", se = FALSE, linetype = "dashed", 
                                color = "black") +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggtitle("Calibration")

g2 <- ggplot(db2, aes(x = p.hat)) +
  geom_histogram(fill = "black", bins = 200) +
  scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous(limits = c(0,4000),breaks = seq(0, 4000, by = 1000)) +
  theme(panel.grid.minor = element_blank())

g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
grid.newpage()
grid.draw(g)


db2 <- phat.list[[4]]


# week 1
t0 <- db2 %>% filter(q6 <= 27)
auc.df[5,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


# week 2
t0 <- db2 %>% filter(q6 > 27)
auc.df[6,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


## calibration plot

g1 <- mutate(db2, bin = ntile(p.hat, 10)) %>% 
  # Bin prediction into 10ths
  group_by(bin) %>%
  mutate(n = n(), # Get ests and CIs
         bin_pred = mean(p.hat), 
         bin_prob = mean(as.numeric(i.sev_died.1) - 1), 
         se = sqrt((bin_prob * (1 - bin_prob)) / n), 
         ul = bin_prob + 1.96 * se, 
         ll = bin_prob - 1.96 * se) %>%
  ungroup() %>%
  ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  geom_pointrange(size = 0.5, color = "black") +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  geom_abline() + geom_smooth(method = "loess", se = FALSE, linetype = "dashed", 
                              color = "black", size = 1) +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggtitle("Calibration")

g2 <- ggplot(db2, aes(x = p.hat)) +
  geom_histogram(fill = "black", bins = 200) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous(limits = c(0,5000),breaks = seq(0, 5000, by = 1000)) +
  theme(panel.grid.minor = element_blank())

g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
grid.newpage()
grid.draw(g)



db2 <- phat.list[[5]]


# week 1
t0 <- db2 %>% filter(q6 <= 27)
auc.df[7,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


# week 2
t0 <- db2 %>% filter(q6 > 27)
auc.df[8,2:4] <- unname(ci.auc(t0$i.sev_died.1, t0$p.hat))


## calibration plot

g1 <- mutate(db2, bin = ntile(p.hat, 10)) %>% 
  # Bin prediction into 10ths
  group_by(bin) %>%
  mutate(n = n(), # Get ests and CIs
         bin_pred = mean(p.hat), 
         bin_prob = mean(as.numeric(i.sev_died.1) - 1), 
         se = sqrt((bin_prob * (1 - bin_prob)) / n), 
         ul = bin_prob + 1.96 * se, 
         ll = bin_prob - 1.96 * se) %>%
  ungroup() %>%
  ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  geom_pointrange(size = 0.5, color = "black") +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  geom_abline() + geom_smooth(method = "loess", se = FALSE, linetype = "dashed", 
                                color = "black", size = 1, span = 0.2) +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggtitle("Calibration")

g2 <- ggplot(db2, aes(x = p.hat)) +
  geom_histogram(fill = "black", bins = 200) +
  scale_x_continuous(limits = c(0, 0.6), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous(limits = c(0,2000),breaks = seq(0, 2000, by = 500)) +
  theme(panel.grid.minor = element_blank())

g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
grid.newpage()
grid.draw(g)