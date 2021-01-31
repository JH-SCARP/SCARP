#1-week risk analysis for progression to severe disease or death in patients hospitalized with COVID-19  
#Prospective Validation (use data from March 5 - July 4, 2020 for development and data from July 5 - December 4, 2020 for validation)

library(tidyverse)
library(magrittr)
library(dplyr)
library(pROC)
library(rfSLAM)
library(splines)
library(grid)
library(gridExtra)


# 1-week risk predictions
# train is the dataframe with data from March 5 - July 4, 2020 in the 6-hour interval format
node.size <- round((dim(train)[1])*0.1)
set.seed(321) 
start_time <- Sys.time()
mymodel.1.full <- rfSLAM::rfsrc(i.sev_died.7 ~., data = train, nodesize = node.size, ntree = 100, na.action = "na.impute", splitrule="poisson.split1", risk.time = train_b$rt_7p, stratifier=train_b$q6, membership = TRUE, bootstrap = "by.user", samp = rfslam.samp.1)
end_time <- Sys.time()

# obtain predicted event rates 
p.cpiu.be.ppl <- risk.adjust.be(rf = mymodel.1.full, status = train_b$i.sev_died.7, rt = train_b$rt_7p, alpha.tm = 0)

train_b$p.hat <- 1-exp(-p.cpiu.be.ppl) # predicted event probabilities 
train_b$ni.sca <- as.numeric(train_b$i.sev_died.7)-1


# calibrate predictions
lr_model = glm(as.factor(ni.sca) ~ ns(p.hat,2)*ns(q6,2),data = train_b, family = binomial)

# predict using test data
new.pred <- predict(mymodel.1.full, newdata = test, na.action = "na.impute", membership = TRUE)

new.p <- risk.adjust.new.be(rf=mymodel.1.full, status = train_b$i.sev_died.7, rt = train_b$rt_7p, pre.mem = new.pred, k = 2, alpha.tm = 0)


test$p.hat <- 1-exp(-new.p)

p.hat.df <- data.frame(p.hat = test$p.hat, q6 = test$q6)

# calibrated predictions
lr_probs = predict(lr_model,  
                   newdata = p.hat.df, 
                   type = "response")

test$p.cal <- lr_probs


##### Calibration ##### 

db2 <- test
db2$p.hat <- db2$p.cal

g1 <- mutate(db2, bin = ntile(p.hat, 10)) %>% 
  # Bin prediction into 10ths
  group_by(bin) %>%
  mutate(n = n(), # Get ests and CIs
         bin_pred = mean(p.hat), 
         bin_prob = mean(as.numeric(i.sev_died.7) - 1), 
         se = sqrt((bin_prob * (1 - bin_prob)) / n), 
         ul = bin_prob + 1.96 * se, 
         ll = bin_prob - 1.96 * se) %>%
  ungroup() %>%
  ggplot(aes(x = bin_pred, y = bin_prob, ymin = ll, ymax = ul)) +
  geom_pointrange(size = 0.5, color = "black") +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
  scale_x_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
  geom_abline() +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,face = "bold")) +
  ggtitle("Calibration")

g2 <- ggplot(db2, aes(x = p.hat)) +
  geom_histogram(fill = "black", bins = 200) +
  scale_x_continuous(limits = c(0, 0.7), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1200),breaks = seq(0, 1200, by = 300)) +
  theme(panel.grid.minor = element_blank())

g <- arrangeGrob(g1, g2, respect = TRUE, heights = c(1, 0.25), ncol = 1)
grid.newpage()
grid.draw(g)

#### AUC ####

auc.df <- data.frame(matrix(NA, nrow = 2, ncol = 4))

names(auc.df) <- c("time", "lci", "auc", "uci")

auc.df$time <- c(1, 2)

# week 1
t0 <- db2 %>% filter(q6 <= 27)
auc.df[1,2:4] <- unname(ci.auc(t0$i.sev_died.7, t0$p.hat))


# week 2
t0 <- db2 %>% filter(q6 > 27)
auc.df[2,2:4] <- unname(ci.auc(t0$i.sev_died.7, t0$p.hat))
