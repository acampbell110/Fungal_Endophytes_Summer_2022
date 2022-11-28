##############################################################################################
#### This script is used to load in the data and run analyses on the first set of surveys ####
##############################################################################################
setwd("/Users/alicampbell/Documents/GitHub/Fungal_Endophytes_Summer_2022/Summer Endophyte Surveys")
## Load necessary packages
library(readxl)
library(`nnet`)
library(lme4)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
## set up important function
invlogit <- function(x){exp(x)/(1+exp(x))}
## Read data in
data_2022 <- read_excel("Huntsville ELVI Survey.xlsx", sheet = 3)
data_2021 <- read_excel("Summer 2021 Seed Surveys.xlsx", sheet = 1)
data_2021$Site <- as.factor(data_2021$Site)
head(data_2021)
## Create a separate subset for each site
hunt_2021 <- subset(data_2021, data_2021$Site == "HUNT 21-1")
wall_2021 <- subset(data_2021, data_2021$Site == "WALL 21-3")
faye_2021 <- subset(data_2021, data_2021$Site == "FAYE21-2")
##################### Create All Models to Select From #####################################
###### Null ######
mod1 <- glm(data_2021$Endo_L~1, family = "binomial")
summary(mod1) ## p value is significant
###### VWC ######
mod2 <- glm(hunt_2021$`peel lib.`~hunt_2021$soil_moist, family = "binomial")
mod3 <- glm(wall_2021$Endo_L~wall_2021$vwc, family = "binomial")
mod4 <- glm(faye_2021$Endo_L~faye_2021$vwc, family = "binomial")
summary(mod2)
summary(mod3)
summary(mod4)
###### pH ######
mod5 <- glm(data_2021$Endo_L~data_2021$pH, family = "binomial")
summary(mod5)
###### Cover ######
mod6 <- glm(data_2021$Endo_L~data_2021$Cover, family = "binomial")
summary(mod6)
###### Site ######
mod7 <- glm(data_2021$Endo_L~data_2021$Site, family = "binomial")
summary(mod7)
###### VWC pH ######
mod8 <- glm(hunt_2021$Endo_L~hunt_2021$vwc + hunt_2021$pH, family = "binomial")
mod9 <- glm(hunt_2021$Endo_L~hunt_2021$vwc * hunt_2021$pH, family = "binomial")
mod10 <- glm(wall_2021$Endo_L~wall_2021$vwc + wall_2021$pH, family = "binomial")
mod11 <- glm(wall_2021$Endo_L~wall_2021$vwc * wall_2021$pH, family = "binomial")
mod12 <- glm(faye_2021$Endo_L~faye_2021$vwc + faye_2021$pH, family = "binomial")
mod13 <- glm(faye_2021$Endo_L~faye_2021$vwc * faye_2021$pH, family = "binomial")
summary(mod8)
summary(mod9)
summary(mod10)
summary(mod11) ## significant
summary(mod12)
summary(mod13)
###### VWC Cover ######
mod14 <- glm(hunt_2021$Endo_L~hunt_2021$vwc + hunt_2021$Cover, family = "binomial")
mod15 <- glm(hunt_2021$Endo_L~hunt_2021$vwc * hunt_2021$Cover, family = "binomial")
mod16 <- glm(wall_2021$Endo_L~wall_2021$vwc + wall_2021$Cover, family = "binomial")
mod17 <- glm(wall_2021$Endo_L~wall_2021$vwc * wall_2021$Cover, family = "binomial")
mod18 <- glm(faye_2021$Endo_L~faye_2021$vwc + faye_2021$Cover, family = "binomial")
mod19 <- glm(faye_2021$Endo_L~faye_2021$vwc * faye_2021$Cover, family = "binomial")
summary(mod14)
summary(mod15)
summary(mod16)
summary(mod17) ## significant
summary(mod18)
summary(mod19)
###### pH Cover ######
mod20 <- glm(data_2021$Endo_L~data_2021$pH + data_2021$Cover, family = "binomial")
mod21 <- glm(data_2021$Endo_L~data_2021$pH * data_2021$Cover, family = "binomial")
summary(mod20)
summary(mod21)
###### VWC pH Cover ######
mod22 <- glm(hunt_2021$Endo_L~hunt_2021$pH + hunt_2021$Cover + hunt_2021$vwc, family = "binomial")
mod23 <- glm(hunt_2021$Endo_L~hunt_2021$pH * hunt_2021$Cover * hunt_2021$vwc, family = "binomial")
mod24 <- glm(wall_2021$Endo_L~wall_2021$pH + wall_2021$Cover + wall_2021$vwc, family = "binomial")
mod25 <- glm(wall_2021$Endo_L~wall_2021$pH * wall_2021$Cover * wall_2021$vwc, family = "binomial")
mod26 <- glm(faye_2021$Endo_L~faye_2021$pH + faye_2021$Cover + faye_2021$vwc, family = "binomial")
mod27 <- glm(faye_2021$Endo_L~faye_2021$pH * faye_2021$Cover * faye_2021$vwc, family = "binomial")
summary(mod22)
summary(mod23)
summary(mod24)
summary(mod25)
summary(mod26)
summary(mod27)
###### Select the model that best addresses the varition in the endoohyte status)
#### No specific site
aictab(list(mod1,mod5,mod6,mod7,mod20,mod21))
## Here the best model appears to be the mod7 which includes only the site to deteermine the prevalence
#### Hunt
aictab(list(mod2,mod8,mod9,mod14,mod15,mod22,mod23))
## Here the best model appears to be mod2 which includes only vwc
#### Wall
aictab(list(mod3,mod10,mod11,mod16,mod17,mod24,mod25))
## Here the best model appears to be mod11 which includes vwc and pH as an interaction
#### Faye
aictab(list(mod4,mod12,mod13,mod18,mod19,mod26,mod27))
## Here the best model appears to be mod4 which includes vwc
################################# ANALYSIS ####################################################
## In all models, the inclusion of vwc appears to be the most important. The soil water content
## appears to be the best predictor of endophyte status. 
## Now I will test if there is a significant difference in the mean vwc of each E status
######## Huntsville 2022
## 
a <- subset(data_2022, data_2022$`peel lib.` == 1)
b <- subset(data_2022, data_2022$`peel lib.` == 0)
t.test(a$soil_moist,b$soil_moist, alternative = "less")
######## Hunstville 2021
##
a <- subset(hunt_2021, hunt_2021$Agrin_Endo == 1)
b <- subset(hunt_2021, hunt_2021$Agrin_Endo == 0)
t.test(a$vwc,b$vwc, alternative = "less")
######## Waller 2021
a <- subset(wall_2021, wall_2021$Agrin_Endo == 1)
b <- subset(wall_2021, wall_2021$Agrin_Endo == 0)
t.test(a$vwc,b$vwc, alternative = "less")
######## Faye 2021
a <- subset(faye_2021, faye_2021$Agrin_Endo == 1)
b <- subset(faye_2021, faye_2021$Agrin_Endo == 0)
t.test(a$vwc,b$vwc, alternative = "less")
##################### Visualize All Best Models #################################################
## Create colors
huntcol <- "deeppink2"
wallcol <- "darkorchid2"
fayecol <- "aquamarine2"
anycol <- "deepskyblue2"
## ggplot version (needs to be remade because the error bars go beyond 0)
a <- ggplot(data = hunt_2021,aes(x = vwc, y = Agrin_Endo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = huntcol,level = 0.5) + 
  theme_classic() + 
  labs(x = "") +   labs(y = "") + 
  ylim(0,1) 
d <- ggplot(data = data_2022,aes(x = soil_moist, y = `peel lib.`)) + 
  geom_point() +
  stat_smooth(method = "lm", col = anycol,level = 0.5) + 
  theme_classic() + 
  labs(x = "") +   labs(y = "") + 
  ylim(0,1) 
b <- ggplot(data = wall_2021,aes(x = vwc, y = Agrin_Endo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = wallcol,level = 0.5) + 
  labs(x = "") +   labs(y = "") + 
  theme_classic() + 
  ylim(0,1) 
c <- ggplot(data = faye_2021,aes(x = vwc, y = Agrin_Endo)) + 
  geom_point() +
  stat_smooth(method = "lm", col = fayecol,level = 0.5) + 
  labs(x = "",y = "")  + 
  theme_classic() + 
  ylim(0,1) 
e <- ggplot() + geom_text(aes(x=0, y=0, label = "Endophyte Status"), 
               parse = TRUE, size = 6, angle = 90, vjust = 1.4) +
  theme_void()
f <- ggplot() + geom_text(aes(x=0, y=0, label = " "), 
               parse = TRUE, size = 6, hjust = -1) +
  theme_void()
g <- ggplot() + geom_text(aes(x=0, y=0, label = "Soil"), 
                            parse = TRUE, size = 6, vjust = -.8) +
  theme_void()
group <- ggarrange(a,b,c,d,
                   labels = c("a)","b)","c)","d)"),
                   ncol = 4, nrow = 1)
png("Soil_Water_Content2.png")
annotate_figure(group,left = textGrob("Endophyte Status",rot = 90, vjust = 1.4, gp = gpar(cex = 1.5)),
                bottom = textGrob("Soil Water Content", vjust = -.2, gp = gpar(cex = 1.5)))
dev.off()

png("Other Envi Variables.png")
a <- ggplot(data = data_2021, aes(x = Cover, y = Endo_L)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = anycol, level = 0.5) + 
  theme_classic()
b <- ggplot(data = data_2021, aes(x = pH, y = Endo_L)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = anycol, level = 0.5) + 
  theme_classic()
c <- ggplot(data = data_2021, aes(x = Cover, y = Endo_L)) + 
  geom_point() + 
  stat_smooth(method = "lm", col = anycol, level = 0.5) + 
  theme_classic()
ggarrange(a,b,c,
          lebels = c("a)","b)","c)"),
          ncol = 3, nrow = 1)
dev.off()

##############################
## FAKE FIG NEEDED

dummy_1 <- seq(1,10, by = .1)
y_1 <- exp(-dummy_1)
y_2 <- exp(-1.3*dummy_1)
y_3 <- exp(-1.6*dummy_1)
y_4 = 8*exp(-.4*dummy_1)
y_5 = -0.1*exp(.4*dummy_1) + 5
y_6 = sin(dummy_1)
y_7 = sin(dummy_1 + 2)
y_8 = dummy_1

png("Coexistence_Mech.png")
par(mar=c(3,3,2,2))
layout(matrix(c(1,2),
              ncol = 2, nrow = 1, byrow = T), heights = c(1)) 
## Storage Effect
plot(dummy_1,y_1,ylim = c(0,.3), type = "l", lwd = 3, col = "red", xlim = c(1,10), 
     main = "                                                                      ",
     xlab = "", ylab = "")
lines(dummy_1,y_2, lwd = 3, col = "orange")
lines(dummy_1,y_3, lwd = 3, col = "gold")
legend("topright",legend = c("strong E","med. E","weak E"), fill = c("red","orange","gold"))
title(ylab="Fitness", line=1.8, cex.lab=1.2)
title(xlab = "Effects of Competition", line = 1.8, cex.lab = 1.2)
title(main = "a)                                           ")
## Spatial Relative Non-linearity
plot(dummy_1,y_5, type = "l", ylim = c(0,6), col = "deeppink", lwd = 3,
     main = "b)                                                                    ",
     xlab = "", ylab = "")
lines(dummy_1, y_4, col = "darkorchid", lwd = 3)
legend("topright",legend = c("species 1","species 2"),fill = c("deeppink","darkorchid"))
title(ylab="Fitness", line=1.8, cex.lab=1.2)
title(xlab = "Effects of Competition", line = 1.8, cex.lab = 1.2)
title(main = "b)                                           ")                                         ")
dev.off()
