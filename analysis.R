## Setup ----
# Load libraries
library(dplyr) #for manipulating data
library(tidyverse) #for data manipulating
library(car) #for stats e.g. Levene test not in the core stats package
library(vegan) #for species diversity indices, NMDS, and computing permanova and anosim
library(betapart) #for nestedness and turnover of communities
library(ggplot2) #for plotting graphs
library(ggpubr) #for plotting boxplots
library(multcompView) #for significant difference letters
library(here) # for wasy file referencing
library(svglite) #for saving figures as .svg
# renv, already installed in RStudio, used for version control 
all_data <- readr::read_csv(
  here::here("Data", "2025-03-19_all-variables-measured.csv"), show_col_types = FALSE
) 
####  Alpha diversity of vegetation ----

#boxplot of grassland vs heathland, species richness
richness_bxp <- ggboxplot(all_data, x = "Habitat", y = "Vegetation Richness", color = "Vegetation", ylab = "Species Richness", palette = c("black", "limegreen"), lwd = 0.75) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
  #remove legend title
  legend.title = element_blank()
) 
#graphing boxplots with bracken split from nonbracken
shannon_bxp <-ggboxplot(all_data, x = "Habitat", y = "Vegetation Shannon", color = "Vegetation", ylab = "Shannon Diversity", palette = c("black", "limegreen"), lwd = 0.75) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
  #remove legend title
  legend.title = element_blank()
) 
#graphing boxplots with bracken split from nonbracken
simpson_bxp <- ggboxplot(all_data, x = "Habitat", y = "Vegetation Simpson", color = "Vegetation", ylab = "Simpson Diversity", palette = c("black", "limegreen"), lwd = 0.75) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
  #remove legend title
  legend.title = element_blank()
) 
#boxplot of grassland vs heathland, species richness
evenness_bxp <- ggboxplot(all_data, x = "Habitat", y = "Vegetation Pielou's Evenness", color = "Vegetation", ylab = "Pielou's Evenness", palette = c("black", "limegreen"), lwd = 0.75) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
  #remove legend title
  legend.title = element_blank()
) 

#three sites have only 1 species, meaning the eveness is NaN.  Do we make this 0, or omit as the code is currently doing?
all_bxp <- ggarrange(richness_bxp, evenness_bxp, shannon_bxp, simpson_bxp, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     common.legend = TRUE)
#show the plot in the Plots window
show(all_bxp)

#save the image
ggsave(path = "Figures", paste0(Sys.Date(), "_4-panel-veg-diversity.svg"), all_bxp, width = 7, height = 5, dpi = 300)


#Vegetation Richness data is normal.

#run the GLM
glm <- glm(all_data$`Vegetation Richness` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)


nested_anova <- aov(all_data$`Vegetation Richness` ~ all_data$Habitat / factor(all_data$Vegetation))
summary(nested_anova)
#check homogeneity of variance
plot(nested_anova, 1)
#levene test.  if p value < 0.05, there is eidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Vegetation Richness` ~ all_data$Habitat / factor(all_data$Vegetation))
#check normality.  
plot(nested_anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
nestaov_residuals <- residuals(object = nested_anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = nestaov_residuals)


#do anova on 
nested_anova <- aov(d$evenness ~ d$Habitat / factor(d$Vegetation))
summary(nested_anova)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(nested_anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)
