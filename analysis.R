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

#### load the master data file ----
all_data <- readr::read_csv(
  here::here("Data", "2025-03-19_all-variables-measured.csv"), show_col_types = FALSE
) 
all_data <- as.data.frame(all_data)
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

####  Beta diversity of vegetation ----

#extract the veg abundances
spe <- all_data[, 28:49]
#assign sample codes to row names
rownames(spe) <- all_data[, 1]


#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#add extra space to the right of the plot so the legend fits
par(mar=c(5, 4, 4, 15), xpd=TRUE)
#plot the NMDS
plot(example_NMDS, col = "white")

#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))


#set the colour for each treatment
#colors=c(rep("#44AA99",5),rep("#117733",5), rep("#88CCEE",5),rep("#332288",5), rep("#AA4499", 5), rep("#882255", 5))
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 

text(-0.7,1.4, paste("Stress = ", round(example_NMDS$stress, 3)))
#point shapes
pchs<- c(rep(16, 5), rep(17, 5), rep(18, 5), rep(15, 5), rep(3, 5), rep(4, 5)) 
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
   #plot the sample IDs on the NMDS, with the colour specific to the treatment
  #  orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
  #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], display = "sites", pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7, air = 0.01)
    
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }
#specify legend manually
#legend(1.7,0.9, legend = c("Grassland Bracken Present", "Grassland Bracken Asbent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"))

legend(1.7,0.9, legend=c("Grassland Bracken Present", "Grassland Bracken Asbent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), pch = c(16, 17, 18, 3, 4))


#save the file using Export -> Save As Image -> Width = 655, Height = 500 

#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- d[,(2:3)]
#run the permanova
veg_permanova <- adonis2(spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray")
veg_permanova
#veg_permanova indicates that Habitat and Vegetation have significant effects, with habitat explainng 26.7% of the variation and Vegetation explaining 36.0 %

#run an anosim - be sure to run with grouping = d$Habitat
m_com <- as.matrix(spe)
ano = anosim(m_com, grouping = d$Habitat, distance = "bray", permutations = 9999)
#check output of anosim
ano
