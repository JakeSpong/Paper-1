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
library(broom) #for curve fitting when analysing WEOC quality
library(scales) # displays 100,000 as 100,000, not as 1e5
library(DHARMa) #to help interpret GLMs

# renv, already installed in RStudio, used for version control 

#### load the master data file ----
all_data <- readr::read_csv(
  here::here("Data", "1) 2025-03-19_all-variables-measured.csv"), show_col_types = FALSE
) 
all_data <- as.data.frame(all_data)

#### Map sample coordinates ----

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggmap)

#to get google satellite layer, Go to the Google Cloud Console.  Create a new project (or use an existing one). Enable the Maps Static API and Geocoding API.  Get your API key.


# Register your Google API key
register_google(key = "AIzaSyA42ZmX_RGv9dTA7PFy4UR0YIjxW3x0rUE")

# Load the UK map from the rnaturalearth package
uk_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == "United Kingdom")

#extract he sample coordinates from the dataframe
sample_coordinates <- data.frame(
  Longitude = all_data$LongitudeE,
  Latitude = all_data$LatitudeN,
  Habitat = all_data$Habitat,
  Vegetation = all_data$Vegetation
)
# Calculate the bounding box surrounding the 30 sample coordinates
longitude_range_sample <- range(sample_coordinates$Longitude)
latitude_range_sample <- range(sample_coordinates$Latitude)

# Plot the UK map with only the bounding box
UK_Hawes <- ggplot(data = uk_map) +
  geom_sf(fill = "lightblue") +  # UK map with light blue fill
  geom_rect(
    aes(xmin = longitude_range_sample[1], xmax = longitude_range_sample[2], ymin = latitude_range_sample[1], ymax = latitude_range_sample[2]), 
    color = "red", fill = NA, size = 2) +  # Add bounding box around the points
  theme_minimal()
#show the map
#show(UK_Hawes)
#save the UK map with Haweswater highlighted
#ggsave("Figures/UK_Haweswater_boundingbox.svg", width = 8, height = 6)


# Get the Google Satellite map for the zoomed-in area based on the 30 sample coordinates
satellite_map <- get_map(location = c(lon = mean(sample_coordinates$Longitude), 
                                      lat = mean(sample_coordinates$Latitude)),
                         zoom = 15,  # Adjust zoom level for the desired zoom
                         maptype = "satellite",  # Use the "satellite" map type
                         source = "google")
#plot the map
sample_coords_map <- ggmap(satellite_map) +
  geom_point(data = sample_coordinates, aes(x = Longitude, y = Latitude, 
                                            color = Habitat, shape = Vegetation), size = 2) +
  scale_color_manual(values = c("Grassland" = "green", "Heathland" = "purple", "Woodland" = "brown")) +  # Color by habitat
  scale_shape_manual(values = c("Bracken Present" = 16, "Bracken Absent" = 17)) +  # Shape by vegetation
  theme_minimal() + 
  theme(legend.position = "right") 
#show the map
#show(sample_coords_map)
# Save the final map to a file
#ggsave("Figures/sample_locations.svg", plot = last_plot(), width = 7.5, height = 6, dpi = 300)


#combine the two maps into a single figure
maps_figure <- ggarrange(UK_Hawes, sample_coords_map, 
                     labels = c("A", "B"),
                     ncol = 2, nrow = 1,
                     #the width of each panel of the multifigure plot
                     widths = c(2,5))
#show the plot in the Plots window
#show(maps_figure)
#save the figure
ggsave("Figures/Maps_panel_figure.svg", plot = last_plot(), width = 7.5, height = 6, dpi = 300)




#### Alpha diversity of vegetation ----

#subset the vegetation species data
veg_spe <- all_data[, 28:49]
#compute the diversity indices
#species richness
all_data$`Vegetation Richness` <- apply(veg_spe[,]>0,1,sum)
#calculate diversity
all_data$`Vegetation Shannon` <- diversity(veg_spe[,], "shannon")
all_data$`Vegetation Simpson` <- diversity(veg_spe[,], "simpson")
all_data$`Vegetation Pielou's Evenness` <- all_data$`Vegetation Shannon` / log(all_data$`Vegetation Richness`)

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

#Veg richness

#run the GLM
glm <- glm(all_data$`Vegetation Richness` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`Vegetation Richness` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Vegetation Richness` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Veg evenness

#run the GLM
glm <- glm(all_data$`Vegetation Pielou's Evenness` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`Vegetation Pielou's Evenness` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Vegetation Pielou's Evenness` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#veg Shannon


#run the GLM
glm <- glm(all_data$`Vegetation Shannon` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`Vegetation Shannon` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Vegetation Shannon` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)




#veg Simpson


#run the GLM
glm <- glm(all_data$`Vegetation Simpson` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`Vegetation Simpson` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Vegetation Simpson` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)

#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#### Beta diversity of vegetation ----

#extract the veg abundances
veg_spe <- all_data[, 28:49]
#assign sample codes to row names
rownames(veg_spe) <- all_data[, 1]


#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(veg_spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#add extra space to the right of the plot so the legend fits
#par(mar=c(5, 4, 4, 15), xpd=TRUE)
#plot the NMDS
plot(example_NMDS, col = "white")


#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))
#set the colour for each treatment
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 
#point shapes
pchs<- c(rep(15, 5), rep(0, 5), rep(16, 5), rep(1, 5), rep(17, 5), rep(2, 5)) 
text(-1.8,1.3, paste("Stress = ", round(example_NMDS$stress, 3)))

for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
   #plot the sample IDs on the NMDS, with the colour specific to the treatment
   # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
  #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }
#specify legend manually for text sample IDs
#legend(1.7,0.9, legend = c("Grassland Bracken Present", "Grassland Bracken Asbent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"))

#save the file using Export -> Save As Image -> Width = 655, Height = 500 
#save the legend as a separate file and combine in Inkscape
plot(example_NMDS, col = "white")
#legend for point codes
legend(-1.7,0.9, legend=c("Grassland Bracken Present", "Grassland Bracken Absent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), pch = c(15, 0,16,1,17,2))


#save the file using Export -> Save As Image -> Width = 655, Height = 500 

#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- all_data[,(3:4)]
#run the permanova
veg_permanova <- adonis2(veg_spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray", by = "terms")
veg_permanova
#veg_permanova indicates that Habitat and Vegetation have significant effects, with habitat explainng 26.7% of the variation and Vegetation explaining 36.0 %

#run an anosim - when grouping by habitat
ano = anosim(as.matrix(veg_spe), grouping = all_data$Habitat, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)
#run an anosim - when grouping by vegetation
ano = anosim(as.matrix(veg_spe), grouping = all_data$Vegetation, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)



#### WEOC quality data processing ----
d <- readr::read_csv(
  here::here("Data", "2) 2023-12-13_WEOC-Quality-Microplate-Reader_Run-2-diluted-10-fold.csv")) 
#trim off file header information
d <- d[10:107,]
#convert to dataframe
d <- as.data.frame(d)
#name columns
names(d) <- d[2,]
colnames(d)[2] <- "Sample ID"
#remove rows of dataframe which contained column names
d <- d[3:98,]
#remove first column
d <- d[,-1]
#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
#remove final 6 rows - these wells were empty
d <- d[1:90,]
#remove second column
d <- d[, -2]
#add character to column names else we can't covnert from wide to long
# adding suffix to column names  
colnames(d) <- paste("new",sep="_",colnames(d)) 
#convert the numbers (which are datatype character for some reason) to integers
d[,2:582] <- sapply(d[,2:582],as.numeric)

#take average of the three technical replicates
# Calculate the averages and add a new row
d_averaged <- d %>%
  group_by(d[1]) %>%
  summarise(across(everything(), mean)) %>%
  mutate(row_type = "average")

#add wavelength (nm) as a column
d_long <- pivot_longer(d_averaged, cols = new_220:new_800, names_to = "Wavelength (nm)", values_to = "Absorbance")
#remove 2nd column
d_long <- d_long[,-2]
#remove the characters we added so we can plot the data
d_long$`Wavelength (nm)` <- gsub("new_","",d_long$`Wavelength (nm)`)
colnames(d_long)[1] <- "Sample ID"
#convert the wavelength numbers (which are datatype character for some reason) to integers
d_long[,2] <- sapply(d_long[,2],as.numeric)
#standardise the absorbance
#read in the concnetraiton data
d_conc <- readr::read_csv(
  here::here("Data", "3) Soil WEOC data means.csv"))
#add the DOM concentration column to our absorbance dataframe
d_stnd <- merge(d_long, d_conc[, c("Sample ID", "NPOC")], by = "Sample ID", all.x = TRUE)
#divide absorbance by concentration to standardize the measurements
d_stnd$`Standardized Absorbance (L mg-1 cm-1)` <- d_stnd$Absorbance/d_stnd$`NPOC`


#save our processed data file
write.csv(d_stnd, "Data//4) Standardised Processed WEOC Absorbance Data.csv", row.names =FALSE)


#### WEOC Wavelength of Interest boxplots ----
#read in the processed, STANDARDISED absorbance data
d <- readr::read_csv(
  here::here("Data", "4) Standardised Processed WEOC Absorbance Data.csv")
) 
# repeat at wavelengths of 250 (aromaticity, apparent molecular weight), 254 (aromaticity), 260 (hydrophobic C content), 265 (relative abundance of functional groups), 272 (aromaticity), 280 (hydrophobic C content, humification index, apparent molecular size), 285 (humification index), 300 (characterization of humic substances), 340 (colour), 350 (apparent molecular size), 365 (aromaticity, apparent molecular weight), 400 (humic substances characterization), 436 (quality indicator), 465 (relative abundance of functional groups)

#list of wavelengths of interest
wavelength_of_interest <- list(254, 250, 254, 260, 265, 272, 280, 285, 300, 340, 350, 365, 400, 436, 465)

#function to plot data for DOM, requires dataframe d (with columns SampleID, Wavelength, Absorbance) and wavelenth of interest
DOMboxplotter <- function(d, wavelength){
  #define the wavelength of interest
  wavelength_of_interest <- wavelength
  #filter data to extract absorbance at wavelength of interest
  abs <- d %>%
    filter(`Wavelength (nm)` == wavelength_of_interest) %>% #filter for specific wavelength
    select(`Sample ID`, `Standardized Absorbance (L mg-1 cm-1)`) #select the relevant columns
  
  #add in habitat and vegetation factors
  abs$Habitat <- c(rep("Grassland",10), rep("Heathland",10),rep("Woodland",10))
  abs$Vegetation <- c(rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5))
  
  abs <- as.data.frame(abs)
  
  #create string for y axis label
  yaxis_label <- paste("Absorbance at ", wavelength_of_interest, "nm")
  
    #plot absorbance
  abs_bxp <- ggboxplot(abs, x = "Habitat", y = "`Standardized Absorbance (L mg-1 cm-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
    labs(y = yaxis_label) + theme(
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
      legend.title = element_blank()
    ) 
  
  #create file name for our plot.  Use paste0 so there are no spaces between each item in the list
  filename <- paste0(Sys.Date(), "_absorbance_at_", wavelength_of_interest, ".svg")
  #save our plot.  As this is a function, we need specify the entire file path
  ggsave(path = "C:/Users/jakef/Documents/York/Paper-1/Figures", filename, abs_bxp)
  
  
  print(wavelength_of_interest)
  #statistical analysis
  #run the GLM
  glm <- glm(abs$`Standardized Absorbance (L mg-1 cm-1)` ~ abs$Habitat*abs$Vegetation)
  summary(glm)
  #out model seems to fit well
  simulationOutput <- simulateResiduals(fittedModel = glm)
  plot(simulationOutput)
  
  
  #two way anova
  anova <- aov(abs$`Standardized Absorbance (L mg-1 cm-1)` ~ abs$Habitat*abs$Vegetation)
  #check homogeneity of variance
  plot(anova, 1)
  #levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
  leveneTest(abs$`Standardized Absorbance (L mg-1 cm-1)` ~ abs$Habitat*abs$Vegetation)
  #check normality.  
  plot(anova, 2)
  #conduct shapiro-wilk test on ANOVA residules
  #extract the residuals
  aov_residuals <- residuals(object = anova)
  #run shapiro-wilk test.  if p > 0.05 the data is normal
  shapiro.test(x = aov_residuals)
  #view the ANOVA results
  summary(anova)
  #tukey's test to identify significant interactions
  tukey <- TukeyHSD(anova)
  print(tukey)
  #compact letter display
  cld <- multcompLetters4(anova, tukey)
  print(cld)
}
#plot the absorbance boxplot at the following given wavelengths
for (wavelength in wavelength_of_interest){
  DOMboxplotter(d, wavelength)
}




#### SUVA (254nm) statistical analysis ----
d <- readr::read_csv(
  here::here("Data", "4) Standardised Processed WEOC Absorbance Data.csv"))

#list of wavelengths of interest
wavelength_of_interest <- list(254)
abs <- d %>%
  filter(`Wavelength (nm)` == wavelength_of_interest) %>% #filter for specific wavelength
  select(`Sample ID`, `Standardized Absorbance (L mg-1 cm-1)`) #select the relevant columns

#add in habitat and vegetation factors
abs$Habitat <- c(rep("Grassland",10), rep("Heathland",10),rep("Woodland",10))
abs$Vegetation <- c(rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5))

abs <- as.data.frame(abs)



#### Panel figure of soil pH, moisture, total soil carbon, and WEOC concentration ----

#plot pH
pH_bxp <- ggboxplot(all_data, x = "Habitat", y = 'pH', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = "Soil pH") + theme(
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

#plot soil moisture
wc_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Water content (% of wet soil mass)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = "Water Content (%)") + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    #remove legend
    legend.position = "none",
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

#plot carbon
tsc_bxp <- ggboxplot(all_data, x = "Habitat", y = '`Total Carbon (g kg-1)`', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = expression("Total Soil Carbon (g kg"^-1*")")) + theme(
    
    #remove x axis label, tickes, labels
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

#plot WEOC
WEOC_bxp <- ggboxplot(all_data, x = "Habitat", y = "`DOC Concentration (mg C g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = expression("WEOC Concentration (mg C g"^-1*")")) + theme(
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

#save 770 wide, 640 high
panel_bxp <- ggarrange(pH_bxp, wc_bxp, tsc_bxp, WEOC_bxp, 
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 2,
                       common.legend = TRUE, legend="top")

show(panel_bxp)
#save the image
ggsave(path = "Figures", paste0(Sys.Date(), "_4-panel-pH-moisture-TC-WEOC.svg"), panel_bxp, width = 8, height = 7, dpi = 300)

#soil pH

#run the GLM
glm <- glm(all_data$pH ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$pH ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$pH ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#soil water content

#run the GLM
glm <- glm(all_data$`Water content (% of wet soil mass)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`Water content (% of wet soil mass)` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Water content (% of wet soil mass)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#soil total carbon

#run the GLM
glm <- glm(all_data$`Total Carbon (g kg-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`Total Carbon (g kg-1)` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Total Carbon (g kg-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#soil WEOC

#run the GLM
glm <- glm(all_data$`DOC Concentration (mg C g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(all_data$`DOC Concentration (mg C g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`DOC Concentration (mg C g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#view the ANOVA results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#### Panel figure of DOM alpha parameter, total soil nitrogen, WEN, C:N ratio ----

#WEOC quality analysis - alpha parameter calculation and plotting 
#quality needs to be analysed; easiest thing to do is: check data for nice curve, fit 2 component exponential decay curve through data.  Get slope from this (spectral slope), and the alpha paramter from this slope
#read in the processed, STANDARDISED absorbance data
d <- readr::read_csv(
  here::here("Data", "4) Standardised Processed WEOC Absorbance Data.csv")) 
#eliminate all absorbance values above 600 nm
d <- d %>% filter(`Wavelength (nm)` <= 600)

# Fit a two-component exponential decay curve through the data for all our curves
fitted <- d %>%
  nest(data = -`Sample ID`) %>%
  mutate(
    fit = map(data, ~nls(`Standardized Absorbance (L mg-1 cm-1)` ~ SSasymp(`Wavelength (nm)`, yf, y0, log_alpha), data = .)),
    tidied = map(fit, tidy),
    augmented = map(fit, augment),
  )
# Produce a table of fit parameters: y0, yf, alpha
table <- fitted %>% 
  unnest(tidied) %>% 
  select(`Sample ID`, term, estimate) %>% 
  spread(term, estimate) %>% 
  mutate(alpha = exp(log_alpha))
#display table of fit parameters
table

#plot each absorbance curve along with the line of best fit
augmented <- fitted %>% 
  unnest(augmented)
qplot(`Wavelength (nm)`, `Standardized Absorbance (L mg-1 cm-1)`, data = augmented, geom = 'point', colour = `Sample ID`) +
  geom_line(aes(y=.fitted))

#now we have extracted the paramters of our lines of best fit, add the descriptors
table$Habitat <- c(rep("Grassland",10), rep("Heathland",10),rep("Woodland",10))
table$Vegetation <- c(rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5))

#boxplot the alpha, which describes the curve ie how quicky we go from low wavelength (high mass C compounds) to high wavelength (low mass C compounds)
alpha_bxp <- ggboxplot(table, x = "Habitat", y = 'alpha', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = "DOM fitted curve alpha parameter") + theme(
    #remove x axis label, ticks, labels
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
    legend.title = element_blank()
  ) 
show(alpha_bxp)
#save our plot
#ggsave(path = "Figures", paste0(Sys.Date(), '_standardised_DOM-curve-alpha-paramter_black-green.svg'), alpha_bxp)




#run the GLM
glm <- glm(table$alpha ~ table$Habitat * table$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA
anova <- aov(table$alpha ~ table$Habitat * table$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(table$alpha ~ table$Habitat * table$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)



#plot nitrogen
tsn_bxp <- ggboxplot(all_data, x = "Habitat", y = '`Total Nitrogen (g kg-1)`', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(x = "Habitat x Vegetation",
       y = expression("Total Soil Nitrogen (g kg"^-1*")")) + theme(
         #remove x axis label, tickes, labels
         #remove x axis label, ticks, labels
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
         legend.title = element_blank()
       ) 

show(tsn_bxp)
#save our plot
#ggsave(path = "Figures", paste0(Sys.Date(), '_total-soil-nitrogen_black-green.svg'), tsn_bxp)

#run the GLM. We sqrt transform the data in order to pass Shapiro Wilk.  Though the results of the ANOVA
#are unchanged
glm <- glm(sqrt(all_data$`Total Nitrogen (g kg-1)`) ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  We sqrt transform the data in order to pass Shapiro Wilk.  Though the results of the ANOVA
#are unchanged
anova <- aov(sqrt(all_data$`Total Nitrogen (g kg-1)`) ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(sqrt(all_data$`Total Nitrogen (g kg-1)`) ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# Water extractable nitrogen boxplots 

#graphing boxplots. Units are ppm aka mg per kg
tnb_bxp <- ggboxplot(all_data, x = "Habitat", y = "`TNb Concentration (mg N g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = expression("TNb Concentration (mg N g"^-1*")")) + theme(
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
    axis.line = element_line(colour = "black", linewidth = 0.75),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
    legend.title = element_blank()
  ) 

show(tnb_bxp)
#save our plot
#ggsave(path = "Figures", paste0(Sys.Date(), "_tnb_black-green.svg"), tnb_bxp)



#run the GLM
glm <- glm(all_data$`TNb Concentration (mg N g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  
anova <- aov(all_data$`TNb Concentration (mg N g-1)` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`TNb Concentration (mg N g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# C:N analysis
#plot carbon:nitrogen ratio
cnr_bxp <- ggboxplot(all_data, x = "Habitat", y = 'CN ratio', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(x = "Habitat x Vegetation",
       y = "C:N ratio") + theme(
         #remove x axis label, tickes, labels
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
         legend.title = element_blank()
       ) 
show(cnr_bxp)
#save our plot
#ggsave(path = "Figures", paste0(Sys.Date(), '_total-soil-carbon-nitrogen-ratio_black-green.svg'), cnr_bxp)

#run the GLM
glm <- glm(all_data$`CN ratio` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  We sqrt transform the data in order to pass Shapiro Wilk.  Though the results of the ANOVA
#are unchanged
anova <- aov(all_data$`CN ratio` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`CN ratio` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)




#save 770 wide, 640 high
panel_bxp <- ggarrange(alpha_bxp, tsn_bxp, tnb_bxp, cnr_bxp, 
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 2,
                       common.legend = TRUE, legend="top")

show(panel_bxp)
#save the image
ggsave(path = "Figures", paste0(Sys.Date(), "_4-panel-alpha-TN-WEN-CNratio.svg"), panel_bxp, width = 8, height = 7, dpi = 300)





#### standard deviations of WEN ----
# Group by habitat and vegetation, then calculate standard deviation of nitrogen
df_summary <- all_data %>%
  group_by(Habitat, Vegetation) %>%
  summarise(sd_nitrogen = sd(`TNb Concentration (mg N g-1)`, na.rm = TRUE))

# View result
print(df_summary)

#### Soil cation boxplots ----

#graphing boxplots. Units are ppm aka mg per kg
Al_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Al (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("Al (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Zn_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Zn (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("Zn (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Mn_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Mn (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("Mn (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Mg_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Mg (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("Mg (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Ca_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Ca (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("Ca (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 


#graphing boxplots. Units are ppm aka mg per kg
Na_bxp <- ggboxplot(all_data, x = "Habitat", y = "`Na (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("Na (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 



#graphing boxplots. Units are ppm aka mg per kg
K_bxp <- ggboxplot(all_data, x = "Habitat", y = "`K (mg g-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y =  expression("K (mg g"^-1*")"))+ theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 

panel_bxp <- ggarrange(Na_bxp, K_bxp, Ca_bxp, Mg_bxp, Zn_bxp, Mn_bxp, Al_bxp,
                       labels = c("A", "B", "C", "D", "E", "F", "G"),
                       ncol = 2, nrow = 4,
                       common.legend = TRUE, legend="top")

show(panel_bxp)
ggsave(path = "Figures", paste0(Sys.Date(), "_all_ions_bxps.svg"), panel_bxp, width = 7, height = 7, dpi= 300)


#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_Al_bxp.svg"), Al_bxp)
#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_Zn_bxp.svg"), Zn_bxp)
#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_Mn_bxp.svg"), Mn_bxp)
#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_Mg_bxp.svg"), Mg_bxp)
#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_Ca_bxp.svg"), Ca_bxp)
#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_Na_bxp.svg"), Na_bxp)
#save our plot.  As this is a function, we need specify the entire file path
#ggsave(path = "C:/Users/jakef/Documents/York/Project 1 Analysis/project-1/figures", paste0(Sys.Date(), "_K_bxp.svg"), K_bxp)

#Na

#run the GLM
glm <- glm(all_data$`Na (mg g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  
anova <- aov(all_data$`Na (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Na (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#K

#run the GLM
glm <- glm(all_data$`K (mg g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  
anova <- aov(all_data$`K (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`K (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Ca

#run the GLM
glm <- glm(sqrt(all_data$`Ca (mg g-1)`) ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  We sqrt transform the data in order to pass Shapiro Wilk.  
anova <- aov(all_data$`Ca (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Ca (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Mg

#run the GLM
glm <- glm(all_data$`Mg (mg g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  
anova <- aov(all_data$`Mg (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Mg (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Zn

#run the GLM
glm <- glm(all_data$`Zn (mg g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA. 
anova <- aov(all_data$`Zn (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Zn (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Mn

#run the GLM
glm <- glm(all_data$`Mn (mg g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  
anova <- aov(all_data$`Mn (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Mn (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Al

#run the GLM
glm <- glm(all_data$`Al (mg g-1)` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)


#two-way ANOVA.  We sqrt transform the data in order to pass Levene .  If we don't, we get p = 0.05 for habitat, so not significant anyway
anova <- aov(all_data$`Al (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Al (mg g-1)` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)
#look at the results
summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)



#### pH vs Al linear modelling ----

#plot the relationship
plot <- ggplot(all_data, aes(x = pH, y = `Al (mg g-1)`)) +
  geom_point(color = "black", size = 2) + # Add points
  geom_smooth(method = "lm", color = "red", se = TRUE) + # Add regression line
  # Line of best fit
  annotate("text", x = 5.3, y = max(all_data$`Al (mg g-1)`), 
           label = paste("y =", 
                         round(coef(model)[2], 2), "x +", 
                         round(coef(model)[1], 2), 
                         "\nR² =", 
                         round(summary(model)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(
    x = "pH",
    y = expression("[Al] (mg g"^-1*")")
  ) +  theme(
    
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
  ) 

show(plot)

ggsave(path = "Figures", paste0(Sys.Date(), "_Al-vs-pH.svg"), plot)

#linear model of the relationship
model <- lm(all_data$`Al (mg g-1)` ~ all_data$pH)
summary(model)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)


#### Soil cation summary statistics ----


## Calculate summary statistics for each ion
summary_stats_ions <- all_data %>%
  gather(key = "`(mg g-1)`", value = "Replicate_Measurement", ends_with("`(mg g -1)`")) %>%
  group_by(Habitat, Vegetation, `(mg g-1)`) %>%
  summarise(
    Mean = mean(Replicate_Measurement),
    SD = sd(Replicate_Measurement),
  )

#export as svg
summary_stats_ions %>%
  readr::write_csv(file = fs::path("Data", "5) Soil cation summary statistics.csv")) #NEEDS FIXING #NEEDS FIXING

#ALL OF THE BELOW HAVE YET TO HAVE NEW STATS DONE ON THEM!


#### Microinvertebrate abundances to a standard depth ----

#extract the morphospecies count data
morpho_spe <- all_data[, 51:436]
#calculate the abundance of morphotypes in each sample
all_data$`Morphotype abundance` <- rowSums(morpho_spe)
#print the total number of morphotypes
print(sum(all_data$`Morphotype abundance`))
#volume of the soil corer
volume <- 0.1 * (pi*(0.025^2))
all_data$CoreVolume <- volume
#number of mesofauna per m3, to a depth of 10cm
all_data$`Individuals per m2 to 10 cm depth` <- ((all_data$`Morphotype abundance`/all_data$CoreVolume)*0.1)/1000


#number of mite morphospecies
all_data <- all_data %>%
  mutate(Mites = rowSums(select(morpho_spe, contains("Mite") | contains("Mstg") | contains("Noth") | contains("Astig") | contains("Box") | contains("Turt")| contains("Tank")| contains("Pter")| contains("Grp")) > 0))
#number of collembola morphospecies
all_data <- all_data %>%
  mutate(Collembola = rowSums(select(morpho_spe, contains("Sprg") | contains("Neel") | contains("Smin") ) > 0))
#other morphospecies
all_data <- all_data %>%
  mutate(`Other Morphospecies` = rowSums(select(morpho_spe, contains("Hyme") | contains("Larv") | contains("Cent") | contains("Mill") | contains("Ant") | contains("Aph")| contains("Wdl")| contains("Ench")| contains("Fly")| contains("Spid")| contains("Tick")| contains("Bug")| contains("Pseu")) > 0))

##analyse mites, springtails and other in each treatment

# Create a list of group codes (or patterns to match)
codes <- c("GBP", "GBA", "HBP", "HBA", "WBP", "WBA")
# Initialize an empty list to store results
group_summary <- list()
# Loop through each code and calculate mean and SD for rows where 'ID' contains the code
for(code in codes) {
  # Filter rows where ID contains the pattern (code)
  group_data <- all_data %>% filter(grepl(code, `Sample ID`))
  
  # Calculate mean and standard deviation for the target_column
  summary_stats <- group_data %>%
    summarise(  #change first argument to column of choice
      mean_value = mean(`Mites`, na.rm = TRUE),
      sd_value = sd(`Mites`, na.rm = TRUE)
    )
  
  # Store the result in the list
  group_summary[[code]] <- summary_stats
}
# Combine the results into a single data frame
df_summary <- bind_rows(group_summary, .id = "group_code")
# View the result
print(df_summary)

#get percentages of mites/springtails of all morphospecies
all_data$`Percentage Mites` = (all_data$Mites/(all_data$Mites + all_data$Collembola +all_data$`Other Morphospecies`))*100
all_data$`Percentage Collembola` = (all_data$Collembola/(all_data$Mites + all_data$Collembola +all_data$`Other Morphospecies`))*100
all_data$`Percentage Other Morphospecies` = (all_data$`Other Morphospecies`/(all_data$Mites + all_data$Collembola +all_data$`Other Morphospecies`))*100

#we are printing multiple lines so use cat()
cat("Min mite percentage: ", min(all_data$`Percentage Mites`), "Max mite percentage: ", max(all_data$`Percentage Mites`),"Min collembola percentage: ", min(all_data$`Percentage Collembola`), "Max collembola percentage: ", max(all_data$`Percentage Collembola`), "Min other morphospecies percentage: ", min(all_data$`Percentage Other Morphospecies`), "Max other morphospecies: ", max(all_data$`Percentage Other Morphospecies`))

#plot abundances
figure <- ggboxplot(all_data, x = "Habitat", y = 'Individuals per m2 to 10 cm depth', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = expression("1000 individuals per m"^2*" soil (to 10 cm depth)")) + theme( #remove x axis label
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
    legend.title = element_blank()
  ) + scale_y_continuous(labels = label_comma())

#display our plot
figure
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), '_morphotype-abundances_black-green.svg'), figure, width = 5, height = 4, dpi = 300)


#run the GLM
glm <- glm(all_data$`Individuals per m2 to 10 cm depth` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

#two-way ANOVA.  
anova <- aov(all_data$`Individuals per m2 to 10 cm depth` ~ as.factor(all_data$Habitat) * as.factor(all_data$Vegetation))
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Individuals per m2 to 10 cm depth` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#figure +
 # annotate("text", x = 1.5, y = 300, label = "a   b    ab", color = "black") 
#save our plot
#ggsave(path = "Figures", paste0(Sys.Date(), '_cld-letters-for-inkscape.pdf'), plot = last_plot(), width = 5, height = 4, dpi = 300, device = "pdf")

#### Alpha diversity (Richess, evenness, Shannon and Simpson Diversity) of mesofauna groups ----

#extract the morphospecies count data
morpho_spe <- all_data[, 51:436]
#species richness
all_data$`Morphotype Richness` <- apply(morpho_spe[,]>0,1,sum)
#calculate diversity
all_data$`Morphotype Shannon` <- diversity(morpho_spe[,], "shannon")
all_data$`Morphotype Simpson` <- diversity(morpho_spe[,], "simpson")
all_data$`Morphotype Pielou's Evenness` <- all_data$`Morphotype Shannon` / log(all_data$`Morphotype Richness`)

#graphing abundances
evenness_bxp <- ggboxplot(all_data, x = "Habitat", y = "Morphotype Pielou's Evenness", color = "Vegetation", ylab = "Pielou's Evenness", palette = c("black", "limegreen"), lwd = 0.75) + theme(
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
  legend.title = element_blank()
) 
#boxplot of grassland vs heathland, species richness
richness_bxp <- ggboxplot(all_data, x = "Habitat", y = "Morphotype Richness", color = "Vegetation", ylab = "Morphotype \n Richness", palette = c("black", "limegreen"), lwd = 0.75) + theme(
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
  legend.title = element_blank()
) 
#graphing boxplots with bracken split from nonbracken
shannon_bxp <-ggboxplot(all_data, x = "Habitat", y = "Morphotype Shannon", color = "Vegetation", ylab = "Shannon Diversity", palette = c("black", "limegreen"), lwd = 0.75) + theme(
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
  legend.title = element_blank()
) 
#graphing boxplots with bracken split from nonbracken
simpson_bxp <- ggboxplot(all_data, x = "Habitat", y = "Morphotype Simpson", color = "Vegetation", ylab = "Simpson Diversity", palette = c("black", "limegreen"), lwd = 0.75) + theme(
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
  legend.title = element_blank()
) 

all_bxp <- ggarrange(richness_bxp, evenness_bxp, shannon_bxp, simpson_bxp, 
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, legend="top")
show(all_bxp)
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), '_alpha-diversity_black-green.svg'), width = 7, height = 5, all_bxp)

#Morphotype richness

#run the GLM
glm <- glm(all_data$`Morphotype Richness` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

#two-way ANOVA.  
anova <- aov(all_data$`Morphotype Richness` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Morphotype Richness` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#Morphotype evenness

#run the GLM
glm <- glm(all_data$`Morphotype Pielou's Evenness` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

#two-way ANOVA.  
anova <- aov(all_data$`Morphotype Pielou's Evenness` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Morphotype Pielou's Evenness` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#Morphotype Shannon

#run the GLM
glm <- glm(all_data$`Morphotype Shannon` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

#two-way ANOVA.  
anova <- aov(all_data$`Morphotype Shannon` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Morphotype Shannon` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

#Morphotype Simpson

#run the GLM
glm <- glm(all_data$`Morphotype Simpson` ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

#two-way ANOVA.  
anova <- aov(all_data$`Morphotype Simpson` ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(all_data$`Morphotype Simpson` ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)




#### NMDS analysis of week 1 only mesofauna abundance data ----

# Load mesofauna abundance data

#file containg all morphotyes, not just springtails and mites
filename <- "6) Week-1-morphotypes.csv"
#csv containing only mites and springtails
#filename <- "mites-springtails-only week-1-mesofauna-corrected.csv"
d <- readr::read_csv(here::here("Data", filename), show_col_types = FALSE) 
#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#remove all empty rows
d <- d[1:30,]
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0
#replace row index with sample names
rownames(d) <- d[,1]
#just the morphospecies counts
spe <- d[,-(1:4)]
spe <- as.matrix(spe)
#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#plot the NMDS
plot(example_NMDS, col = "white")
#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))
#set the colour for each treatment
#colors =c(rep("#44AA99",5),rep("#117733",5), rep("#88CCEE",5),rep("#332288",5), rep("#AA4499", 5), rep("#882255", 5)) 
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 
#shapes for point codes
pchs<- c(rep(15, 5), rep(0, 5), rep(16, 5), rep(1, 5), rep(17, 5), rep(2, 5))
#display the stress
text(-0.6,2.2, paste("Stress = ", round(example_NMDS$stress, 3)))
#visualise the points and ellipses
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment 
    
    #ellipses are of the kind se = standard error (or sd = standard deviation) at a 95% confidence
    ordiellipse(example_NMDS$point[grep(i,treat),], kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }


#save the file using Export -> Save As Image -> Width = 655, Height = 500 

#save the legend as its own plot, to be modified in Inkscape later
#plot the NMDS
plot(example_NMDS, col = "white")
#legend for text codes
#legend(-3.5,-1, legend = c("Grassland Bracken", "Grassland Non-bracken", "Heathland Bracken", "Heathland Non-bracken", "Woodland Bracken", "Woodland Non-Bracken"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), bty = "n")
#legend for point codes
legend(-1.7,0.9, legend=c("Grassland Bracken Present", "Grassland Bracken Absent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), pch = c(15, 0,16,1,17,2))

#ggsave(path = "Figures", paste0(Sys.Date(), '_morphospecies-nmds.svg'), last_plot(), width = 8, height = 6, dpi = 300)

# do PERMANOVA analysis


#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- all_data[,(3:4)]
#run the permanova
morph_permanova <- adonis2(spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray", by = "terms")
morph_permanova


#run an ANOSIM. The ANOSIM test is similar to an ANOVA hypothesis test, but it uses a dissimilarity matrix as input instead of raw data. It is also non-parametric, meaning it doesn’t assume much about your data (like normal distribution etc), so it’s a good bet for often-skewed microbial abundance data. As a non-parametric test, ANOSIM uses ranked dissimilarities instead of actual distances, and in this way it’s a very nice complement to an NMDS plot. The main point of the ANOSIM test is to determine if the differences between two or more groups are significant.
#run an anosim - when grouping by habitat
ano = anosim(as.matrix(spe), grouping = all_data$Habitat, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)
#run an anosim - when grouping by vegetation
ano = anosim(as.matrix(spe), grouping = all_data$Vegetation, permutations = 9999, distance = "bray")
# When interpreting these results you want to look at the ANOSIM statistic R and the Significance values. A Significance value less than 0.05 is generally considered to be statistically significant, and means the null hypothesis can be rejected. “The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups” (GUSTAME). In other words, the higher the R value, the more dissimilar your groups are in terms of microbial community composition.
ano
plot(ano)

#### NMDS analysis of week 1 + 2 mesofauna abundance data ----

#only the mite and springtail morphospecies
#filename <- "mites-springtails-only_week-1+2-mesofauna-corrected.csv"
#all morphospecies
spe <- all_data[,51:436]
#replace row index with sample names
rownames(spe) <- all_data[,1]
spe <- as.matrix(spe)

#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#plot the NMDS
plot(example_NMDS, col = "white")


#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))
#set the colour for each treatment
#colors =c(rep("#44AA99",5),rep("#117733",5), rep("#88CCEE",5),rep("#332288",5), rep("#AA4499", 5), rep("#882255", 5)) 
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 
#shapes for point codes
pchs<- c(rep(15, 5), rep(0, 5), rep(16, 5), rep(1, 5), rep(17, 5), rep(2, 5))
#display the stress
text(-0.8,1.4, paste("Stress = ", round(example_NMDS$stress, 3)))
#visualise the points and ellipses
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }



#overlay mesofauna morphotypes (instrinsic variables)
meso.intfit <- envfit(example_NMDS, spe, permutations = 999)
dev.new()
ordiplot(example_NMDS, type = "n", main = "intrinsic variables")
#orditorp(example_NMDS, display = "sites", labels = F, pch = c(16, 8, 17, 18) [as.numeric(env$`CN ratio`)], col = c("green", "blue", "orange", "black") [as.numeric(env$`CN ratio`)], cex = 1)
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }

plot(meso.intfit, col = "black", cex = 0.7)




#overlay environmental variables
env <- all_data[, 12:20]
rownames(env) <- all_data[, 1]
meso.envfit <- envfit(example_NMDS, env, permutations = 999)

ordiplot(example_NMDS, type = "n")
#orditorp(example_NMDS, display = "sites", labels = F, pch = c(16, 8, 17, 18) [as.numeric(env$`CN ratio`)], col = c("green", "blue", "orange", "black") [as.numeric(env$`CN ratio`)], cex = 1)
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }

plot(meso.envfit, col = "black", cex = 0.7)


#save the file using Export -> Save As Image -> Width = 655, Height = 500 

# do PERMANOVA analysis
#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- all_data[,(3:4)]
#run the permanova
morph_permanova <- adonis2(spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray", by = "terms")
morph_permanova


#run an ANOSIM. The ANOSIM test is similar to an ANOVA hypothesis test, but it uses a dissimilarity matrix as input instead of raw data. It is also non-parametric, meaning it doesn’t assume much about your data (like normal distribution etc), so it’s a good bet for often-skewed microbial abundance data. As a non-parametric test, ANOSIM uses ranked dissimilarities instead of actual distances, and in this way it’s a very nice complement to an NMDS plot. The main point of the ANOSIM test is to determine if the differences between two or more groups are significant.
#run an anosim - when grouping by habitat
ano = anosim(as.matrix(spe), grouping = all_data$Habitat, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)
#run an anosim - when grouping by vegetation
ano = anosim(as.matrix(spe), grouping = all_data$Vegetation, permutations = 9999, distance = "bray")
# When interpreting these results you want to look at the ANOSIM statistic R and the Significance values. A Significance value less than 0.05 is generally considered to be statistically significant, and means the null hypothesis can be rejected. “The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups” (GUSTAME). In other words, the higher the R value, the more dissimilar your groups are in terms of microbial community composition.
ano
plot(ano)



#### Morphotype nestedness and turnover ----

#just the morphospecies counts
species_matrix <- all_data[,(51:436)]
rownames(species_matrix) <-all_data[,1]
species_matrix <- as.matrix(species_matrix)
# Convert the matrix to a presence-absence (binary) matrix if not already
species_binary <- ifelse(species_matrix > 0, 1, 0)

#run nestedness/turnover analysis on each of the treatments
#grassland bracken present
gb <- species_binary[(1:5),]
#grassland bracken absent
gn <- species_binary[(6:10),]
#heathland bracken present
hb <- species_binary[(11:15),]
#heathland bracken absent
hn <- species_binary[(16:20),]
#woodland bracken present
wb <- species_binary[(21:25),]
#woodland bracken absent
wn <- species_binary[(26:30),]

#get betapart object
gb.core <- betapart.core(gb)
gn.core <- betapart.core(gn)
hb.core <- betapart.core(hb)
hn.core <- betapart.core(hn)
wb.core <- betapart.core(wb)
wn.core <- betapart.core(wn)

# multiple site measures
gb.multi <- beta.multi(gb.core, index.family = "jaccard")
gn.multi <- beta.multi(gn.core, index.family = "jaccard")
hb.multi <- beta.multi(hb.core, index.family = "jaccard")
hn.multi <- beta.multi(hn.core, index.family = "jaccard")
wb.multi <- beta.multi(wb.core, index.family = "jaccard")
wn.multi <- beta.multi(wn.core, index.family = "jaccard")

#compile the measures into a table. JTU = value of the turnover component, measured as turnover fraction of Jaccard dissimilarity.  JNE = value of the nestedness component, measured as nestedness-resultant faction of Jaccard dissimilarity.  .JAC = value of the overall beta diversity, measured as the Jaccard dissimilarity

#rownames
habitats <- c("Grassland", "Grassland", "Heathland", "Heathland", "Woodland", "Woodland")
bracken <- c("Present", "Absent","Present", "Absent","Present", "Absent")
#put all the dissimilarity measures into a table
# Combine lists into a data frame, each list becomes a row
beta_table <- as.data.frame(rbind(gb.multi, gn.multi, hb.multi, hn.multi, wb.multi, wn.multi))
#ensure all data is numeric
beta_table <- as.data.frame(lapply(beta_table, function(x) as.numeric(as.character(x))))
# Combine with new columns on the left
beta_table <- cbind(Habitat = habitats, Bracken = bracken, beta_table)
beta_table <- as.data.frame(beta_table)



# # sampling across a given number of sites, with "samples" being the number of random samples used to calculate the distribution of dissimilarity measures
# gb.samp <- beta.sample(gb.core, index.family = "jaccard",
#                        sites=5, samples=100)
# gn.samp <- beta.sample(gn.core, index.family = "jaccard",
#                        sites=5, samples=100)
# hb.samp <- beta.sample(hb.core, index.family = "jaccard",
#                        sites=5, samples=100)
# hn.samp <- beta.sample(hn.core, index.family = "jaccard",
#                        sites=5, samples=100)
# wb.samp <- beta.sample(wb.core, index.family = "jaccard",
#                        sites=5, samples=100)
# wn.samp <- beta.sample(wn.core, index.family = "jaccard",
#                        sites=5, samples=100)
# 
# # plotting the distributions of components
# dist.gb <- gb.samp$sampled.values
# dist.gn <- gn.samp$sampled.values
# dist.hb <- hb.samp$sampled.values
# dist.hn <- hn.samp$sampled.values
# dist.wb <- wb.samp$sampled.values
# dist.wn <- wn.samp$sampled.values
# 
# #Multi-site dissimilarities
# #plots total dissimilarity for each treatment
# plot(density(dist.gb$beta.JAC), col = "#999999", xlim = c(0, 2), ylim = c(0, 2), xlab=expression(beta[JAC]), main='', lwd=3)
# lines(density(dist.gn$beta.JAC), col = "#E69F00", lwd=3) 
# lines(density(dist.hb$beta.JAC), col = "#56B4E9", lwd=3) 
# lines(density(dist.hn$beta.JAC), col = "#009E73", lwd=3) 
# lines(density(dist.wb$beta.JAC), col = "#CC79A7", lwd=3) 
# lines(density(dist.wn$beta.JAC), col = "#0072B2", lwd=3) 
# #legend(1.5, 2.0, legend = c("Grassland Bracken", "Grassland Non-bracken", "Heathland Bracken", "Heathland Non-bracken", "Woodland Bracken", "Woodland Non-Bracken"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7","#0072B2"), bty = "n")
# 
# #plot nestedness separately due to large difference in values compared to turnover and total dissimilarity
# plot(density(dist.gb$beta.JNE), col = "#999999", xlim = c(0, 0.1), ylim = c(0, 100), xlab=expression(beta[JNE]), main='', lty = 3, lwd=3)
# lines(density(dist.gn$beta.JNE), col = '#E69F00',lty=3, lwd=2) 
# lines(density(dist.hb$beta.JNE), col = '#56B4E9',lty=3, lwd=2) 
# lines(density(dist.hn$beta.JNE), col = '#009E73',lty=3, lwd=2)
# lines(density(dist.wb$beta.JNE), col = '#CC79A7',lty=3, lwd=2) 
# lines(density(dist.wn$beta.JNE), col = '#0072B2',lty=3, lwd=2) 
# 
# plot(density(dist.gb$beta.JTU), col = "#999999", xlim = c(0, 2), ylim = c(0, 2), xlab=expression(beta[JTU]), main='', lty = 2,lwd=3)
# #plots turnover-resultant dissimilarity for grassland bracken sites
# lines(density(dist.gn$beta.JTU), col = '#E69F00',lty=2, lwd=2)
# lines(density(dist.hb$beta.JTU), col = '#56B4E9',lty=2, lwd=2)
# lines(density(dist.hn$beta.JTU), col = '#009E73',lty=2, lwd=2)
# lines(density(dist.wb$beta.JTU), col = '#CC79A7',lty=2, lwd=2)
# lines(density(dist.wn$beta.JTU), col = '#0072B2',lty=2, lwd=2)
# 


# 
# #pairwise
# pair.gb <- beta.pair(gb, index.family = "jaccard")
# # plotting clusters
# dist.s <- gb.samp$sampled.values
# plot(hclust(pair.gb$beta.jtu, method="average"), hang=-1, main='', sub='', xlab='')
# title(xlab=expression(beta[jtu]), line=0.3)
# plot(hclust(pair.gb$beta.jne, method="average"), hang=-1, main='', sub='',, xlab='')
# title(xlab=expression(beta[jne]), line=0.3)
# 
# 
# #run dissimilarity on bracken sites vs non bracken sites
# 
# bracken <- species_binary[c(1,2,3,4,5,11,12,13,14,15, 21,22,23,24,25), ]
# nonbracken <- species_binary[c(6,7,8,9,10, 16, 17, 18, 19, 20, 26,27,28,29,30),]
# pair.bracken<- beta.pair(bracken, index.family = "jaccard")
# pair.nonbracken<- beta.pair(nonbracken, index.family = "jaccard")
# # plotting clusters
# #bracken
# plot(hclust(pair.bracken$beta.jtu, method="average"), hang=-1, main='', sub='', xlab='')
# title(xlab=expression(beta[jtu - bracken]), line=0.3)
# plot(hclust(pair.bracken$beta.jne, method="average"), hang=-1, main='', sub='',, xlab='')
# title(xlab=expression(beta[jne - bracken]), line=0.3)
# #non-bracken
# plot(hclust(pair.nonbracken$beta.jtu, method="average"), hang=-1, main='', sub='', xlab='')
# title(xlab=expression(beta[jtu] - nonbracken), line=0.3)
# plot(hclust(pair.nonbracken$beta.jne, method="average"), hang=-1, main='', sub='',, xlab='')
# title(xlab=expression(beta[jne - nonbracken]), line=0.3)
# 
# #the correlation between the two clustered datasets (bracken and non-bracken)
# #cor(cophenetic(as.dendrogram(hclust(pair.nonbracken$beta.jtu, method="average"))), cophenetic(as.dendrogram(hclust(pair.bracken$beta.jne, method="average"))))
# #get betapart object
# species_binary.core <- betapart.core(species_binary)
# 
# ## USING SORENSEN DISSIMILARITY
# 
# # multiple site measures
# species_binary.multi <- beta.multi(species_binary.core, index.family = "sorensen")
# # sampling across a given number of sites, with "samples" being the number of random samples used to calculate the distribution of dissimilarity measures
# species_binary.samp <- beta.sample(species_binary.core, index.family = "sorensen",
#                                    sites=30, samples=100)
# # plotting the distributions of components
# dist.s <- species_binary.samp$sampled.values
# plot(density(dist.s$beta.SOR), xlim = c(0, 3), ylim = c(0, 70), xlab='Beta
# diversity (Sorensen)', main='', lwd=3)
# lines(density(dist.s$beta.SNE), lty=1, lwd=2) 
# lines(density(dist.s$beta.SIM), lty=2, lwd=2)
# 
# # pairwise  using sorensen dissimilarity
# pair.s <- beta.pair(species_binary, index.family = "sorensen")
# # plotting clusters
# dist.s <- species_binary.samp$sampled.values
# plot(hclust(pair.s$beta.sim, method="average"), hang=-1, main='', sub='', xlab='')
# title(xlab=expression(beta[sim]), line=0.3)
# plot(hclust(pair.s$beta.sne, method="average"), hang=-1, main='', sub='',, xlab='')
# title(xlab=expression(beta[sne]), line=0.3)
# 
# 
# 
# 
# ## USING JACCARD DISSIMILARITY
# 
# # multiple site measures
# species_binary.multi <- beta.multi(species_binary.core, index.family = "jaccard")
# # sampling across a given number of sites, with "samples" being the number of random samples used to calculate the distribution of dissimilarity measures
# species_binary.samp <- beta.sample(species_binary.core, index.family = "jaccard",
#                                    sites=30, samples=100)
# 
# # plotting the distributions of components
# dist.s <- species_binary.samp$sampled.values
# plot(density(dist.s$beta.JAC), xlim = c(0, 3), ylim = c(0, 150), xlab='Beta
# diversity (Jaccard)', main='', lwd=3)
# lines(density(dist.s$beta.JNE), lty=1, lwd=2) 
# lines(density(dist.s$beta.JTU), lty=2, lwd=2)
# 
# #pairwise
# pair.s <- beta.pair(species_binary, index.family = "jaccard")
# # plotting clusters
# dist.s <- species_binary.samp$sampled.values
# plot(hclust(pair.s$beta.jtu, method="average"), hang=-1, main='', sub='', xlab='')
# title(xlab=expression(beta[jtu]), line=0.3)
# plot(hclust(pair.s$beta.jne, method="average"), hang=-1, main='', sub='',, xlab='')
# title(xlab=expression(beta[jne]), line=0.3)

#### Nematode abundances boxplots ----

#plot nematodes
nem_bxp <- ggboxplot(all_data, x = "Habitat", y = "Nematodes per g dry soil", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75)  +
  labs(y = expression("Nematodes (individuals g dry soil "^-1*")"))  + theme(
    #remove x axis label, tickes, labels
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
    legend.title = element_blank()
  ) 

#show the boxplot
show(nem_bxp)
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_nematode-abundances_black-green.svg"), nem_bxp)



#run the GLM
glm <- glm(sqrt(all_data$`Nematodes per g dry soil`) ~ all_data$Habitat*all_data$Vegetation)
summary(glm)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = glm)
plot(simulationOutput)

#two-way ANOVA.  
anova <- aov(sqrt(all_data$`Nematodes per g dry soil`) ~ all_data$Habitat * all_data$Vegetation)
summary(anova)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(sqrt(all_data$`Nematodes per g dry soil`) ~ all_data$Habitat * all_data$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)





#### Nematodes vs water moisture linear modelling ----

#plot nematode abundance vs soil moisture. Do with updated moisture numbers (% fresh soil mass, not graviemtric)

#linear model of the relationship
model <- lm(all_data$`Nematodes per g dry soil` ~ all_data$`Water content (% of wet soil mass)`)
summary(model)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)

#plot the relationship
plot <- ggplot(all_data, aes(x = `Water content (% of wet soil mass)`, y = `Nematodes per g dry soil`)) +
  geom_point(color = "black", size = 2) + # Add points
  geom_smooth(method = "lm", color = "red", se = TRUE) + # Add regression line
  # Line of best fit
  annotate("text", x = 50, y = max(all_data$`Water content (% of wet soil mass)`), 
           label = paste("y =", 
                         round(coef(model)[2], 2), "x", 
                         round(coef(model)[1], 2), 
                         "\nR² =", 
                         round(summary(model)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = "Water Content (%)", y = "Nematodes per g dry soil") +  theme(
    
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
  ) 

show(plot)

ggsave(path = "Figures", paste0(Sys.Date(), "_nematodes-vs-soil-water-content.svg"), plot)




#subset bracken absent and present separately
bracken_absent_data <- all_data[all_data$Vegetation == "Bracken Absent", ]
bracken_present_data <- all_data[all_data$Vegetation == "Bracken Present", ]
#linear models not a good fit for comparsiong nematodes vs soil water.  Try GAM instead

#bracken absent model
model_ba <- lm(bracken_absent_data$`Nematodes per g dry soil` ~ bracken_absent_data$`Water content (% of wet soil mass)`)
summary(model_ba)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = model_ba)
plot(simulationOutput)
#plot the relationship
plot <- ggplot(bracken_absent_data, aes(x = `Water content (% of wet soil mass)`, y = `Nematodes per g dry soil`)) +
  geom_point(color = "black", size = 2) + # Add points
  geom_smooth(method = "lm", color = "red", se = TRUE) + # Add regression line
  # Line of best fit
  annotate("text", x = 50, y = max(bracken_absent_data$`Water content (% of wet soil mass)`), 
           label = paste("y =", 
                         round(coef(model_ba)[2], 2), "x", 
                         round(coef(model_ba)[1], 2), 
                         "\nR² =", 
                         round(summary(model_ba)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = "Water Content (%)", y = "Nematodes per g dry soil") +  theme(
    
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
  ) 

show(plot)




#bracken absent model
model_bp <- lm(bracken_present_data$`Nematodes per g dry soil` ~ bracken_present_data$`Water content (% of wet soil mass)`)
summary(model_bp)
#out model seems to fit well
simulationOutput <- simulateResiduals(fittedModel = model_bp)
plot(simulationOutput)
#plot the relationship
plot <- ggplot(bracken_present_data, aes(x = `Water content (% of wet soil mass)`, y = `Nematodes per g dry soil`)) +
  geom_point(color = "black", size = 2) + # Add points
  geom_smooth(method = "lm", color = "red", se = TRUE) + # Add regression line
  # Line of best fit
  annotate("text", x = 50, y = max(bracken_present_data$`Water content (% of wet soil mass)`), 
           label = paste("y =", 
                         round(coef(model_bp)[2], 2), "x", 
                         round(coef(model_bp)[1], 2), 
                         "\nR² =", 
                         round(summary(model_bp)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = "Water Content (%)", y = "Nematodes per g dry soil") +  theme(
    
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
  ) 

show(plot)





#plot the relationship
plot <- ggplot(all_data, aes(x = `Water content (% of wet soil mass)`, y = `Nematodes per g dry soil`)) +
  geom_point(color = "black", size = 2) + # Add points
  geom_smooth(method = "lm", color = "red", se = TRUE) + # Add regression line
  # Line of best fit
  annotate("text", x = 50, y = max(all_data$`Water content (% of wet soil mass)`), 
           label = paste("y =", 
                         round(coef(model)[2], 2), "x", 
                         round(coef(model)[1], 2), 
                         "\nR² =", 
                         round(summary(model)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = "Water Content (%)", y = "Nematodes per g dry soil") +  theme(
    
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
  ) 

show(plot)


library(Hmisc) #needed for the regression to run, splitting bracken absent from bracken present

p <- ggplot(all_data,aes(`Water content (% of wet soil mass)`, `Nematodes per g dry soil`, colour = Vegetation)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm')  +
  scale_colour_manual(name = "Vegetation", 
                      values = c("black", "limegreen")) + labs(x = "Water Content (%)", y = "Nematodes per g dry soil") +  
  # Line of best fit for bracken present
  annotate("text", x = 50, y = max(all_data$`Water content (% of wet soil mass)` - 100), 
           label = paste("y =", 
                         round(coef(model_bp)[2], 2), "x", 
                         round(coef(model_bp)[1], 2), 
                         "\nR² =", 
                         round(summary(model_bp)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "limegreen") +
  # Line of best fit for bracken absent
  annotate("text", x = 50, y = max(all_data$`Water content (% of wet soil mass)`+ 30), 
           label = paste("y =", 
                         round(coef(model_ba)[2], 2), "x", 
                         round(coef(model_ba)[1], 2), 
                         "\nR² =", 
                         round(summary(model_ba)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") +
  
  theme(
                        
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
                        legend.title = element_blank()
                      ) 

show(p)

ggsave(path = "Figures", paste0(Sys.Date(), "_nematodes-vs-soil-water-content_bracken-absent-vs-present.svg"), p, width = 7, heigh = 5, dpi = 300)




#### Effect of soil moisture on number of morphotypes and individuals extracted ----

#dataframe of soil moisture
moist <- all_data[,-(5:11)]
moist <- moist[,-(6:429)]


#read in the week 1 extracts
filename <- "6) Week-1-morphotypes.csv"
#csv containing only mites and springtails
#filename <- "mites-springtails-only week-1-mesofauna-corrected.csv"
wk1 <- readr::read_csv(
  here::here("Data", filename), show_col_types = FALSE
) 
#order samples by ID alphabetically
wk1 <- arrange(wk1, wk1["Sample ID"])
wk1 <- as.data.frame(wk1)
#remove all empty rows
wk1 <- wk1[1:30,]
#replace null (empty excell cell) with "0"
wk1[is.na(wk1)] <- 0
#replace row index with sample names
rownames(wk1) <- wk1[,1]
#make sure our variables are coded as factors
wk1$Habitat <- factor(wk1$Habitat, levels = c("Grassland", "Heathland", "Woodland"), labels = c("Grassland", "Heathland", "Woodland"))
#just the morphospecies counts
wk1spe <- wk1[,-(1:4)]
wk1spe <- as.matrix(wk1spe)
#number of morphospecies extracted after 7 days
#species richness
wk1$Richness_1 <- apply(wk1spe[,]>0,1,sum)
#number of individuals extracted after 7 days
wk1$Individuals_1 <- rowSums(wk1spe)


#order samples by ID alphabetically
wk2 <- all_data[, 1:4]
wk2 <- as.data.frame(wk2)
#make sure our variables are coded as factors
wk2$Habitat <- factor(wk2$Habitat, levels = c("Grassland", "Heathland", "Woodland"), labels = c("Grassland", "Heathland", "Woodland"))


#just the morphospecies counts of week 1 + 2 extracts
wk2spe <- all_data[,51:436]
rownames(wk2spe) <- all_data[,1]
#number of morphospecies extracted after 7 days
#species richness
wk2$Richness_2 <-apply(wk2spe[,]>0,1,sum)
#number of individuals extracted after 7 days
wk2$Individuals_2 <- rowSums(wk2spe)

#create dataframe containing moisture, and richness/individuals from 7 and 14 day extracts
moist <- merge(moist, wk1[, c("Sample ID", "Richness_1", "Individuals_1")], by = "Sample ID", all.x = TRUE)
#create dataframe containing moisture, and richness/individuals from 7 and 14 day extracts
moist <- merge(moist, wk2[, c("Sample ID", "Richness_2", "Individuals_2")], by = "Sample ID", all.x = TRUE)

#reshape the data into long format
# Example of reshaping (assuming your columns are named 'Richness_1' and 'Richness_2' for species richness at 1 and 2 weeks)
moist_long <- moist %>%
  pivot_longer(cols = c("Richness_1", "Richness_2", "Individuals_1", "Individuals_2"),
               names_to = c(".value", "Week"),
               names_pattern = "(.*)_(.*)")

#number of morphospecies and individuals gained from an additional week of extraction
moist$RichnessGained <- moist$Richness_2 - moist$Richness_1
moist$IndividualsGained <- moist$Individuals_2 - moist$Individuals_1


#is there a difference ebtween soil moisture and addition morphospecies found in week 2?
anv <- aov(moist$RichnessGained ~ moist$`Water content (% of wet soil mass)` + moist$Habitat/moist$Vegetation)
summary(anv)

#is there a relationship between soil moisture and addition morphospecies found in week 2?
anv <- aov(moist$IndividualsGained ~ moist$`Water content (% of wet soil mass)` + moist$Habitat/moist$Vegetation)
summary(anv)


#this is the way rob said to analyse this, but no fucking idea how to interpret the results
library(lme4) #for linear modelling
library(Matrix) #for linear modelling to work
# Mixed effects model for Species Richness
richness_model <- lmer(Richness ~ `Water content (% of wet soil mass)` * Week + (1 | `Sample ID`), data = moist_long)
# Summary of the model
summary(richness_model)



# Mixed effects model for Species Richness
individual_model <- lmer(Individuals ~ `Water content (% of wet soil mass)` * Week + (1 | `Sample ID`), data = moist_long)
# Summary of the model
summary(individual_model)


#### run PCA, to see what is driving mesofauna diversity ----

#generate dataframe ew shall run PCA on.  Ensure veg and morphotype alpha diversity indices
#have been computed and added to the dataframe using the relevant tabs

library('FactoMineR') #includes functions needed to vizualize outputs of the PCA

#remove morphospecies count data
d <- all_data[, -(51:436)]
#remove non-numerical data
numerical_d <- d[, 7:9]
numerical_d <- cbind(numerical_d, d[, 12:58])
#remove column as it is all 0s
numerical_d <- numerical_d[, -26]
#remove non-numerical data
numerical_d <- na.omit(numerical_d)
# Remove columns that are entirely zeros
numerical_d <- numerical_d[, colSums(numerical_d != 0) > 0]
#remove veg diversity indices
numerical_d <- numerical_d[, -c(40:43)]
#normalize the data
data_normalized <- scale(numerical_d)
#compute the PCA
pca_result <- prcomp(data_normalized)
# View the results
summary(pca_result)  # Summary of the PCA (explained variance, etc.)

# Principal components
pca_result$x  # The scores (the projections of the original data onto the PCs)

# The loadings (eigenvectors).
pca_result$rotation  # The loadings (coefficients of the principal components)

# Variance explained by each principal component
pca_result$sdev^2 / sum(pca_result$sdev^2)  # Proportion of variance explained


#needed for fviz_eig() function
library(factoextra)
#generate scree plot using fviz_eig() function.  This shows the eigenvalues from highest to lowest, i.e. from the components which explain the most variance to the components which explain the least
fviz_eig(pca_result, addlabels = TRUE)


# Plotting the PCA results (optional)
# First PC vs. Second PC
plot(pca_result$x[, 1], pca_result$x[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA - PC1 vs PC2")


# Biplot to visualize the PCA results.  We are vizualising the similarities and differences between samples, and shows the impact of each attribute on each of the principal components.  Variables that are grouped together are positively correlated with one another.  The further the distance between the variable and the origin, the better represented the variaible is.  Variables that are negatively correlated are displayed to the opposite side of the biplot's origin.
biplot(pca_result)

#now determine the variable's contribution to principal components.  This representation is called the Cos2, and corresponds to the square cosine.  A low value means the variable is not perfectly represented by that component, whilst a high value means a good representation of the variable on that component.
fviz_cos2(pca_result, choice = "var", axes = 1:2)

dev.new()
#combine biplot and attribute importance.  Attributes with similar cos2 scores will have similar colours.
pca <- fviz_pca_var(pca_result, col.var = "cos2", 
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
#save the PCA
ggsave(path = "Figures", paste0(Sys.Date(), '_PCA.svg'), width = 14, height = 10, pca)
