# Analysis Synapse results 

# ANOVA 

# Gephyrin and vGat:----- 

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Create the dataset
data <- data.frame(
  animal_number = c("SL0171", "SL2389", "SL3437", "SL6617"),
  Control = c(521, 490, 604, 326),
  Perilesional = c(127, 405, 417, 497),
  Proximal = c(319, 336, 404, 501),
  Distal = c(301, 446, 328, 321)
)

# Reshape the data to long format
data_long <- pivot_longer(data, 
                          cols = Control:Distal, 
                          names_to = "condition", 
                          values_to = "value")

# Perform one-way ANOVA
anova_result <- aov(value ~ condition, data = data_long)

# Summary of ANOVA
summary(anova_result)

# Post-hoc test (Tukey's HSD) for pairwise comparisons
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Visualization: Boxplot for Gephyrin and vGAT
gephyrin_vgat_plot <- ggplot(data_long, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(width = 0.3, outlier.shape = 16, outlier.size = 2, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7, shape = 21, size = 3, color = "black") +
  scale_fill_manual(
    values = c(
      "Control" = "grey",
      "Perilesional" = "#4DAF4A",
      "Proximal" = "#E41A1C",
      "Distal" = "#377EB8"
    )
  ) +
  labs(
    title = "Colocalisation Count: Gephyrin and vGAT by Condition",
    x = "Condition",
    y = "Colocalisation Count per ROI"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

# Save the plot
ggsave(
  "C:/Users/yperstat/Documents/Studium/IDB_Master/Master_thesis/Neuron_subtype_Stainings/gephyrin_vgat_results_by_condition.png",
  gephyrin_vgat_plot, width = 8, height = 6, dpi = 600
)


# Homer and vglut: ----
# Load necessary libraries
library(lme4)
library(lmerTest)
library(tidyr)
library(dplyr)
library(emmeans)
library(ggplot2)

# Create the dataset
data <- data.frame(
  animal_number = c("SL0171", "SL2389", "SL3437", "SL6617"),
  Control = c(307, 252, NA, 102),
  Perilesional = c(336, 382, 429, NA),
  Proximal = c(258, 412, 83, 200),
  Distal = c(413, 387, 42, 123)
)

# Reshape the data to long format
data_long <- pivot_longer(data, 
                          cols = Control:Distal, 
                          names_to = "condition", 
                          values_to = "value")

# Remove rows with missing values
data_long <- na.omit(data_long)

# Fit the Linear Mixed-Effects Model
lmm <- lmer(value ~ condition + (1 | animal_number), data = data_long)

# Summary of the model
summary(lmm)

# Perform pairwise comparisons
pairwise_results <- emmeans(lmm, pairwise ~ condition, adjust = "tukey")
print(pairwise_results)

# Visualization: Boxplot for Homer1 and vGLUT
homer_vglut_plot <- ggplot(data_long, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(width = 0.3, outlier.shape = 16, outlier.size = 2, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7, shape = 21, size = 3, color = "black") +
  scale_fill_manual(
    values = c(
      "Control" = "grey",
      "Perilesional" = "#4DAF4A",
      "Proximal" = "#E41A1C",
      "Distal" = "#377EB8"
    )
  ) +
  labs(
    title = "Colocalisation Count: Homer1 and vGLUT by Condition",
    x = "Condition",
    y = "Colocalisation Count per ROI"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

# Save the plot
ggsave(
  "C:/Users/yperstat/Documents/Studium/IDB_Master/Master_thesis/Neuron_subtype_Stainings/homer_vglut_results_by_condition.png",
  homer_vglut_plot, width = 8, height = 6, dpi = 600
)

# Bassoon: ----
# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

# Create the dataset
data <- data.frame(
  animal_ID = c("SL6940", "SL6940", "SL6940", "SL6940", "SL0171", "SL2390", "SL2390", "SL2390", "SL2390", "SL2390", 
                "SL6617", "SL6617", "SL6617", "SL6617", "SL6617"),
  region = c("Proximal", "middle", "Distal", "Perilesional", "Control", "Control", "Proximal", "middle", 
             "Distal", "Perilesional", "Control", "Proximal", "middle", "Distal", "Perilesional"),
  count = c(607, 536, 552, 757, 267, 586, 784, 809, 890, 740, 874, 643, 684, 640, 619)
)

# Remove the "middle" region
data_filtered <- data %>% filter(region != "middle")

# Fit the Linear Mixed-Effects Model
lmm <- lmer(count ~ region + (1 | animal_ID), data = data_filtered)

# Summary of the model
summary(lmm)

# Post-hoc pairwise comparisons
pairwise_results <- emmeans(lmm, pairwise ~ region, adjust = "tukey")
print(pairwise_results)

# Visualization: Bassoon Results by Region
bassoon_plot <- ggplot(data_filtered, aes(x = region, y = count, fill = region)) +
  geom_boxplot(width = 0.3, outlier.shape = 16, outlier.size = 2, alpha = 0.5, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7, shape = 21, size = 3, color = "black") +
  scale_fill_manual(
    values = c(
      "Control" = "grey",
      "Perilesional" = "#4DAF4A",
      "Proximal" = "#E41A1C",
      "Distal" = "#377EB8"
    )
  ) +
  labs(
    title = "Bassoon Results by Region",
    x = "Condition",
    y = "Particle Count per ROI"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

# Save the plot
ggsave(
  "C:/Users/yperstat/Documents/Studium/IDB_Master/Master_thesis/Neuron_subtype_Stainings/bassoon_results_by_region.png",
  bassoon_plot, width = 8, height = 6, dpi = 600
)
