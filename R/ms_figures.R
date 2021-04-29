########################
# Main figures
########################

# Figure 4A
moduleColors <- read.csv("../results/files/moduleColors.csv")
df <- data.frame(dc = moduleColors) %>% dplyr::count(dc)
mergedModuleSizes <- ggplot(df , aes(x = reorder(dc, n),
                                     y = n,
                                     fill = dc)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = as.character(df$dc)) +
  ylim(c(0, 7500)) +
  theme_minimal() +
  # Remove panel borders and grid lines
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Number of genes") +
  xlab("") +
  geom_text(aes(label = n),
            color = "black",
            size = 10/.pt,
            family = "sans",
            fontface = "bold",
            hjust = -0.1,
            vjust = 0.5) +
  theme(legend.position = 'None',
        axis.title = element_text(size = 10,
                                  family = "sans",
                                  face = "bold"),
        axis.text = element_text(size = 10,
                                 family = "sans",
                                 face = "bold"),
        aspect.ratio = 1) # square plot area
# save svg version
svg("mergedModuleSizes.svg", width = 10, height = 10)
mergedModuleSizes
dev.off()

# Figure 4B and C - Preparation
source("../R/own_functions.R")
g4reac <- read.csv("../results/files/g4reac.csv")
data <- module_enrich_pathway("blue")
# Set path for subfolders that contain figures for each module
my_path <- paste(getwd(), "/", figures, sep = "")
my_path <- getwd()
# Create folders for these modules
create_color_folder(my_path, "blue")
create_color_folder(my_path, "turquoise")

# Figure 4B
complete_visualization(color = "blue",
                       nr_categories = 15,
                       type = "dotplot",
                       path = my_path)

# Figure 4C
complete_visualization(color = "turquoise",
                       nr_categories = 15,
                       type = "dotplot",
                       path = my_path)


########################
# Supplementary Figures
########################
