args <- commandArgs(trailingOnly = TRUE)

###Arguments
matrix_file_path <- args[1]
dist_method <- args[2] ## manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, mahalanobis, chisq, chord, hellinger, aitchison, or robust.aitchison
num_centers <- as.integer(args[3])
confidence <- as.numeric(args[4])

pkgs <- c('vegan', 'ggplot2', 'tidyr', 'dplyr', 'ellipse', 'tibble')
for (pkg in pkgs) {
  suppressMessages({
    if (!require(pkg, character.only=TRUE)) {
      install.packages(pkg)
      suppressMessages(library(pkg, character.only=TRUE))
    }
  })
}

allowed_methods <- c('manhattan', 'euclidean', 'canberra', 'clark', 'bray', 
                     'kulczynski', 'jaccard', 'gower', 'altGower', 'morisita', 
                     'horn', 'mountford', 'raup', 'binomial', 'chao', 'cao', 
                     'mahalanobis', 'chisq', 'chord', 'hellinger', 'aitchison', 
                     'robust.aitchison')

# Check if the chosen method is allowed
if (!dist_method %in% allowed_methods) {
  stop(paste("Invalid distance method:", dist_method,
             "\nPlease choose from the following methods:\n", 
             paste(allowed_methods, collapse = ", ")),
       call. = FALSE)
}

print(paste("Data File:", args[1]))
print(paste("Method:", args[2]))
print(paste("Number of Centers:", args[3]))
print(paste("Confidence:", args[4]))

##Read in data
data <- read.csv(matrix_file_path)

##Transpose data
df <- setNames(data.frame(t(data[,-1])), data[,1])

##Set name of column one vars
df <- df %>% rownames_to_column(var="vars")

dist_matrix <- vegdist(df[,-1], method=dist_method)

mva <- cmdscale(dist_matrix, eig=TRUE, k=2)

eigvals <- mva$eig
percent_var <- eigvals / sum(eigvals) * 100

clusters <- kmeans(mva$points, centers = num_centers)

df <- data.frame(vars=df$vars, mva$points, cluster=factor(clusters$cluster))

percent_var1 <- percent_var[1]
percent_var2 <- percent_var[2]

# Function to add ellipses to the plot
add_ellipses <- function(plot, df, confidence) {
  tryCatch({
    plot + 
      stat_ellipse(aes(fill=cluster, group=cluster), type = 'norm', level=confidence, 
                   geom = 'polygon', alpha = I(0.05), show.legend=F) +
      stat_ellipse(aes(group=cluster), type = 'norm', level=confidence, linetype = 4, 
                   color="black", linewidth=1, 
                   geom = 'path', show.legend=FALSE)
  }, error = function(e) {
    message("Failed to add confidence ellipses: ", e$message)
    return(plot)
  })
}

# Initial plot without ellipses
plot <- ggplot(df, aes(x=X1, y=X2)) +
  theme_bw() +
  geom_point(aes(color=cluster, alpha=I(0.75)), size=7, show.legend=FALSE) +
  geom_text(aes(label=vars), vjust=-1, hjust=1, size=4) +
  labs(x = paste0("PCo1 (", round(percent_var1, 2), "%)"), 
       y = paste0("PCo2 (", round(percent_var2, 2), "%)")) +
  theme(
    axis.title.x = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    axis.text.x = element_text(size=25),
    axis.text.y = element_text(size=25),
    legend.text = element_text(size=25),
    legend.key.size = unit(1, 'cm'),
    legend.title = element_text(size=25))

# Attempt to add ellipses
plot <- add_ellipses(plot, df, confidence)

# Save the plot
ggsave(filename="MVA.png", plot=plot, width=10, height=10, dpi=300)
ggsave(filename="MVA.pdf", plot=plot, width=10, height=10)