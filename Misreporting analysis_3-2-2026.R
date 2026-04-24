## Jing Li ##
## Misreporting_2019 data ##
## 3/2/2026 ##

setwd("C:/Users/sammi/Dropbox/Jing Li/Dissertation/Dissertation Writeup/Misreporting chapter/Educational and Psychological Measurement/Analysis")
# Hanning (Smoothing) of Person Response Functions: https://www.rasch.org/rmt/rmt264b.htm

## read the data
library(readxl)
d1 = as.data.frame(read_excel("2019data.xlsx"))
head(d1)
table(d1$hh_fscat)
table(d1$adult_fscat)
table(d1$child_fscat)  # N=2582

## covert polytomous to dichotomous: 0 =Low,,,0; 1 =High,,,1+2+3
d1_poly= d1[,c(1,3, 4, 5, 6:23)]
d1_dicho <- as.data.frame(lapply(d1_poly[ , 5:22], function(x) ifelse(x %in% c(1, 2, 3), 1, x)))
dim(d1_dicho) # N=2582 adult with children; dichotomous 
head(d1_dicho)

## run a joint maximum likelihood estimation of the Rasch model
library(TAM)
mod1= tam.jml(resp=d1_dicho) # calibrate dif to average 0
summary(mod1)

#############################
# add variables to the data #
#############################
person_fit <- as.data.frame(round(tam.personfit(mod1),2)) # person infit and outfit
person_abilities <- as.data.frame(mod1$theta)
person_dat <- as.data.frame(cbind(d1$hhid,
  d1_dicho[, c("I1", "I2", "I3", "I11", "I12", "I5", "I4", "I6", "I7","I13", "I8", 
              "I9", "I10", "I14", "I15", "I16", "I17", "I18"),], round(person_abilities,2)
  ,person_fit[, c("infitPerson", "outfitPerson")]
))

colnames(person_dat)[colnames(person_dat) == "infitPerson"] <- "Infit"
colnames(person_dat)[colnames(person_dat) == "outfitPerson"] <- "Outfit"
colnames(person_dat)[colnames(person_dat) == "d1$hhid"] <- "ID"
colnames(person_dat)[colnames(person_dat) == "mod1$theta"] <- "ability"
head(person_dat)


# add fit category to the data (using INfit categorize them; using outfit categorize them)
person_dat$Infitcat <- ifelse(person_dat$Infit >= 2, "D",
                              ifelse((person_dat$Infit >= 1.5 & person_dat$Infit < 2), "C",
                                     ifelse((person_dat$Infit >= 0.5 & person_dat$Infit < 1.5), "A",
                                          ifelse(person_dat$Infit < 0.5, "B", NA))))

person_dat$Outfitcat <- ifelse(person_dat$Outfit >= 2, "D",
                              ifelse((person_dat$Outfit >= 1.5 & person_dat$Outfit < 2), "C",
                                     ifelse((person_dat$Outfit >= 0.5 & person_dat$Outfit < 1.5), "A",
                                            ifelse(person_dat$Outfit < 0.5, "B", NA))))
head(person_dat)


person_dat_misfit <- subset(person_dat,
                            Infitcat %in% c("C","D") | Outfitcat %in% c("C","D")) # N=114 (C or D)

person_dat_misfit$pecent= rowMeans(person_dat_misfit[, 2:19], na.rm = TRUE) # add in percentage getting correct
person_dat_misfit_50= subset(person_dat_misfit, person_dat_misfit$pecent==.5) # 50% and misfit; N=8

#####################
#   clustering      #
#####################

# ----- Clustering ----- 
library(cluster)
library(dplyr)
library(ggplot2)
library(ggrepel)

## ----- Build response matrix (people x items) ------
resp_mis <- as.matrix(person_dat_misfit[, 2:19])
resp_mis <- apply(resp_mis, 2, function(x) as.numeric(as.character(x)))
resp_mis <- as.matrix(resp_mis)
rownames(resp_mis) <- person_dat_misfit$ID

## ----- Distance matrix (scaled Manhattan = proportion mismatches) -----
d_mis <- dist(resp_mis, method = "manhattan") / ncol(resp_mis)

## ----- Distance matrix (scaled Manhattan = proportion mismatches) -----
d_mis <- dist(resp_mis, method = "manhattan") / ncol(resp_mis)


k_grid <- 2:min(8, nrow(resp_mis)-1) # try 2-8 clusters

## ----- pam() pairwise distances between people -----
# PAM us medoids: actual data points that best represent each cluster
sil <- sapply(k_grid, function(k){
  cl <- pam(d_mis, k=k)$clustering # k-medoids; 
  mean(silhouette(cl, d_mis)[, 3]) # compute silhouette scores; Average silhouette for this k
})

k_best <- k_grid[which.max(sil)]
pam_fit <- pam(d_mis, k=k_best)
person_dat_misfit$cluster_misfit <- pam_fit$clustering
table(person_dat_misfit$cluster_misfit)

d1_misfit_5clus <- person_dat_misfit[person_dat_misfit$cluster_misfit == 1, ] # response

## ----- compare K=2-5 ------
compare_df <- data.frame(
  K = k_grid,
  Avg_Silhouette = sil
)

compare_2_5 <- subset(compare_df, K %in% 2:5)
print(compare_2_5)

## ----- compare K = 2:5 at the individual-cluster level -----
cluster_perf_list <- lapply(2:5, function(k) {
  
  # fit PAM
  pam_k <- pam(d_mis, k = k)
  
  # silhouette for each person
  sil_k <- silhouette(pam_k$clustering, d_mis)
  sil_k <- as.data.frame(sil_k)
  
  # summarize by cluster
  cluster_summary <- aggregate(
    sil_width ~ cluster,
    data = sil_k,
    FUN = function(x) c(
      n = length(x),
      mean = mean(x),
      sd = sd(x),
      min = min(x),
      max = max(x),
      neg_prop = mean(x < 0)
    )
  )
  
  # make it clean
  cluster_summary2 <- data.frame(
    K = k,
    Cluster = cluster_summary$cluster,
    n = cluster_summary$sil_width[, "n"],
    Mean_Silhouette = cluster_summary$sil_width[, "mean"],
    SD_Silhouette = cluster_summary$sil_width[, "sd"],
    Min_Silhouette = cluster_summary$sil_width[, "min"],
    Max_Silhouette = cluster_summary$sil_width[, "max"],
    Prop_Negative = cluster_summary$sil_width[, "neg_prop"]
  )
  
  return(cluster_summary2)
})

cluster_perf_df <- do.call(rbind, cluster_perf_list)
print(cluster_perf_df)

## ----- 1) Fit PAM with K = 5 (d_mis must already exist as a dist object)-----
pam_fit5 <- pam(d_mis, k = 5)

# (optional) store cluster labels
person_dat_misfit$cluster_misfit5 <- pam_fit5$clustering
table(person_dat_misfit$cluster_misfit5)

## ----- 2) 2D embedding from distance matrix (MDS)  (same xy is fine, but recompute for clarity) -----
xy <- cmdscale(d_mis, k = 2)

plot_df5 <- data.frame(
  ID      = as.character(person_dat_misfit$ID),
  Dim1    = xy[, 1],
  Dim2    = xy[, 2],
  Cluster = factor(pam_fit5$clustering)
)

# (optional but recommended) ensure IDs align with distance labels
labs <- attr(d_mis, "Labels")
if (!is.null(labs) && length(labs) == nrow(plot_df5)) plot_df5$ID <- labs

# 3) Medoids (actual observed individuals)
med_df5 <- plot_df5[pam_fit5$id.med, , drop = FALSE]

# 4) Build segment data (connect each point to its cluster medoid)
med_lookup5 <- med_df5[, c("Cluster", "Dim1", "Dim2")]
names(med_lookup5) <- c("Cluster", "mx", "my")

seg_df5 <- plot_df5
seg_df5$mx <- med_lookup5$mx[match(seg_df5$Cluster, med_lookup5$Cluster)]
seg_df5$my <- med_lookup5$my[match(seg_df5$Cluster, med_lookup5$Cluster)]
stopifnot(!anyNA(seg_df5$mx), !anyNA(seg_df5$my))

# 5) Plot: points + connection lines + ellipses + medoids
p5 <- ggplot(plot_df5, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_segment(
    data = seg_df5,
    aes(x = mx, y = my, xend = Dim1, yend = Dim2, color = Cluster),
    inherit.aes = FALSE,
    alpha = 0.45,          # a bit lighter since 5 clusters = more lines
    linewidth = 0.55,
    lineend = "round",
    show.legend = FALSE
  ) +
  stat_ellipse(aes(group = Cluster, color = Cluster),
               linewidth = 1, level = 0.95) +
  geom_point(
    data = med_df5,
    aes(Dim1, Dim2),
    inherit.aes = FALSE,
    shape = 4,
    size = 6,
    stroke = 1.8,
    color = "black"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Misfit Person Clustering (PAM, K = 5)",
    x = "Dimension 1 (MDS)",
    y = "Dimension 2 (MDS)",
    color = "Cluster"
  )

print(p5)

## Quality of 5 clusters 
library(cluster)
cl <- person_dat_misfit$cluster_misfit
sil_obj <- silhouette(cl, d_mis)

# overall quality
mean(sil_obj[, 3]) # global clustering quality
# > 0.50 → good separation
# 0.25–0.50 → moderate / usable but overlapping
# < 0.25 → weak structure (clusters overlap a lot)

tapply(sil_obj[, 3], cl, mean) # Mean silhouette within each cluster； Helps you find bad clusters

# visualize silhouette
plot(sil_obj, border = NA)
abline(v = 0.25, col="red", lty=2) # The red vertical line at 0.25 is your chosen “minimum acceptable” reference.

# Mean silhouette width (you already computed it): 
# > 0.50: good separation
# 0.25-0.50: moderate
# < 0.25: weak / likely overlapping

##################
#   4 clusters   #
##################
# ----- PAM clustering (K = 4) + visualization with IDs -----
pam_fit <- pam(d_mis, k = 4)
person_dat_misfit$cluster_misfit <- pam_fit$clustering

# 2D embedding from SAME distance (MDS)
xy <- cmdscale(d_mis, k = 2)

# Build plot_df in the SAME order as d_mis / pam_fit$clustering
plot_df <- data.frame(
  ID      = as.character(person_dat_misfit$ID),
  Dim1    = xy[, 1],
  Dim2    = xy[, 2],
  Cluster = factor(person_dat_misfit$cluster_misfit)
)

# Medoids (indices are in the same order as pam_fit$clustering)
med_idx <- pam_fit$id.med
med_df  <- plot_df[med_idx, , drop = FALSE]

# ---- IMPORTANT: create med_lookup (you missed this line) ----
med_lookup <- med_df[, c("Cluster", "Dim1", "Dim2")]
names(med_lookup) <- c("Cluster", "mx", "my")

# Safer than merge(): keep row order and map correct medoid coords by cluster
seg_df <- plot_df
seg_df$mx <- med_lookup$mx[match(seg_df$Cluster, med_lookup$Cluster)]
seg_df$my <- med_lookup$my[match(seg_df$Cluster, med_lookup$Cluster)]
# 4) Plot: points + ellipses + medoids + ID labels
# p <- ggplot(plot_df, aes(Dim1, Dim2, color = Cluster)) +
#   geom_point(size = 2.8, alpha = 0.85) +
#   stat_ellipse(aes(group = Cluster, color = Cluster),
#                linewidth = 1, level = 0.95) +
#   geom_point(
#     data = med_df,
#     aes(Dim1, Dim2),
#     inherit.aes = FALSE,
#     shape = 4,
#     size = 6,
#     stroke = 1.8,
#     color = "black"
#   ) +
#   theme_classic(base_size = 13) +
#   labs(
#     title = "Misfit Person Clustering (PAM, K = 4)",
#     x = "Dimension 1 (MDS)",
#     y = "Dimension 2 (MDS)",
#     color = "Cluster"
#   )
# print(p)

p <- ggplot(plot_df, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_segment(
    data = seg_df,
    aes(x = mx, y = my, xend = Dim1, yend = Dim2, color = Cluster),
    inherit.aes = FALSE,
    alpha = 0.40,
    linewidth = 0.45,
    lineend = "round",
    show.legend = FALSE
  ) +
  stat_ellipse(aes(group = Cluster, color = Cluster),
               linewidth = 1, level = 0.95) +
  geom_point(
    data = med_df,
    aes(Dim1, Dim2),
    inherit.aes = FALSE,
    shape = 4, size = 6, stroke = 1.8, color = "black"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Misfit Person Clustering (PAM, K = 4)",
    x = "Dimension 1 (MDS)",
    y = "Dimension 2 (MDS)",
    color = "Cluster"
  )

print(p)

# K=4 quality
sil4 <- silhouette(pam_fit$clustering, d_mis)
mean(sil4[,3])                       # compare to 0.3267 (K=5)
tapply(sil4[,3], pam_fit$clustering, mean)  # cluster-wise
plot(sil4, border = NA)
abline(v = 0.25, col="red", lty=2)

# another K=4 cluster accuracy
sil_by_cluster <- by(sil4[, 3], sil4[, 1], function(x) {
  c(
    n = length(x),
    mean_sil = mean(x),
    median_sil = median(x),
    min_sil = min(x),
    pct_neg = mean(x < 0),        # "likely mis-assigned" points
    pct_gt_025 = mean(x > 0.25)   # "reasonably good" points (your red line)
  )
})

sil_summary <- do.call(rbind, sil_by_cluster)
sil_summary

# metroids
pam_fit$medoids
medoid_points <- plot_df[pam_fit$medoids,]
medoid_points
medoid_ids <- c(96156, 47615, 71340, 5928) 
# Filter rows from person_dat_misfit where ID matches 
medoid_data <- person_dat_misfit %>% filter(ID %in% medoid_ids) 
medoid_data

# Loop over ONLY the medoid individuals
ordered_items <- mod1$item1[order(mod1$item1$xsi), c(1, 3)]# dif from low to high
n_persons <- nrow(medoid_data)

for (i in 1:n_persons) {
  # Extract responses (columns 2 to 19)
  responses <- as.numeric(medoid_data[i, 2:19])
  # Person ID
  person_id <- medoid_data[i, "ID"]
  # Title showing ID + raw response string
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  # Initialize smoothed responses
  smoothed_responses <- responses
  n_items <- length(responses)
  # Hanning smoothing (9 iterations)
  for (iteration in 1:9) {
    temp_responses <- smoothed_responses
    for (j in 2:(n_items - 1)) {
      smoothed_responses[j] <- (temp_responses[j - 1] +
                                  2 * temp_responses[j] +
                                  temp_responses[j + 1]) / 4
    }
  }
  
  # Plot
  plot(ordered_items$xsi, smoothed_responses,
       type = "n",
       ylim = c(0, 1),
       xlim = c(0, 10),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = "black", lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = "black", pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = "red", pch = 19, cex = 1)
}

## plot the non-parametric PRF based on 4-cluster
pam_fit <- pam(d_mis, k = 4)
person_dat_misfit$cluster_misfit <- pam_fit$clustering

clus1data <- person_dat_misfit[person_dat_misfit$cluster_misfit == 1, ]
clus2data <- person_dat_misfit[person_dat_misfit$cluster_misfit == 2, ]
clus3data <- person_dat_misfit[person_dat_misfit$cluster_misfit == 3, ]
clus4data <- person_dat_misfit[person_dat_misfit$cluster_misfit == 4, ]

### ----- plot cluster 1 -----
# response columns (keep your original 2:19)
resp_cols <- 2:19

# optional: if ordered_items$item matches your column names, reorder columns to match difficulty order
if ("item" %in% names(ordered_items)) {
  item_names <- as.character(ordered_items$item)
  if (all(item_names %in% names(clus1data))) {
    resp_cols <- match(item_names, names(clus1data))
  }
}

# smoothing function (same logic as your loop)
hanning_smooth <- function(x, n_iter = 9) {
  x <- as.numeric(x)
  n <- length(x)
  sm <- x
  for (iter in seq_len(n_iter)) {
    tmp <- sm
    for (j in 2:(n - 1)) {
      sm[j] <- (tmp[j - 1] + 2 * tmp[j] + tmp[j + 1]) / 4
    }
  }
  sm
}

# 3x2 panels (will cycle pages automatically in the plot window)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))

n_persons <- nrow(clus1data)

for (i in seq_len(n_persons)) {
  
  # extract responses + ID
  responses <- as.numeric(clus1data[i, resp_cols])
  person_id <- as.character(clus1data[i, "ID"])
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  # smooth
  smoothed_responses <- hanning_smooth(responses, n_iter = 9)
  
  # plot
  plot(ordered_items$xsi, smoothed_responses,
       type = "n",
       ylim = c(0, 1),
       xlim = c(0, 10),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = "black", lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = "black", pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = "red", pch = 19, cex = 1)
}

# reset (optional)
par(mfrow = c(1, 1))


### ------ plot cluster 2 -------
# response columns (keep your original 2:19)
resp_cols <- 2:19

# optional: if ordered_items$item matches your column names, reorder columns to match difficulty order
if ("item" %in% names(ordered_items)) {
  item_names <- as.character(ordered_items$item)
  if (all(item_names %in% names(clus2data))) {
    resp_cols <- match(item_names, names(clus2data))
  }
}

# smoothing function (same logic as your loop)
hanning_smooth <- function(x, n_iter = 9) {
  x <- as.numeric(x)
  n <- length(x)
  sm <- x
  for (iter in seq_len(n_iter)) {
    tmp <- sm
    for (j in 2:(n - 1)) {
      sm[j] <- (tmp[j - 1] + 2 * tmp[j] + tmp[j + 1]) / 4
    }
  }
  sm
}

# 3x2 panels (will cycle pages automatically in the plot window)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))

n_persons <- nrow(clus2data)

for (i in seq_len(n_persons)) {
  
  # extract responses + ID
  responses <- as.numeric(clus2data[i, resp_cols])
  person_id <- as.character(clus2data[i, "ID"])
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  # smooth
  smoothed_responses <- hanning_smooth(responses, n_iter = 9)
  
  # plot
  plot(ordered_items$xsi, smoothed_responses,
       type = "n",
       ylim = c(0, 1),
       xlim = c(0, 10),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = "black", lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = "black", pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = "red", pch = 19, cex = 1)
}

# reset (optional)
par(mfrow = c(1, 1))


### -------- plot cluster 3 --------
# response columns (keep your original 2:19)
resp_cols <- 2:19

# optional: if ordered_items$item matches your column names, reorder columns to match difficulty order
if ("item" %in% names(ordered_items)) {
  item_names <- as.character(ordered_items$item)
  if (all(item_names %in% names(clus3data))) {
    resp_cols <- match(item_names, names(clus3data))
  }
}

# smoothing function (same logic as your loop)
hanning_smooth <- function(x, n_iter = 9) {
  x <- as.numeric(x)
  n <- length(x)
  sm <- x
  for (iter in seq_len(n_iter)) {
    tmp <- sm
    for (j in 2:(n - 1)) {
      sm[j] <- (tmp[j - 1] + 2 * tmp[j] + tmp[j + 1]) / 4
    }
  }
  sm
}

# 3x2 panels (will cycle pages automatically in the plot window)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))

n_persons <- nrow(clus3data)

for (i in seq_len(n_persons)) {
  
  # extract responses + ID
  responses <- as.numeric(clus3data[i, resp_cols])
  person_id <- as.character(clus3data[i, "ID"])
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  # smooth
  smoothed_responses <- hanning_smooth(responses, n_iter = 9)
  
  # plot
  plot(ordered_items$xsi, smoothed_responses,
       type = "n",
       ylim = c(0, 1),
       xlim = c(0, 10),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = "black", lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = "black", pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = "red", pch = 19, cex = 1)
}

# reset (optional)
par(mfrow = c(1, 1))

### ------ plot cluster 4 -------
# response columns (keep your original 2:19)
resp_cols <- 2:19

# optional: if ordered_items$item matches your column names, reorder columns to match difficulty order
if ("item" %in% names(ordered_items)) {
  item_names <- as.character(ordered_items$item)
  if (all(item_names %in% names(clus4data))) {
    resp_cols <- match(item_names, names(clus4data))
  }
}

# smoothing function (same logic as your loop)
hanning_smooth <- function(x, n_iter = 9) {
  x <- as.numeric(x)
  n <- length(x)
  sm <- x
  for (iter in seq_len(n_iter)) {
    tmp <- sm
    for (j in 2:(n - 1)) {
      sm[j] <- (tmp[j - 1] + 2 * tmp[j] + tmp[j + 1]) / 4
    }
  }
  sm
}

# 3x2 panels (will cycle pages automatically in the plot window)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))

n_persons <- nrow(clus4data)

for (i in seq_len(n_persons)) {
  
  # extract responses + ID
  responses <- as.numeric(clus4data[i, resp_cols])
  person_id <- as.character(clus4data[i, "ID"])
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  # smooth
  smoothed_responses <- hanning_smooth(responses, n_iter = 9)
  
  # plot
  plot(ordered_items$xsi, smoothed_responses,
       type = "n",
       ylim = c(0, 1),
       xlim = c(0, 10),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = "black", lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = "black", pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = "red", pch = 19, cex = 1)
}

# reset (optional)
par(mfrow = c(1, 1))


# =========================================================
# Misfit Person Clustering (PAM, K = 3)
# Complete workflow with:
#   1. PAM fit
#   2. 2D MDS plot
#   3. Rotated 2D MDS plot
#   4. 3D MDS plot
#   5. Cluster-colored lines from medoid to points in 3D
#   6. Medoid lookup tables
#   7. Silhouette diagnostics
#   8. Medoid PRF plots (smoothed)
# =========================================================

# ----- Packages -----
library(plotly)

# ----- Fit PAM clustering -----
## ----- K = 3 -----
pam_fit3 <- pam(d_mis, k = 3)

person_dat_misfit$cluster_misfit3 <- pam_fit3$clustering
table(person_dat_misfit$cluster_misfit3)

# ----- MDS embedding -----
## ----- 2D MDS -----
xy <- cmdscale(d_mis, k = 2)

plot_df3 <- data.frame(
  ID      = as.character(person_dat_misfit$ID),
  Dim1    = xy[, 1],
  Dim2    = xy[, 2],
  Cluster = factor(pam_fit3$clustering)
)

## ----- Ensure IDs match distance labels -----
labs <- attr(d_mis, "Labels")
if (!is.null(labs) && length(labs) == nrow(plot_df3)) {
  plot_df3$ID <- labs
}

## ----- 3D MDS -----
xyz <- cmdscale(d_mis, k = 3)

plot_df3d <- data.frame(
  ID      = plot_df3$ID,
  Dim1    = xyz[, 1],
  Dim2    = xyz[, 2],
  Dim3    = xyz[, 3],
  Cluster = factor(pam_fit3$clustering)
)

# ----- Extract medoids -----
## ----- 2D medoids -----
med_df3 <- plot_df3[pam_fit3$id.med, , drop = FALSE]

## ----- 3D medoids -----
med_df3d <- plot_df3d[pam_fit3$id.med, , drop = FALSE]

## ----- Medoid IDs -----
medoid_ids3 <- plot_df3$ID[pam_fit3$id.med]
medoid_ids3

## ----- Medoid lookup table (2D) -----
medoid_lookup3 <- data.frame(
  Cluster    = factor(pam_fit3$clustering[pam_fit3$id.med]),
  Medoid_Row = pam_fit3$id.med,
  Medoid_ID  = plot_df3$ID[pam_fit3$id.med],
  Dim1       = plot_df3$Dim1[pam_fit3$id.med],
  Dim2       = plot_df3$Dim2[pam_fit3$id.med]
)

medoid_lookup3 <- medoid_lookup3[
  order(as.numeric(as.character(medoid_lookup3$Cluster))), 
]
print(medoid_lookup3)

## ----- Medoid lookup table (3D) -----
medoid_lookup3d <- data.frame(
  Cluster    = factor(pam_fit3$clustering[pam_fit3$id.med]),
  Medoid_Row = pam_fit3$id.med,
  Medoid_ID  = plot_df3d$ID[pam_fit3$id.med],
  Dim1       = plot_df3d$Dim1[pam_fit3$id.med],
  Dim2       = plot_df3d$Dim2[pam_fit3$id.med],
  Dim3       = plot_df3d$Dim3[pam_fit3$id.med]
)

medoid_lookup3d <- medoid_lookup3d[
  order(as.numeric(as.character(medoid_lookup3d$Cluster))), 
]
print(medoid_lookup3d)

# ----- Build segment data for 2D plot -----
## ----- Lookup table -----
med_lookup3 <- med_df3[, c("Cluster", "Dim1", "Dim2")]
names(med_lookup3) <- c("Cluster", "mx", "my")

## ----- Attach medoid coordinates -----
seg_df3 <- plot_df3
seg_df3$mx <- med_lookup3$mx[match(seg_df3$Cluster, med_lookup3$Cluster)]
seg_df3$my <- med_lookup3$my[match(seg_df3$Cluster, med_lookup3$Cluster)]

stopifnot(!anyNA(seg_df3$mx), !anyNA(seg_df3$my))

# =========================================================
# 2D Visualization
# =========================================================

# ----- Standard 2D MDS plot -----
p3 <- ggplot(plot_df3, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_segment(
    data = seg_df3,
    aes(x = mx, y = my, xend = Dim1, yend = Dim2, color = Cluster),
    inherit.aes = FALSE,
    alpha = 0.30,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_point(
    data = med_df3,
    aes(Dim1, Dim2),
    inherit.aes = FALSE,
    shape = 4,
    size = 6,
    stroke = 1.8,
    color = "black"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Misfit Person Clustering (PAM, K = 3)",
    x = "Dimension 1 (MDS)",
    y = "Dimension 2 (MDS)",
    color = "Cluster"
  )

print(p3)

# =========================================================
# Rotated 2D Plot
# =========================================================

# ----- Rotation function -----
rotate_2d <- function(x, y, angle_deg = 35) {
  theta <- angle_deg * pi / 180
  R <- matrix(
    c(cos(theta), -sin(theta),
      sin(theta),  cos(theta)),
    nrow = 2, byrow = TRUE
  )
  rot <- cbind(x, y) %*% R
  data.frame(x_rot = rot[, 1], y_rot = rot[, 2])
}

## ----- Apply rotation -----
angle_deg <- 35
rot <- rotate_2d(plot_df3$Dim1, plot_df3$Dim2, angle_deg)

plot_df3$Dim1_rot <- rot$x_rot
plot_df3$Dim2_rot <- rot$y_rot

## ----- Update medoid coordinates -----
med_df3 <- plot_df3[pam_fit3$id.med, ]

med_lookup3_rot <- med_df3[, c("Cluster", "Dim1_rot", "Dim2_rot")]
names(med_lookup3_rot) <- c("Cluster", "mx", "my")

seg_df3_rot <- plot_df3
seg_df3_rot$mx <- med_lookup3_rot$mx[match(seg_df3_rot$Cluster, med_lookup3_rot$Cluster)]
seg_df3_rot$my <- med_lookup3_rot$my[match(seg_df3_rot$Cluster, med_lookup3_rot$Cluster)]

stopifnot(!anyNA(seg_df3_rot$mx), !anyNA(seg_df3_rot$my))

# ----- Rotated plot -----
p3_rot <- ggplot(plot_df3, aes(Dim1_rot, Dim2_rot, color = Cluster)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_segment(
    data = seg_df3_rot,
    aes(x = mx, y = my, xend = Dim1_rot, yend = Dim2_rot, color = Cluster),
    inherit.aes = FALSE,
    alpha = 0.30,
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_point(
    data = med_df3,
    aes(Dim1_rot, Dim2_rot),
    inherit.aes = FALSE,
    shape = 4,
    size = 6,
    stroke = 1.8,
    color = "black"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = paste0("Misfit Person Clustering (PAM, K = 3, Rotated ", angle_deg, "�)"),
    x = "Rotated Dimension 1",
    y = "Rotated Dimension 2",
    color = "Cluster"
  )

print(p3_rot)

# =========================================================
# 3D Visualization
# =========================================================

# ----- Cluster colors (match ggplot palette) -----
cluster_colors <- c(
  "1" = "#F8766D",
  "2" = "#00BA38",
  "3" = "#619CFF"
)

# ----- Interactive 3D MDS -----
fig3d <- plot_ly(
  plot_df3d,
  x = ~Dim1,
  y = ~Dim2,
  z = ~Dim3,
  color = ~Cluster,
  colors = cluster_colors,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4)
)

# ----- Add cluster-colored lines from medoid to each point -----
## ----- Medoid lookup map -----
medoid_map3d <- medoid_lookup3d
rownames(medoid_map3d) <- as.character(medoid_map3d$Cluster)

for (i in seq_len(nrow(plot_df3d))) {
  
  cl <- as.character(plot_df3d$Cluster[i])
  med <- medoid_map3d[cl, ]
  
  fig3d <- fig3d %>%
    add_trace(
      x = c(med$Dim1, plot_df3d$Dim1[i]),
      y = c(med$Dim2, plot_df3d$Dim2[i]),
      z = c(med$Dim3, plot_df3d$Dim3[i]),
      type = "scatter3d",
      mode = "lines",
      line = list(
        width = 1,
        color = cluster_colors[cl]
      ),
      showlegend = FALSE,
      inherit = FALSE
    )
}

# ----- Add medoid markers and labels -----
for (i in seq_len(nrow(medoid_lookup3d))) {
  
  cl <- medoid_lookup3d$Cluster[i]
  
  fig3d <- fig3d %>%
    add_trace(
      x = medoid_lookup3d$Dim1[i],
      y = medoid_lookup3d$Dim2[i],
      z = medoid_lookup3d$Dim3[i],
      type = "scatter3d",
      mode = "markers+text",
      marker = list(
        size = 2,
        symbol = "x",
        color = "black"
      ),
      text = paste0("C", cl),
      textposition = "top center",
      textfont = list(color = "black"),
      showlegend = FALSE,
      inherit = FALSE
    )
}

fig3d <- fig3d %>%
  layout(
    title = "Misfit Person Clustering (PAM, K = 3, 3D MDS)",
    scene = list(
      xaxis = list(title = "Dimension 1 (MDS)"),
      yaxis = list(title = "Dimension 2 (MDS)"),
      zaxis = list(title = "Dimension 3 (MDS)")
    )
  )

fig3d

# =========================================================
# Silhouette Diagnostics
# =========================================================

# ----- Compute silhouette -----
sil3 <- silhouette(pam_fit3$clustering, d_mis)

## ----- Overall silhouette -----
mean(sil3[, 3])

## ----- Cluster means -----
tapply(sil3[, 3], pam_fit3$clustering, mean)

## ----- Silhouette plot -----
plot(sil3, border = NA, main = "Silhouette Plot (PAM, K = 3)")
abline(v = 0.25, col = "red", lty = 2)

## ----- Cluster summary -----
sil_by_cluster3 <- by(sil3[, 3], sil3[, 1], function(x) {
  c(
    n = length(x),
    mean_sil = mean(x),
    median_sil = median(x),
    min_sil = min(x),
    pct_neg = mean(x < 0),
    pct_gt_025 = mean(x > 0.25)
  )
})

sil_summary3 <- do.call(rbind, sil_by_cluster3)
sil_summary3

# =========================================================
# Medoid PRF Visualization
# =========================================================

# ----- Extract medoid respondents -----
medoid_data3 <- person_dat_misfit %>%
  filter(as.character(ID) %in% as.character(medoid_ids3)) %>%
  mutate(ID = as.character(ID))

# ----- Order items by difficulty -----
ordered_items <- mod1$item1[order(mod1$item1$xsi), c(1, 3)]
colnames(ordered_items) <- c("item", "xsi")

resp_cols <- names(person_dat_misfit)[2:19]

if (all(as.character(ordered_items$item) %in% names(person_dat_misfit))) {
  resp_cols <- as.character(ordered_items$item)
}

# ----- Hanning smoothing -----
## ----- Smoothing function -----
hanning_smooth <- function(x, n_iter = 9) {
  sm <- x
  n <- length(x)
  for (iter in 1:n_iter) {
    tmp <- sm
    for (j in 2:(n - 1)) {
      sm[j] <- (tmp[j - 1] + 2 * tmp[j] + tmp[j + 1]) / 4
    }
  }
  sm
}

# ----- Plot smoothed PRFs for medoids -----
par(mfrow = c(1, nrow(medoid_data3)))

for (i in seq_len(nrow(medoid_data3))) {
  
  responses <- as.numeric(medoid_data3[i, resp_cols])
  smoothed <- hanning_smooth(responses)
  
  title_i <- paste0("ID ", medoid_data3$ID[i], ": ", paste(responses, collapse = ""))
  
  plot(
    ordered_items$xsi, smoothed,
    type = "n",
    ylim = c(0, 1),
    xlab = "Item Difficulty",
    ylab = "Probability",
    main = title_i
  )
  
  lines(ordered_items$xsi, smoothed, lwd = 2)
  points(ordered_items$xsi, smoothed, pch = 22, cex = 1.4)
  points(ordered_items$xsi, responses, col = "red", pch = 19)
}

par(mfrow = c(1, 1))





# ----- Try 2 clusters (PAM) + visualization with IDs + silhouette -----

## ----- 1) Fit PAM with K = 2 -----
pam_fit2 <- pam(d_mis, k = 2)

# (optional) store cluster labels
person_dat_misfit$cluster_misfit2 <- pam_fit2$clustering
table(person_dat_misfit$cluster_misfit2)

## ----- 2) 2D embedding from distance matrix (MDS) -----
xy <- cmdscale(d_mis, k = 2)

plot_df2 <- data.frame(
  ID      = as.character(person_dat_misfit$ID),
  Dim1    = xy[, 1],
  Dim2    = xy[, 2],
  Cluster = factor(pam_fit2$clustering)
)

# (optional but recommended) ensure IDs align with distance labels
labs <- attr(d_mis, "Labels")
if (!is.null(labs) && length(labs) == nrow(plot_df2)) plot_df2$ID <- labs

## ----- 3) Medoids (actual observed individuals) -----
med_df2 <- plot_df2[pam_fit2$id.med, , drop = FALSE]
table(med_df2$Cluster)

## ----- 4) Build segment data (connect each point to its cluster medoid) -----
med_lookup2 <- med_df2[, c("Cluster", "Dim1", "Dim2")]
names(med_lookup2) <- c("Cluster", "mx", "my")

seg_df2 <- plot_df2
seg_df2$mx <- med_lookup2$mx[match(seg_df2$Cluster, med_lookup2$Cluster)]
seg_df2$my <- med_lookup2$my[match(seg_df2$Cluster, med_lookup2$Cluster)]
stopifnot(!anyNA(seg_df2$mx), !anyNA(seg_df2$my))

## ----- 5) Plot: points + connection lines + ellipses + medoids -----
p2 <- ggplot(plot_df2, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_segment(
    data = seg_df2,
    aes(x = mx, y = my, xend = Dim1, yend = Dim2, color = Cluster),
    inherit.aes = FALSE,
    alpha = 0.30,
    linewidth = 0.45,
    lineend = "round",
    show.legend = FALSE
  ) +
  # stat_ellipse(aes(group = Cluster, color = Cluster),
  #              linewidth = 1, level = 0.95) +
  geom_point(
    data = med_df2,
    aes(Dim1, Dim2),
    inherit.aes = FALSE,
    shape = 4,
    size = 6,
    stroke = 1.8,
    color = "black"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Misfit Person Clustering (PAM, K = 2)",
    x = "Dimension 1 (MDS)",
    y = "Dimension 2 (MDS)",
    color = "Cluster"
  )

print(p2)

## ----- 6) Cluster-specific silhouette performance for K = 2 -----
sil2 <- silhouette(pam_fit2$clustering, d_mis)
sil2_df <- as.data.frame(sil2)

sil_by_cluster2 <- lapply(split(sil2_df$sil_width, sil2_df$cluster), function(x) {
  data.frame(
    n = length(x),
    mean_sil = mean(x),
    median_sil = median(x),
    min_sil = min(x),
    pct_neg = mean(x < 0),
    pct_gt_025 = mean(x > 0.25)
  )
})

sil_summary2 <- do.call(rbind, sil_by_cluster2)
sil_summary2$Cluster <- rownames(sil_summary2)
rownames(sil_summary2) <- NULL

sil_summary2 <- sil_summary2[, c("Cluster", "n", "mean_sil", "median_sil", "min_sil", "pct_neg", "pct_gt_025")]
print(sil_summary2)


# add medoids (X only, no label)
fig2_3d <- fig2_3d %>%
  add_trace(
    data = med_df2d,
    x = ~Dim1,
    y = ~Dim2,
    z = ~Dim3,
    type = "scatter3d",
    mode = "markers",
    name = "Medoids",
    showlegend = TRUE,
    hoverinfo = "skip",
    marker = list(
      size = 7,
      color = "black",
      symbol = "x"
    )
  ) %>%
  layout(
    title = "Misfit Person Clustering (PAM, K = 2, 3D MDS)",
    scene = list(
      xaxis = list(title = "Dimension 1"),
      yaxis = list(title = "Dimension 2"),
      zaxis = list(title = "Dimension 3")
    )
  )
fig2_3d

## ------ K = 2 medoids, nonparametric smoothed PRF ------
# look at 2-cluster medoids
pam_fit2$medoids
medoid_points2 <- plot_df2[pam_fit2$medoids, ]
medoid_points2

medoid_ids2 <- plot_df2$ID[pam_fit2$id.med]
medoid_ids2

# Filter rows from person_dat_misfit where ID matches
library(dplyr)

medoid_data2 <- person_dat_misfit %>%
  filter(as.character(ID) %in% as.character(medoid_ids2)) %>%
  mutate(ID = as.character(ID))

# Item difficulty from low to high
ordered_items <- mod1$item1[order(mod1$item1$xsi), c(1, 3)]

# response columns
resp_cols <- names(person_dat_misfit)[2:19]

if (all(as.character(ordered_items$item) %in% names(person_dat_misfit))) {
  resp_cols <- as.character(ordered_items$item)
}

n_persons <- nrow(medoid_data2)

par(mfrow = c(1, 2))

for (i in 1:n_persons) {
  
  # Extract responses
  responses <- as.numeric(medoid_data2[i, resp_cols])
  
  # Person ID
  person_id <- medoid_data2$ID[i]
  
  # Title showing ID + raw response string
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  # Initialize smoothed responses
  smoothed_responses <- responses
  n_items <- length(responses)
  
  # Hanning smoothing (9 iterations)
  for (iteration in 1:9) {
    temp_responses <- smoothed_responses
    for (j in 2:(n_items - 1)) {
      smoothed_responses[j] <- (temp_responses[j - 1] +
                                  2 * temp_responses[j] +
                                  temp_responses[j + 1]) / 4
    }
  }
  
  # Plot
  plot(ordered_items$xsi, smoothed_responses,
       type = "n",
       ylim = c(0, 1),
       xlim = c(min(ordered_items$xsi), max(ordered_items$xsi)),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = "black", lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = "black", pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = "red", pch = 19, cex = 1)
}

par(mfrow = c(1, 1))


## K = 2 medoids, parametric PRF
medoid_ids2 <- plot_df2$ID[pam_fit2$id.med]
medoid_ids2

# Pull the medoid rows from person-level data
library(dplyr)

medoid_data2 <- person_dat_misfit %>%
  filter(as.character(ID) %in% as.character(medoid_ids2)) %>%
  mutate(ID = as.character(ID))

# Item difficulties in difficulty order
ordered_items <- mod1$item1[order(mod1$item1$xsi), c(1, 3)]
colnames(ordered_items) <- c("item", "xsi")
item_difficulty <- ordered_items$xsi

# Make sure response columns are aligned with ordered item names
resp_cols <- names(person_dat_misfit)[2:19]

if (all(as.character(ordered_items$item) %in% names(person_dat_misfit))) {
  resp_cols <- as.character(ordered_items$item)
}

# Plot PRFs for each medoid
par(mfrow = c(1, 2))

for (i in seq_len(nrow(medoid_data2))) {
  
  responses <- as.numeric(medoid_data2[i, resp_cols])
  person_id <- as.character(medoid_data2$ID[i])
  theta_i   <- as.numeric(medoid_data2$ability[i])  # must exist
  
  # Rasch parametric PRF
  prf_i <- plogis(theta_i - item_difficulty)
  
  # Title = ID + raw response string
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  plot(item_difficulty, prf_i,
       type = "n",
       ylim = c(0, 1),
       xlim = c(min(item_difficulty), max(item_difficulty)),
       xlab = "Item Difficulty",
       ylab = "Probability",
       main = response_title)
  
  # expected PRF
  lines(item_difficulty, prf_i, col = "black", lwd = 2)
  points(item_difficulty, prf_i, col = "black", pch = 22, cex = 1.4)
  
  # observed 0/1 responses
  points(item_difficulty, responses, col = "red", pch = 19, cex = 0.9)
}

par(mfrow = c(1, 1))

## ----- Cluster-specific silhouette summary for K = 2, 3, 4, 5 -----

k_vals <- 2:5

sil_summary_all <- lapply(k_vals, function(k) {
  
  # fit PAM
  pam_k <- pam(d_mis, k = k)
  
  # silhouette values
  sil_k <- silhouette(pam_k$clustering, d_mis)
  sil_k_df <- as.data.frame(sil_k)
  
  # summarize by cluster
  sil_by_cluster_k <- lapply(split(sil_k_df$sil_width, sil_k_df$cluster), function(x) {
    data.frame(
      n = length(x),
      mean_sil = mean(x),
      median_sil = median(x),
      min_sil = min(x),
      pct_neg = mean(x < 0),
      pct_gt_025 = mean(x > 0.25)
    )
  })
  
  sil_summary_k <- do.call(rbind, sil_by_cluster_k)
  sil_summary_k$Cluster <- rownames(sil_summary_k)
  rownames(sil_summary_k) <- NULL
  
  sil_summary_k$K <- k
  
  sil_summary_k <- sil_summary_k[, c("K", "Cluster", "n", "mean_sil", "median_sil", "min_sil", "pct_neg", "pct_gt_025")]
  sil_summary_k
})

sil_summary_all <- do.call(rbind, sil_summary_all)

# optional: clean row names
rownames(sil_summary_all) <- NULL

print(sil_summary_all)


















##################################
# plot person response function ##
##################################
ordered_items <- mod1$item1[order(mod1$item1$xsi), c(1, 3)]# dif from low to high

################
# person 8692  # 
################
## non-parametric
responses8692 <- as.numeric(person_dat_misfit[21, 2:19])  # Ensure responses are numeric
response_title8692 <- paste(responses8692, collapse = "")   # Create the response title

smoothed_responses <- responses8692 # Hanning smoothing function 
n_items <- length(responses8692)
for (iteration in 1:9) {
  temp_responses <- smoothed_responses
  for (j in 2:(n_items-1)) {
    smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
  }
}

# plot the smoothed responses
plot(ordered_items$xsi, smoothed_responses, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty',       # x-axis: dif
     ylab = 'Probability',  # y-axis: prob
     main = response_title8692)

lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2) # add line
points(ordered_items$xsi, smoothed_responses, col = 'black',  pch = 22, cex = 1.5) # add smooth response
points(ordered_items$xsi, responses8692, col = 'red', pch = 19, cex = 1) # observed response
#-----------------------------------------------------------------------------------------
## parametric PRF 
# Extract observed responses
responses8692 <- as.numeric(person_dat_misfit[21, 2:19])
response_title8692 <- paste(responses8692, collapse = "")   

# Get person ability and item difficulty
theta21 <- person_dat_misfit$ability[21]
item_difficulty <- ordered_items$xsi

# Compute PRF: predicted probabilities from Rasch model
prf21 <- exp(theta21 - item_difficulty) / (1 + exp(theta21 - item_difficulty))

# Plot regular PRF (not smoothed)
plot(item_difficulty, prf21, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty', 
     ylab = 'Probability', 
     main = response_title8692)

# Add PRF line and square markers
lines(item_difficulty, prf21, col = 'black', lwd = 2)
points(item_difficulty, prf21, col = 'black', pch = 22, cex = 1.5)

# Add actual response points
points(item_difficulty, responses8692, col = 'red', pch = 19, cex = 1)

###############
# person 3014 # 
###############
# par(mfrow = c(1, 1), mar = c(2, 2, 2, 1))  
responses3014 <- as.numeric(person_dat_misfit[10, 2:19])  # Ensure responses are numeric
response_title3014 <- paste(responses3014, collapse = "")   # Create the response title

smoothed_responses <- responses3014 
n_items <- length(responses3014)

for (iteration in 1:9) {
  temp_responses <- smoothed_responses
  
  for (j in 2:(n_items-1)) {
    smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
  }
}

# plot the smoothed responses
plot(ordered_items$xsi, smoothed_responses, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty',       # x-axis: dif
     ylab = 'Probability',  # y-axis: prob
     main = response_title3014)

lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2) # add line
points(ordered_items$xsi, smoothed_responses, col = 'black',  pch = 22, cex = 1.5) # add smooth response
points(ordered_items$xsi, responses3014, col = 'red', pch = 19, cex = 1) # observed response
#-------------------------------------------------------
## parametric PRF
# Extract observed responses
responses3014 <- as.numeric(person_dat_misfit[10, 2:19])  # Ensure responses are numeric
response_title3014 <- paste(responses3014, collapse = "")   # Create the response title 

# Get person ability and item difficulty
theta10 <- person_dat_misfit$ability[10]
item_difficulty <- ordered_items$xsi

# Compute PRF: predicted probabilities from Rasch model
prf10 <- exp(theta10 - item_difficulty) / (1 + exp(theta10 - item_difficulty))

# Plot regular PRF (not smoothed)
plot(item_difficulty, prf10, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty', 
     ylab = 'Probability', 
     main = response_title3014)

# Add PRF line and square markers
lines(item_difficulty, prf10, col = 'black', lwd = 2)
points(item_difficulty, prf10, col = 'black', pch = 22, cex = 1.5)

# Add actual response points
points(item_difficulty, responses3014, col = 'red', pch = 19, cex = 1)

################
#  person 128  #
################
## non-parametric 1.28,  (statistically fit, but graphically no)
responses128 <- as.numeric(person_dat[9, 2:19])  # Ensure responses are numeric
response_title128 <- paste(responses128, collapse = "")   # Create the response title

smoothed_responses <- responses128
n_items <- length(responses128)
for (iteration in 1:8) {
  temp_responses <- smoothed_responses
  
  for (j in 2:(n_items-1)) {
    smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
  }
}

# plot the smoothed responses
plot(ordered_items$xsi, smoothed_responses, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty',       # x-axis: dif
     ylab = 'Probability',  # y-axis: prob
     main = response_title128)

lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2) # add line
points(ordered_items$xsi, smoothed_responses, col = 'black',  pch = 22, cex = 1.5) # add smooth response
points(ordered_items$xsi, responses128, col = 'red', pch = 19, cex = 1) # observed response
#-----------------------------------------------------------------------------------------
## parametric PRF
# Extract observed responses
responses128 <- as.numeric(person_dat[9, 2:19])  # Ensure responses are numeric
response_title128 <- paste(responses128, collapse = "")   # Create the response title

# Get person ability and item difficulty
theta128 <- person_dat$ability[9]
item_difficulty <- ordered_items$xsi

# Compute PRF: predicted probabilities from Rasch model
prf128 <- exp(theta128 - item_difficulty) / (1 + exp(theta128 - item_difficulty))

# Plot regular PRF (not smoothed)
plot(item_difficulty, prf128, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty', 
     ylab = 'Probability', 
     main = response_title128)

# Add PRF line and square markers
lines(item_difficulty, prf128, col = 'black', lwd = 2)
points(item_difficulty, prf128, col = 'black', pch = 22, cex = 1.5)

# Add actual response points
points(item_difficulty, responses128, col = 'red', pch = 19, cex = 1)

###############
#  person 1256  #
###############
## non-parametric 0.49, 0.26
responses1256 <- as.numeric(person_dat[73, 2:19])  # Ensure responses are numeric
response_title1256 <- paste(responses1256, collapse = "")   # Create the response title

smoothed_responses <- responses1256
n_items <- length(responses1256)
for (iteration in 1:8) {
  temp_responses <- smoothed_responses
  
  for (j in 2:(n_items-1)) {
    smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
  }
}

# plot the smoothed responses
plot(ordered_items$xsi, smoothed_responses, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty',       # x-axis: dif
     ylab = 'Probability',  # y-axis: prob
     main = response_title1256)

lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2) # add line
points(ordered_items$xsi, smoothed_responses, col = 'black',  pch = 22, cex = 1.5) # add smooth response
points(ordered_items$xsi, responses1256, col = 'red', pch = 19, cex = 1) # observed response
#-----------------------------------------------------------------------------------------
## parametric PRF
# Extract observed responses
responses1256 <- as.numeric(person_dat[73, 2:19])  # Ensure responses are numeric
response_title1256 <- paste(responses1256, collapse = "")   # Create the response title

# Get person ability and item difficulty
theta1256 <- person_dat$ability[73]
item_difficulty <- ordered_items$xsi

# Compute PRF: predicted probabilities from Rasch model
prf1256 <- exp(theta1256 - item_difficulty) / (1 + exp(theta1256 - item_difficulty))

# Plot regular PRF (not smoothed)
plot(item_difficulty, prf1256, 
     type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
     xlab = 'Item Difficulty', 
     ylab = 'Probability', 
     main = response_title1256)

# Add PRF line and square markers
lines(item_difficulty, prf1256, col = 'black', lwd = 2)
points(item_difficulty, prf1256, col = 'black', pch = 22, cex = 1.5)

# Add actual response points
points(item_difficulty, responses1256, col = 'red', pch = 19, cex = 1)

#############
#  loop it  #
#############
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))  

n_persons <- nrow(person_dat_misfit)
for (i in 1:n_persons) {
  
  # Extract responses
  responses <- as.numeric(person_dat_misfit[i, 2:19])  # Ensure responses are numeric
  person_id <- person_dat_misfit[i, "ID"]  # Get the person's ID
  response_title <- paste(person_id, ": ", paste(responses, collapse = ""))
  
  smoothed_responses <- responses  
  n_items <- length(responses)
  
  # iterate 9 times
  for (iteration in 1:9) {
    temp_responses <- smoothed_responses
    for (j in 2:(n_items-1)) {
      smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
    }
  }
  
  # plot it
  plot(ordered_items$xsi, smoothed_responses, 
       type = 'n', 
       ylim = c(0, 1), 
       xlim = c(0, 10),  
       xlab = 'Item Difficulty',       # x-axis: dif
       ylab = 'Probability',  # y-axis: prob
       main = response_title)
  
  lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = 'black', pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = 'red', pch = 19, cex = 1)
}
##################################################
### add in item-measure correlation
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))
correlations <- numeric(nrow(person_dat_misfit))

item_difficulty <- ordered_items$xsi
n_persons <- nrow(person_dat_misfit)

for (i in 1:n_persons) {
  
  responses <- as.numeric(person_dat_misfit[i, 2:19])
  person_id <- person_dat_misfit[i, "ID"]
  response_title <- paste0(person_id, ": ", paste(responses, collapse = ""))
  
  # smoothing function
  smoothed_responses <- responses  
  n_items <- length(responses)
  
  # 9 iteration
  for (iteration in 1:9) {
    temp_responses <- smoothed_responses
    for (j in 2:(n_items-1)) {
      smoothed_responses[j] <- (temp_responses[j-1] + 2 * temp_responses[j] + temp_responses[j+1]) / 4
    }
  }
  
  # item-measure correlation
  correlations[i] <- cor(responses, item_difficulty, use = "complete.obs")
  
  # plots
  response_title <- paste0("ID: ", person_id, 
                           " | Corr: ", round(correlations[i], 2), 
                           "\nResp: ", paste(responses, collapse = ""))
  
  # plots
  plot(item_difficulty, smoothed_responses, 
       type = 'n', ylim = c(0, 1), xlim = c(min(item_difficulty), max(item_difficulty)),  
       xlab = 'Item locations', 
       ylab = 'Pr(x=1)', 
       main = response_title, cex.main = 0.8)
  
  lines(item_difficulty, smoothed_responses, col = 'black', lwd = 2)
  points(item_difficulty, smoothed_responses, col = 'black', pch = 22, cex = 1.5)
  points(item_difficulty, responses, col = 'red', pch = 19, cex = 1) }


person_dat_misfit$correlation_with_xsi <- round(correlations, digits = 2)
head(person_dat_misfit)
person_dat_misfit$correlation_with_xsi

###################################################
## PRF
par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))

# prepare statistic I need
correlations <- numeric(nrow(person_dat_misfit))
item_difficulty <- ordered_items$xsi
n_persons <- nrow(person_dat_misfit)
n_items <- length(item_difficulty)
person_theta <- person_dat_misfit$ability  # Rasch ability 

for (i in 1:n_persons) {
  
  # response, person ID
  responses <- as.numeric(person_dat_misfit[i, 2:(1 + n_items)])
  person_id <- person_dat_misfit[i, "ID"]
  
  # PRF calculation, Rasch model
  theta <- person_theta[i]
  prf <- exp(theta - item_difficulty) / (1 + exp(theta - item_difficulty))
  
  # item-measure correlation
  correlations[i] <- cor(responses, item_difficulty, use = "complete.obs")
  
  # plot title
  response_title <- paste0("ID: ", person_id, 
                           " | Corr: ", round(correlations[i], 2), 
                           "\nResp: ", paste(responses, collapse = ""))
  
  # ??? ONLY ONE plot() call, with correct syntax
  plot(item_difficulty, prf, 
       type = 'n', 
       ylim = c(0, 1), 
       xlim = c(0, 10),   # ??? ????????????????????????
       xlab = 'Item Difficulty',       
       ylab = 'Probability',  
       main = response_title)
  
  # black line and square 
  lines(item_difficulty, prf, col = 'black', lwd = 1.5)
  points(item_difficulty, prf, col = 'black', pch = 22, bg = "white", cex = 1.3)
  
  # observed response: red dots(0/1)
  points(item_difficulty, responses, col = 'red', pch = 19, cex = 1.2)
}
#-----------------------------------------------------------
# response-measure correlation 
person_dat_misfit$correlation_with_xsi <- round(correlations, digits = 2)

head(person_dat_misfit)
person_dat_misfit$correlation_with_xsi
#-----------------------------------------------------------























####################################################
## use 50% correct data: person_dat_misfit_50
n_persons <- nrow(person_dat_misfit_50)
for (i in 1:n_persons) {
  
  # Extract responses
  responses <- as.numeric(person_dat_misfit_50[i, 2:19])  # Ensure responses are numeric
  person_id <- person_dat_misfit_50[i, "ID"]  # Get the person's ID
  response_title <- paste(person_id, ": ", paste(responses, collapse = ""))
  
  smoothed_responses <- responses  
  n_items <- length(responses)
  
  # iterate 20 times
  for (iteration in 1:9) {
    temp_responses <- smoothed_responses
    for (j in 2:(n_items-1)) {
      smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
    }
  }
  
  # plot it
  plot(ordered_items$xsi, smoothed_responses, 
       type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
       xlab = 'Item locations', 
       ylab = 'Pr(x=1)', 
       main = response_title)  # Use the current person's response title
  
  lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = 'black', pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = 'red', pch = 19, cex = 1)
}






############################
#  group by the sum score  #  
############################
### same sum score, look at their response function
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 17 
#  2  5 12 16  7 11  4  6  8  9  7  5  3  5  2  1

person_dat_misfit$sumscore= rowSums(person_dat_misfit[,2:19])

par(mfrow = c(3, 2), mar = c(2, 2, 2, 1))
sumscore_3 <- subset(person_dat_misfit, sumscore == 3) #dataset for sumscore = 3
n_persons <- nrow(sumscore_3)

# Loop through each person's responses in the subset
for (i in 1:n_persons) {
  
  # Extract responses and ID for person i
  responses <- as.numeric(sumscore_3[i, 2:19])  # Ensure responses are numeric
  person_id <- sumscore_3[i, "ID"]  # Get the person's ID
  
  # Create the response title with the person's ID and response pattern
  response_title <- paste(person_id, ": ", paste(responses, collapse = ""))
  
  # Initialize smoothed responses with original values
  smoothed_responses <- responses  
  n_items <- length(responses)
  
  # Apply the smoothing method 20 times
  for (iteration in 1:9) {
    temp_responses <- smoothed_responses
    for (j in 2:(n_items-1)) {
      smoothed_responses[j] <- (temp_responses[j-1] + 2*temp_responses[j] + temp_responses[j+1]) / 4
    }
  }
  
  # Plot the smoothed responses for each person in the subset
  plot(ordered_items$xsi, smoothed_responses, 
       type = 'n', ylim = c(0, 1), xlim = c(0, 10),  
       xlab = 'Item locations', 
       ylab = 'Pr(x=1)', 
       main = response_title)  # Use the current person's response title
  
  lines(ordered_items$xsi, smoothed_responses, col = 'black', lwd = 2)
  points(ordered_items$xsi, smoothed_responses, col = 'black', pch = 22, cex = 1.5)
  points(ordered_items$xsi, responses, col = 'red', pch = 19, cex = 1)
}

#########################################################################
# item responses correlated with item difficulty (misreporting persons) #
#########################################################################
correlations <- numeric(nrow(person_dat_misfit))

# Loop through each person's responses in the dataset
for (i in 1:nrow(person_dat_misfit)) {
  
  responses <- as.numeric(person_dat_misfit[i, 2:19])
  correlations[i] <- cor(responses, mod1$xsi, use = "complete.obs")
}

person_dat_misfit$correlation_with_xsi <- round(correlations, digits = 2)
View(person_dat_misfit)
person_dat_misfit$correlation_with_xsi

# Question: for all persons with positive correlation:
# sort by the correlation and do the plot by order: positive to negative
# personfit plot: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/PerFit/PerFit.pdf
# by sum score or by correlation








































###################################################
## Conditional Maximum Likelihood (CML) Approach ##
###################################################
library(eRm)
rasch_model_cml= RM(d1_dicho,sum0 = TRUE)
summary(rasch_model_cml)

item_params <- coef(rasch_model_cml) # item parameter
plotjointICC(rasch_model_cml, item.subset =  1:18, cex = .6)

plotPImap(rasch_model_cml, cex.gen = .55)



## plot it
par(mfrow = c(3, 3), mar = c(2, 2, 2, 1))  # Adjust margins to make everything fit neatly

for(i in 1:103) {
  ability <- person_dat_misfit$ability[i]
  probabilities <- rasch_prf(ability, ordered_items$xsi)
  
  responses <- as.numeric(person_dat_misfit[i, 2:19])
  response_title <- paste(paste(responses, collapse=""))  # Create a single string of responses
  
  plot(ordered_items$xsi, probabilities, type = 'l', lwd = 2, col = 'black', ylim = c(0, 1),
       xlab = 'Item Difficulty (Location)', ylab = 'Probability of Correct Response', 
       main = response_title)  # Use response string as plot title
  
  # Add points for the theoretical probabilities
  points(ordered_items$xsi, probabilities, col = 'black', pch = 17)
  
  # Add points for the observed responses
  points(ordered_items$xsi, responses, col = 'red', pch = 19)
}

########################################################################


# probability function 
rasch_prf <- function(theta, b) {
  prob <- exp(theta - b) / (1 + exp(theta - b))
  return(prob)
}
