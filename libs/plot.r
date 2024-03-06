#!/usr/bin/env Rscript
library(ggplot2)

args <- commandArgs(T)

# df_reads
df_stats <- read.table(file = paste0(args[1], ".stats"), col.names = c("Key", "Value"))
df_reads <- df_stats[c(2, 3, 5, 4), ]
df_reads$Key <- c("All reads", "Valid reads", "gRNA reads", "Unknow reads")
df_reads$Key <- factor(df_reads$Key, levels = c("All reads", "Valid reads", "gRNA reads", "Unknow reads"))
names(df_reads)[names(df_reads) == "Value"] <- "Counts"

# df_gRNAs(only gRNAs, no unknow)
df_all <- read.table(file = paste0(args[1], ".percentage"), sep = "\t", header = T)
df_gRNAs <- df_all[df_all$gene_id != "unknow", ]
df_gRNAs_detected <- df_all[(df_all$gene_id != "unknow") & (df_all$counts != 0), ]
df_detected <- df_all[df_all$counts != 0, ]
t <- sum(df_gRNAs$counts)
df_gRNAs$gRNAs_percentage <- (df_gRNAs$counts / t) * 100
t_detected <- sum(df_gRNAs_detected$counts)
df_gRNAs_detected$gRNAs_percentage <- (df_gRNAs_detected$counts / t_detected) * 100

# plot1 reads
reads <- ggplot(df_reads, aes(Key, Counts)) +
  geom_col(width = 0.7, color = "white", fill = "lightblue") +
  geom_text(aes(label = Counts), vjust = -0.2) +
  xlab("") +
  scale_y_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12)
  )
ggsave("reads.pdf", reads)

# plot2(1/2) frequency
frequency <- ggplot(df_gRNAs, aes(x = reorder(gene_id, gRNAs_percentage), y = gRNAs_percentage)) +
  theme_bw() +
  ylab("Frequency(%)") +
  scale_y_continuous(labels = scales::comma) +
  geom_point(size = 0.4) +
  geom_hline(aes(yintercept = mean(gRNAs_percentage)), linetype = 2, color = "red", show.legend = TRUE) +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage, 0.75))) +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage, 0.25))) +
  geom_hline(aes(yintercept = median(gRNAs_percentage))) +
  scale_fill_discrete(limits = c("Average(red), Quartile(black)")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom"
  )
frequency <- frequency + aes(x = seq(1, length(gRNAs_percentage))) + scale_x_continuous() + xlab("Each kind of gRNAs")

# add text
x_range <- layer_scales(frequency)$x$range$range
y_range <- layer_scales(frequency)$y$range$range
Mean <- mean(df_gRNAs$gRNAs_percentage)
Var <- var(df_gRNAs$gRNAs_percentage)
CV <- sqrt(Var) / Mean
frequency <- frequency + geom_text(
  x = x_range[2] * 0.8,
  y = y_range[2] * 0.92,
  aes(label = paste0("Var = ", Var, "\nCV = ", CV))
)
ggsave("frequency.png", frequency)

# plot2(2/2) frequency_detected
frequency_detected <- ggplot(df_gRNAs_detected, aes(x = reorder(gene_id, gRNAs_percentage), y = gRNAs_percentage)) +
  theme_bw() +
  ylab("Frequency(%)") +
  scale_y_continuous(labels = scales::comma) +
  geom_point(size = 0.4) +
  geom_hline(aes(yintercept = mean(gRNAs_percentage)), linetype = 2, color = "red", show.legend = TRUE) +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage, 0.75))) +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage, 0.25))) +
  geom_hline(aes(yintercept = median(gRNAs_percentage))) +
  scale_fill_discrete(limits = c("Average(red), Quartile(black)")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom"
  )
frequency_detected <- frequency_detected + aes(x = seq(1, length(gRNAs_percentage))) + scale_x_continuous() + xlab("Each kind of gRNAs")

# add text
x_range <- layer_scales(frequency_detected)$x$range$range
y_range <- layer_scales(frequency_detected)$y$range$range
Mean <- mean(df_gRNAs_detected$gRNAs_percentage)
Var <- var(df_gRNAs_detected$gRNAs_percentage)
CV <- sqrt(Var) / Mean
frequency_detected <- frequency_detected + geom_text(
  x = x_range[2] * 0.8,
  y = y_range[2] * 0.92,
  aes(label = paste0("Var = ", Var, "\nCV = ", CV))
)
ggsave("frequency_detected.png", frequency_detected)


# plot3(1/2) frequency_distribution
frequency_distribution <- ggplot(df_gRNAs, aes(x = gRNAs_percentage)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = mean(gRNAs_percentage)), linetype = 2) +
  xlab("Frequency(%)") +
  ylab("Counts") +
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave("frequency_distribution.pdf", frequency_distribution)

# plot3(2/2) frequency_distribution_detected
frequency_distribution_detected <- ggplot(df_gRNAs_detected, aes(x = gRNAs_percentage)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = mean(gRNAs_percentage)), linetype = 2) +
  xlab("Frequency(%)") +
  ylab("Counts") +
  scale_x_continuous(labels = scales::comma) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
ggsave("frequency_distribution_detected.pdf", frequency_distribution_detected)

# plot4 accumulative_unknow_percentage
accumulative_unknow_percentage <- ggplot(df_detected, aes(x = seq(1, length(percentage)), y = accumulative_unknow_percentage)) +
  geom_point(size = 0.2) +
  xlab("Accumulative kinds of sequences") +
  ylab("Unknow percentage") +
  theme_bw()
ggsave("accumulative_unknow_percentage.pdf", accumulative_unknow_percentage)
