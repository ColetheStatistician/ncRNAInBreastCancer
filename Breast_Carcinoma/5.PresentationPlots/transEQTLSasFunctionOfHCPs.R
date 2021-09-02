library(ggplot2)
results <- readRDS("iteratePlotData.rds")

# Get num transEQTLs
num_transEQTLs <- rep(as.numeric(), length(results))
for (i in 1:length(results)) {
  num_transEQTLs[i] <- results[[i]]$trans$neqtls
}

# make df to plot
df <- data.frame("transEQTLs" = num_transEQTLs,
                 "num_HCPs" = 1:length(num_transEQTLs))

pdf("transEQTLsByHCP.pdf")
ggplot(data = df, aes(x = num_HCPs, y = transEQTLs)) +
  geom_line(lwd=2, color = "lightblue") +
  geom_point(size=2, color="gold") +
  xlab("Number of HCPs") +
  ylab("Significant trans-EQTLs") +
  theme(aspect.ratio = 1) +
  theme_bw() + 
  scale_y_continuous(limits = c(58000, 65000),
                     breaks = seq(58000, 65000, 1000),
                     minor_breaks = seq(58000, 65000, 500)) +
  scale_x_continuous(limits = c(1, 8),
                     breaks = 1:8)
dev.off()
