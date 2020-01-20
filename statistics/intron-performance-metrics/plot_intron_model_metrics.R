prec <- c(0.24793, 0.27636, 0.32897, 0.38282)
acc <- c(0.841, 0.8687, 0.9065, 0.9286)
recall <- c(0.9145, 0.8407, 0.65547, 0.4644)
deg <- c(5,6,7,8)
metrics <- data.frame(degree=deg, recall=recall, precision=prec, accuracy=acc)

metrics.df <- melt(metrics, id.vars="deg")
ggplot(metrics.df, aes(deg,value, col=variable)) +  geom_point() + geom_line()