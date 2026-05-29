library(spatstat.geom)     
library(spatstat.explore)  
library(spatstat.random)   
library(spatstat.model)

data(gorillas, package = "spatstat.data")
X1 = gorillas
X1_major = unmark(X1[X1$marks$group == "major"])
X1_minor = unmark(X1[X1$marks$group == "minor"])


library(ggplot2)
win <- as.polygonal(gorillas$window)
win_df <- do.call(rbind, lapply(seq_along(win$bdry), function(i) {
  data.frame(
    x = win$bdry[[i]]$x,
    y = win$bdry[[i]]$y,
    id = i
  )
}))

major_df <- data.frame(x = X1_major$x, y = X1_major$y, group = "Major")
minor_df <- data.frame(x = X1_minor$x, y = X1_minor$y, group = "Minor")
points_df <- rbind(major_df, minor_df)


plot.res = 600
png("gorilla_groups.png", width = 6.2 * plot.res, height = 4.5 * plot.res, res = 600)
ggplot() +
  geom_polygon(data = win_df, aes(x = x, y = y, group = id),
               fill = NA, color = "black") +
  geom_point(data = points_df, aes(x = x, y = y, color = group), size = 1, alpha = 0.5) +
  labs(
    title = "Gorillas nests with groups",
    x = "X coordinate",
    y = "Y coordinate",
    color = "Group"
  ) +
  scale_color_manual(values = c("Major" = "black", "Minor" = "red")) +
  coord_equal() +
  theme_minimal()
dev.off()




