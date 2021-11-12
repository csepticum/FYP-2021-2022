WEB5986 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH3gcov/WEB5986.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

LOY73 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH3gcov/LOY7-23-36.gcov', 
                    header = FALSE,  
                    sep = "\t", 
                    as.is = TRUE, 
                    comment.char = "#", 
                    quote = "", 
                    col.names = c("rand",	"x0",	"y0"))


beg <- 1153243
star <- 1173243
sto <- 1319461
end <- 1339461

makeUHdf <- function(df){
  avg <- mean(df$y0)
  UHout <- df[beg:end,]
  UHout$y0 <- (UHout$y0)/avg
  return(UHout)
}

utiplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "tomato3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), star, sto, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("UTI89 reads aligned against MG1655 genome (UH3 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

mgplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "skyblue3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), star, sto, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("MG1655 reads aligned against MG1655 genome (UH3 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

plotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "grey38") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), star, sto, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("WEB5986 reads aligned against MG1655 genome (UH3 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

# out6 <- ggplot() +
#   geom_line(data = uh6a, aes(x = x0, y = y0, color = "WEB5958"), alpha = 0.8) +
#   geom_line(data = uh6b, aes(x = x0, y = y0, color = "WEB6026"), alpha = 0.8) +
#   geom_line(data = uh6c, aes(x = x0, y = y0, color = "WEB6028"), alpha = 0.8) +
#   scale_x_continuous(name = "MG1655 Coordinates", breaks = c(beg, star, sto, end)) +
#   ylab("Sequencing Depth") +
#   ggtitle("UH6 Clones Sequencing Depth") +
#   theme_light() +
#   theme(
#     axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
#     legend.title = element_blank()
#   )
# show(out6)

uti3 <- utiplotter(makeUHdf(uti))
mg3 <- mgplotter(makeUHdf(mg))
uh3a <- plotter(makeUHdf(LOY73))
uh3b <- plotter(makeUHdf(WEB5986))

out3 <- ggarrange(uti3, mg3, uh3a, uh3b,
                  nrow = 4, ncol =1)
show(out3)