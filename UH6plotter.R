WEB5958 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH6gcov/WEB5958.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6026 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH6gcov/WEB6026.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6028 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH6gcov/WEB6028.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))
beg <- 4480715
star <- 4500715
sto <- 4561491
end <- 4581491

makeUHdf <- function(df){
  avg <- mean(df$y0)
  UHout <- df[beg:end,]
  UHout$y0 <- (UHout$y0)/avg
  return(UHout)
}

utiplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "tomato3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), star, 4536278, sto, max(df$x0))) +
    geom_vline(xintercept = 4536278, color = "forestgreen") +
    ylab("Normalized Sequencing Depth") +
    ggtitle("UTI89 reads aligned against MG1655 genome (UH6 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

mgplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "skyblue3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), star, 4536278, sto, max(df$x0))) +
    geom_vline(xintercept = 4536278, color = "forestgreen") +
    ylab("Normalized Sequencing Depth") +
    ggtitle("MG1655 reads aligned against MG1655 genome (UH6 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

plotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "grey38") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), star, 4536278, sto, max(df$x0))) +
    geom_vline(xintercept = 4536278, color = "forestgreen") +
    ylab("Normalized Sequencing Depth") +
    ggtitle("WEB6026 reads aligned against MG1655 genome (UH6 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "red", "black")),
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

uti6 <- utiplotter(makeUHdf(uti))
mg6 <- mgplotter(makeUHdf(mg))
uh6a <- plotter(makeUHdf(WEB6026))

out6 <- ggarrange(uti6, mg6, uh6a,
                  nrow = 3, ncol =1)
show(out6)