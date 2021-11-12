mg <- read.table('/Users/sylvester/Desktop/GenomeVizTest/MG1655.gcov', 
                 header = FALSE,  
                 sep = "\t", 
                 as.is = TRUE, 
                 comment.char = "#", 
                 quote = "", 
                 col.names = c("rand",	"x0",	"y0"))


uti <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UTI89.gcov', 
                  header = FALSE,  
                  sep = "\t", 
                  as.is = TRUE, 
                  comment.char = "#", 
                  quote = "", 
                  col.names = c("rand",	"x0",	"y0"))

WEB5965 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5965.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

LOY7 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/LOY7-23-4.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5979 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5979.gcov', 
                   header = FALSE,  
                   sep = "\t", 
                   as.is = TRUE, 
                   comment.char = "#", 
                   quote = "", 
                   col.names = c("rand",	"x0",	"y0"))

WEB5982 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5982.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5980 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5980.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5981 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5981.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5977 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5977.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5978 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH2gcov/WEB5978.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))


#given the gcov alignment, takes out the coordinates i want and return as a dataframe ready for ggplotting
makeUHdf <- function(df, star, sto){
  avg <- mean(df$y0)
  star <- star - 20000
  sto <- sto + 20000
  UHout <- df[star:sto,]
  UHout$y0 <- (UHout$y0)/avg
  return(UHout)
}

utiplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "tomato3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 566333, 582018, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("UTI89 reads aligned against MG1655 genome (UH2 region)") +
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
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 566333, 582018, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("MG1655 reads aligned against MG1655 genome (UH2 region)") +
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
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(df$x0), 566333, 582018, max(df$x0))) +
    ylab("Normalized Sequencing Depth") +
    ggtitle("LOY7-23-4 reads aligned against MG1655 genome (UH2 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

uti1 <- makeUHdf(uti, 566333,582018)
mg1 <- makeUHdf(mg, 566333, 582018)
uh2a <- makeUHdf(WEB5979, 566333, 582018)
uh2b <- makeUHdf(WEB5982, 566333, 582018)
uh2c <- makeUHdf(WEB5980, 566333, 582018)
uh2d <- makeUHdf(WEB5981, 566333, 582018)
uh2e <- makeUHdf(WEB5977, 566333, 582018)
uh2f <- makeUHdf(WEB5978, 566333, 582018)
uh2g <- makeUHdf(LOY7, 566333, 582018)

out2 <- ggplot() +
  #geom_line(data = uti1, aes(x = x0, y = y0, color = "UTI89"), alpha = 0.8) +
  #geom_line(data = mg1, aes(x = x0, y = y0, color = "MG1655"), alpha = 0.8) +
  geom_line(data = uh2a, aes(x = x0, y = y0, color = "WEB5979"), alpha = 0.8) +
  geom_line(data = uh2b, aes(x = x0, y = y0, color = "WEB5982"), alpha = 0.8) +
  geom_line(data = uh2c, aes(x = x0, y = y0, color = "WEB5980"), alpha = 0.8) +
  geom_line(data = uh2d, aes(x = x0, y = y0, color = "WEB5981"), alpha = 0.8) +
  geom_line(data = uh2e, aes(x = x0, y = y0, color = "WEB5977"), alpha = 0.8) +
  geom_line(data = uh2f, aes(x = x0, y = y0, color = "WEB5978"), alpha = 0.8) +
  geom_line(data = uh2g, aes(x = x0, y = y0, color = "LOY7-23-4"), alpha = 0.8) +
  scale_x_continuous(name = "MG1655 Coordinates", breaks = c(min(uh2a$x0),566333, 582018, max(uh2a$x0))) +
  ylab("Normalized Sequencing Depth") +
  ggtitle("UH2 Clones Sequencing Depth") +
  theme_light() +
  theme(
    axis.text.x = element_text(face = c("plain", "bold", "bold", "plain"), colour = c("black", "red", "red", "black")),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
show(out2)

uti2 <- utiplotter(makeUHdf(uti, 566333, 582018))
mg2 <- mgplotter(makeUHdf(mg, 566333,582018))
loy7 <- plotter(makeUHdf(LOY7, 566333, 582018))

out <- ggarrange(uti2, mg2, loy7,
                 nrow = 3, ncol =1)
show(out)
