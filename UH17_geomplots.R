mg <- read.table('/Users/sylvester/Desktop/GenomeVizTest/MG.gcov', 
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

WEB5960 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB5960.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5962 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB5962.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5965 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB5965.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB5969 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB5969.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6028 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB6028.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6030 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB6030.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6031 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB6031.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

WEB6036 <- read.table('/Users/sylvester/Desktop/GenomeVizTest/UH1gcov/WEB6036.gcov', 
                      header = FALSE,  
                      sep = "\t", 
                      as.is = TRUE, 
                      comment.char = "#", 
                      quote = "", 
                      col.names = c("rand",	"x0",	"y0"))

#given the gcov alignment, takes out the coordinates i want and return as a dataframe ready for ggplotting
makeUHdf <- function(df){
  avg <- mean(df$y0)
  UH7 <- df[4640441:4641652,]
  UH1 <- df[1:1392,]
  UH7$x0 <- c(1:nrow(UH7))
  UH1$x0 <- c(1213:2604)
  out <- rbind(UH7, UH1)
  out$y0 <- (out$y0)/avg
  return(out)
}


utiplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "tomato3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(1, 1001, 1213, 1514, 1604, 2604), labels = c("4640441", "4641441", "1", "302", "392", "1392")) +
    geom_vline(xintercept = 1212, color = "forestgreen") +
    ylab("Normalized Sequencing Depth") +
    ggtitle("UTI89 reads aligned against MG1655 genome (UH1 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "forestgreen", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

mgplotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "skyblue3") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(1, 1001, 1213, 1514, 1604, 2604), labels = c("4640441", "4641441", "1", "302", "392", "1392")) +
    geom_vline(xintercept = 1212, color = "forestgreen") +
    ylab("Normalized Sequencing Depth") +
    ggtitle("MG1655 reads aligned against MG1655 genome (UH1 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "forestgreen", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

plotter <- function(df){
  out <- ggplot() +
    geom_line(data = df, aes(x = x0, y = y0), color = "grey38") +
    scale_x_continuous(name = "MG1655 Coordinates", breaks = c(1, 1001, 1213, 1514, 1604, 2604), labels = c("4640441", "4641441", "1", "302", "392", "1392")) +
    geom_vline(xintercept = 1212, color = "forestgreen") +
    ylab("Normalized Sequencing Depth") +
    ggtitle("WEB5960 reads aligned against MG1655 genome (UH1 region)") +
    theme_light() +
    theme(
      axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "forestgreen", "red", "black")),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  show(out)
  return(out)
}

uti1 <- makeUHdf(WEB5962)
uti2 <- makeUHdf(WEB5965)
uti3 <- makeUHdf(WEB5969)
uti4 <- makeUHdf(WEB6028)
uti5 <- makeUHdf(WEB6030)
uti6 <- makeUHdf(WEB6031)
uti7 <- makeUHdf(WEB6036)

out <- ggplot() +
  geom_line(data = uti1, aes(x = x0, y = y0, color = "WEB5962"), alpha = 0.8) +
  geom_line(data = uti2, aes(x = x0, y = y0, color = "WEB5965"), alpha = 0.8) +
  geom_line(data = uti3, aes(x = x0, y = y0, color = "WEB5969"), alpha = 0.8) +
  geom_line(data = uti4, aes(x = x0, y = y0, color = "WEB6028"), alpha = 0.8) +
  geom_line(data = uti5, aes(x = x0, y = y0, color = "WEB6030"), alpha = 0.8) +
  geom_line(data = uti6, aes(x = x0, y = y0, color = "WEB6031"), alpha = 0.8) +
  geom_line(data = uti7, aes(x = x0, y = y0, color = "WEB6036"), alpha = 0.8) +
  scale_x_continuous(name = "MG1655 Coordinates", breaks = c(1, 1001, 1213, 1514, 1604, 2604), labels = c("4640441", "4641441", "1", "302", "392", "1392")) +
  ylab("Normalized Sequencing Depth") +
  ggtitle("UH1 Clones Sequencing Depth") +
  theme_light() +
  theme(
    axis.text.x = element_text(face = c("plain", "bold", "bold", "bold", "bold", "plain"), colour = c("black", "red", "forestgreen", "forestgreen", "red", "black")),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
show(out)

uti1 <- makeUHdf(uti)
uti1 <- utiplotter(uti1)
mg1 <- makeUHdf(mg)
mg1 <- mgplotter(mg1)
uh1a <- makeUHdf(WEB5960)
uh1a <- plotter(uh1a)

out <- ggarrange(uti1, mg1, uh1a,
                 nrow = 3, ncol = 1)
show(out)
