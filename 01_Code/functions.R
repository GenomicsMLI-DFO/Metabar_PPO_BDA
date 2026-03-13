
pcoa.fig <- function(dist.dna, MOTUS.info, TAXO){
  pcoa.dist <- wcmdscale(dist.dna, eig = TRUE)
  vect.pcoa <- scores(pcoa.dist, choices = c(1:2))
  df.pcoa <- data.frame(QueryAccVer = rownames(vect.pcoa), vect.pcoa)
  df.motus <- left_join(df.pcoa, MOTUS.info)
  df.motus$TAXO <- as.factor(df.motus[,TAXO])
  pcoa.align.obj1 <- ggplot(df.motus) +
    geom_point(aes(x = Dim1, y = Dim2, col = TAXO, size = count_postmetbaR),
               alpha = 0.7) +
    labs(color = TAXO) +
    theme_classic()
  return(pcoa.align.obj1)
}

pcoa.fig2 <- function(dist.dna, MOTUS.info, TAXO){
  vect.pcoa <- dist.dna[,1:2]
  df.pcoa <- data.frame(QueryAccVer = rownames(vect.pcoa), vect.pcoa)
  df.motus <- left_join(df.pcoa, MOTUS.info)
  df.motus$TAXO <- as.factor(df.motus[,TAXO])
  pcoa.align.obj1 <- ggplot(df.motus) +
    geom_point(aes(x = PC_1, y = PC_2, col = TAXO, size = count_postmetbaR),
               alpha = 0.7) +
    labs(color = TAXO) +
    theme_classic()
  return(pcoa.align.obj1)
}
