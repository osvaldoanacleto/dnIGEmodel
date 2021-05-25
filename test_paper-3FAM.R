#####################################################################
############################## 3FAM #################################
#####################################################################
#-------------------#
# NOTACAO DO ARTIGO #
#-------------------#

# T.indiv = Total number of individuals  (N)
# N = Number of Families---------------- (n.family)
# n = Family Size -----------------------(size.family)
# n.w = Group size ----------------------(group_size)
# n.g = numbers of groups ---------------(ngroups)
# n_g.f = Numbers of groups per family---(ng.family)
# n.b = number of blocks ----------------(n.block)
# n_f.b = number of family per block ----(nf.block)
# n_g.b = number of groups per blocks ---(ng.block)
#----------------#
# VALORES '3FAM' #
#----------------#
N = 5760
n.family = 72
size.family = 80 # size family
group_size = 60
ngroups = 96
ng.family = 4
n.block = 8 #completos
nf.block = 9
ng.block = 12

72*80
#-------------------#
# Gerando os 'Pais' #
#-------------------#

pai <- rep(NA, N)

k<-1
for (j in 1:n.family) {
  pai[seq(k,k+(size.family-1),1)] <- rep(j, size.family)
  k <- k+size.family
}


pai 

#-------------------------------------------#
# Gerando as 'Maes' para cada um dos 'Pais' #
#-------------------------------------------#

mae <- rep(NA, N)

i <-1
for(p in 1:n.family){
  for (m in 1:size.family) {
    mae[i] <- paste(paste("mae", m, sep = "_"), p, sep = "_pai_")
    i <- i+1
  }
}

mae

pais <- as.data.frame(cbind(pai, mae))

head(pais)
dim(pais)
#----------------------------------#
# Gerando um filho para cada 'Mae' # 
# Meio-irmÃ£os                      #
#----------------------------------#

filhos <- rep(NA, N)

nf <- 1
for (f in 1:(N)) {
  filhos[f] <- paste(paste("filho", nf, sep = "_"), mae[f], sep = "_")
  
  if (nf == size.family){
    nf <- 1
  }else{
    nf <- nf+1
  }
}

filhos


familia <- as.data.frame(cbind(pais, filhos))
head(familia)

#--------------------------------------------------------#
# Inserinfo a VariÃ¡vel 'groups' no dataset               #
# Esta variÃ¡vel serÃ¡ determinada pelo algoritmo a seguir #
#--------------------------------------------------------#

familia$groups <- rep(NA, N)

#----------------------------------#
# FORMAR BLOCKS                    #
# Gerando os Grupos AleatÃ³riamente #
# apÃ³s definidos os blocos         #
#----------------------------------#

for (block in 1:ceiling(n.block)) {
  
  # Identificando o range dos pais
  p1 <- nf.block*(block-1) +1 
  p2 <- nf.block*block
  
  fam <- as.numeric(familia[familia$pai %in% seq(p1, p2), "pai"])
  max_fam <- max(fam)
  min_fam <- min(fam)
  
  # Criando combinaÃ§Ãµes de familias
  comb_pais <- combn(seq(min_fam, max_fam, 1), 3)
  
  comb_pais <- comb_pais[,c(1, 14, 23, 28, 36, 43, 47, 53, 56, 61, 73, 76)]
  
  # SeleÃ§Ã£o aleatoria para cada grupo
  for (jj in 1:ncol(comb_pais)) {
    combinacao <- comb_pais[, jj]
    
    indiv1 <- sample(rownames(familia[familia$pai == combinacao[1] & is.na(familia$groups), ]), group_size/3, replace = F)
    indiv2 <- sample(rownames(familia[familia$pai == combinacao[2] & is.na(familia$groups), ]), group_size/3, replace = F)
    indiv3 <- sample(rownames(familia[familia$pai == combinacao[3] & is.na(familia$groups), ]), group_size/3, replace = F)
    
    familia[c(indiv1, indiv2, indiv3), "groups"] <- 
      paste(paste(combinacao[1], combinacao[2], sep = "|"), combinacao[3], sep = "|")
    
  }
}


unique(familia$groups)
