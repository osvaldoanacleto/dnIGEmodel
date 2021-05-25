##########################################
####    ALGORITMO PARA FORMAR GROUPS   ###
##########################################

##---------------------------------------------------------------##
## DATASET EXEMPLO do PAPER-  MEIO-IRMAO ##
##---------------------------------------------------------------##

# Semente
set.seed(1234)

#-------------------#
# NOTACAO DO ARTIGO #
#-------------------#

# T.indiv = Total number of individuals  (N)
# N = Number of Families---------------- (n.family)
# n = Family Size -----------------------(size.family)
# n.w = Group size ----------------------(group_size)
# n.g = numbers of groups ---------------(ngroups)
# n_g.f = Numbers of groups per family---(ng.family)
# n.b = number of blocks ----------------(n.blocks)
# n_f.b = number of family per block ----(nf.block)
# n_g.b = number of groups per blocks ---(ng.block)


#----------------#
# VALORES '2FAM' #
#----------------#
N = 5760
n.family = 40
size.family = 150 # size family
group_size = 60
ngroups = 96
ng.family= 5
n.blocks = 6 # 6 completos e o 7 incompleto
nf.block = 6
ng.block = 15
#--------- Maiores Famílias
maior.family=floor(n.blocks)*nf.block 
size.maior.family=size.family
#---------menores famílias
menor.family=(n.family)-(maior.family) 
size.menor.family=((ngroups*group_size)-(n.family-menor.family)*size.family)/(menor.family)
5760/150
#-------------------#
# Gerando os 'Pais' #
#-------------------#

pai <- rep(NA, N)
#n_pai <- N/n.family

k<-1
for (j in 1:(n.family-menor.family)) {
  pai[seq(k,k+(size.family-1),1)] <- rep(j, size.family)
  k <- k+size.family
}
for (j in (n.family-3):n.family) {
  pai[seq(k,k+(size.menor.family-1),1)] <- rep(j, size.menor.family)
  k <- k+size.menor.family
}

length(pai) 

#-------------------------------------------#
# Gerando as 'Maes' para cada um dos 'Pais' #
#-------------------------------------------#

mae <- rep(NA, N)

i <-1
for(p in 1:n.family){
  if (p <= (n.family-menor.family)){
    for (m in 1:size.family) {
      mae[i] <- paste(paste("mae", m, sep = "_"), p, sep = "_pai_")
      i <- i+1
    }
  }else{
    for (m in 1:size.menor.family) {
      mae[i] <- paste(paste("mae", m, sep = "_"), p, sep = "_pai_")
      i <- i+1
    }
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
  
  if (pais[f, 1] <=  n.family- menor.family){
    if (nf == size.family){
      nf <- 1
    }else{
      nf <- nf+1
    }
  }else{
    if (nf == size.menor.family){
      nf <- 1
    }else{
      nf <- nf+1
    }
  }
}

filhos


familia <- as.data.frame(cbind(pais, filhos))
head(familia)
dim(familia)

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

n.blocks= N/(ng.block*group_size)
n.blocks

for (block in 1:ceiling(n.blocks)) {
  
  # Identificando o range dos pais
  p1 <- nf.block*(block-1) +1 
  p2 <- nf.block*block
  
  fam <- as.numeric(familia[familia$pai %in% seq(p1, p2), "pai"])
  max_p <- max(fam)
  min_p <- min(fam)
  
  # Criando combinaÃ§Ãµes de familias
  comb_pais <- combn(seq(min_p, max_p, 1), 2)
  
  # SeleÃ§Ã£o aleatoria para cada grupo
  for (jj in 1:ncol(comb_pais)) {
    combinacao <- comb_pais[, jj]
    
    indiv1 <- sample(rownames(familia[familia$pai == combinacao[1] & is.na(familia$groups), ]), 0.5*group_size, replace = F)
    indiv2 <- sample(rownames(familia[familia$pai == combinacao[2] & is.na(familia$groups), ]), 0.5*group_size, replace = F)
    
    
    familia[c(indiv1, indiv2), "groups"] <- paste(combinacao[1], combinacao[2], sep = "|")
    
  }
}

# TESTE
for (gpr in unique(familia$groups)) {
  print(table(familia[familia$groups == gpr,]$pai))
}


table(familia[c((N-359):N),]$pai)
table(familia[c((N-359):(N-360)),]$pai)


familia[c((N-180):N),]

unique(familia$groups)

