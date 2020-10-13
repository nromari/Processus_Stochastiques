#Simualtions en Biologie - Méthodes stochastiques M2BI Université Paris Diderot
#-----------------------------AHMAR-ERRAS Noura--------------------------------
#---------------------------BARRY Fatoumata Binta-------------------------------

require(MASS)

#Exercice 2 : Methodes Monte Carlo
#*********************************

#Il s'agit ici de comparer les methodes d'integration Monte Carlo.
#Les méthodes de Monte Carlo sont des modèles de simulation stochastiques.
#Trois methodes seront evaluees : MC "tirage noir ou blanc", MC simple, et MC suivant l'importance.

#1) La fonction à simuler g(x)
#-----------------------------

#Afin de mettre en place efficacement les modeles de simulation il faut bien connaître la fonction.
#En effet pour simuler une variable aleataoire continue il est necessaire de calculer son integrale.

g = function(x) {                 #g est une fonction exponentielle.
  return((exp(x)-1)/(exp(1)-1))
}

#Visualisation de g sur un interval [0-2]
val_x = seq(from=0, to=2, len=1000)
val_y = g(val_x)
plot(val_x, val_y, type="l", main="g(x) sur l'interval [0-2]")      #L'allure de la courbe est bien exponentielle.

#Verification de la valeur de l'integrale de g sur [0-2]
I_True = (exp(2)-3)/(exp(1)-1)      #On integre g de 0 à 2
print(I_True)
#La valeur réelle de l'intégrale de g sur [0-2] est 2.554328 ce qui correspond à l'air sous g de 0 à 2.
abline(v=2, lty=2, col="red")

integrate(g, lower=0, upper=2)
#la valeur de l'integrale de g pour l'interval [0-2] est de 2.554328
#g(2) = integ(g)[0-2] --> On retrouve bien la même valeur.

#2) Methode de Monte Carlo par tirage noir ou blanc sur l'interval [0-2]
#-----------------------------------------------------------------------

#Le principe est basé sur le tirage au hasard des points.
#Il faut donc définir la surface et échantillonner l'espace de manière uniforme
#Si l'intégrale correspond à un rectangle les X et Y seront tirés de manière uniforme
#La délimitation du rectangle passe par le choix d'un majorant

#Choix du majorant : g(2) est le majorant optimale car sur l'interval choisit g(2) est le plus petit majorant possible.

MC_BN = function(n) {
  m = g(2)    #choix du majorant
  u = runif(n, min = 0, max = 2)    #Tirage uniforme dans g sur l'intervale [0-2]
  v = runif(n, min = 0, max = m)    #Tirage uniforme dans le rectangle sur l'intervale [0-m]
  s = m*2       #La loi uniforme est entre [0-1]. On applique un facteur 2 au majorant car on évalue sur l'intervale [0-2]
  ns = sum(v<g(u))    #nombre de succes = nombre de tirage en dessous de g(u)
  return(ns/n*s)
}

MC_BN(1000) #A lancer plusieurs fois
# Nous lançons 3 fois : 2.535868 , 2.751529, 2.587924
#Les valeurs obtenues sont proches de la valeur attendue (pour rappel valeur vraie de l'intégrale = 2.554328)

#3) Methode de Monte Carlo par echantillonage simple sur l'interval [0-2]
#-----------------------------------------------------------------------

#La méthode de Monte Carlo dite simple est une methode d'echantillonnage par décomposition de la fonction g
#g est alors exprimé comme le produit d'une fonction quelquonque et sa fonction de densité
#Dans ce cas, le calcul d'intégrale s'approche du calcul de l'Espérance.
#Dans le cas simple, la variable à échantillonner suit une loi uniforme

MC_S = function(n) {
  x = runif(n, min=0, max=2)
  return(2*sum(g(x))/n)
}
#le 2 correspond au dunif de la fonction sur l'intervale [0-2]. Ici on peut le sortir, et le multiplier par l'espérance de g(x)

MC_S(1000) #A lancer plusieurs fois
#Nous lançons 3 fois : 2.684689 , 2.559881, 2.60517
#Les valeurs obtenues sont proches de la valeur attendue (pour rappel valeur vraie de l'intégrale = 2.554328)

#4) Methode de Monte Carlo par echantillonnage suivant l'importance sur l'interval [0-2]
#-------------------------------------------------------------------------------------

#Cette méthode diffère de la méthode simple par la méthode d'échantillonage.
#La variable à échantillonner ne suit pas une loie uniforme.
#Dans ce cas il faut echantillonner pour se rapprocher au mieux de g.
#Cette méthode permet d'echantillonner dans les zones qui ont le plus d'importance.

#Pour construire l'echantillonnage, il faut créer une fonction qui permet de tirer et une fonction qui permet d'avoir la densité.

#La fonction pour le tirage(r_X) doit être du type bêta et va échantilloner entre 0 et 1 (pas [0-2] car impossible avec lois de type bêta).
#La solution pour rester dans l'interval [0-2]  c'est de tout multiplier par 2 or de la fonction.
r_X = function(n) {
  return(2*rbeta(n, shape1 = 2, shape2 = 1))# On créé le Y=aX avec a=2
}

#fonction de densité de probabilité de x
d_X = function(x) {
  return(1/2*dbeta(x/2, shape1 = 2, shape2 = 1))# Il faut diviser par 2 à cause de la transformation Y=aX car dbeta ne peut pas accepter des valeurs hor de [0-1]
}

#Validation graphique de r_X et d_X
truehist(r_X(100000), main="Monte Carlo suivant l'importance \n Loi beta de paramètres (2,1)", xlab="variables aléatoires", ylab="Probabilités", col="yellow")
val_x = seq(from=0, to=2, len=1000)
val_y = d_X(val_x)
lines(val_x, val_y, col="red")

#Integration
MC_I_1 = function(n) {
  x = r_X(n)   #tirage dans la loi beta de la variable aléatoire
  xx = g(x)/d_X(x)  #xx est ici la fonction quelconque. Pour l'exprimer on divise g par sa densité de probabilité
  return(mean(xx))  #I=E[xx]
}
MC_I_1(1000) #A lancer plusieurs fois
#Nous lançons 3 fois : 2.59096 , 2.558245, 2.538063
#Les valeurs obtenues sont proches de la valeur attendue (pour rappel valeur vraie de l'intégrale = 2.554328)

#5) Methode de Monte Carlo par echantillonnage suivant l'importance sur l'interval [0-2] - Lois beta V2
#------------------------------------------------------------------------------------------------------

#Modifions les paramètres de la loi beta pour proposer une autre fonction de densité
r_X_V2 = function(n) {
  return(2*rbeta(n, shape1 = 2.7, shape2 = 1.5))
}
d_X_V2 = function(x) {
  return(1/2*dbeta(x/2, shape1 = 2.7, shape2 = 1.5))
}

#Validation graphique de r_X_V2 et d_X_V2
truehist(r_X_V2(100000), main="Monte Carlo suivant l'importance \n Lois beta de paramètres (2.7,1.5)", xlab="variables aléatoires", ylab="Probabilités", col="yellow")
val_x = seq(from=0, to=2, len=1000)
val_y = d_X_V2(val_x)
lines(val_x, val_y, col="red")

#Integration
MC_I_2 = function(n) {
  x = r_X_V2(n) 
  xx = g(x)/d_X_V2(x) 
  return(mean(xx))
}
MC_I_2(1000) #A lancer plusieurs fois
#Nous lançons 3 fois : 2.586437 , 2.588031, 2.561564
#Les valeurs obtenues sont proches de la valeur attendue (pour rappel valeur vraie de l'intégrale = 2.554328)

#6) Comparer l'efficacité des différentes méthodes
#--------------------------------------------------

#Les méthodes seront évaluées sur la même taille d'échantillon n(petit), et sur le même nombre de répétition nrep (grand)
#Evaluation graphique et par comparaison des mse (mean square error) et des moyennes
par(mfrow=c(2,2))
n=100
n_rep=1000000
val_x = seq(from=0, to=4, len=1000)

#Monte Carlo par tirage blanc ou noir (MC_BN)
res_MC_BN = rep(NA, n_rep)
for (i in 1:n_rep) {
  res_MC_BN[i] = MC_BN(n)
}

truehist(xlim=range(I_True), res_MC_BN, col = "yellow", main="Intégration de g estimée par Monte Carlo \n tirage blanc/noir",xlab="variables aléatoires", ylab="Probabilités")
abline(v=I_True, lty=2, col="red")

mse_MC_BN = round(mean((res_MC_BN-I_True)^2),2)
moy_MC_BN = round(mean(res_MC_BN),2)
print(moy_MC_BN)
print(mse_MC_BN)
#La methode de Monte Carlo par tirage blanc/noir appliquée à g a un mse de 0.12 et une moeyenne de 2.55

#Monte Carlo par échantillonnage simple (MC_S)
res_MC_S = rep(NA, n_rep)
for (i in 1:n_rep) {
  res_MC_S[i] = MC_S(n)
}

truehist(xlim=range(I_True), res_MC_S, col = "yellow", main="Intégration de g estimée par Monte Carlo \n echantillonnage simple", xlab="variables aléatoires", ylab="Probabilités")
abline(v=I_True, lty=2, col="red")

mse_MC_S = round(mean((res_MC_S-I_True)^2),2)
moy_MC_S = round(mean(res_MC_S),2)
print(moy_MC_S)
print(mse_MC_S)
#La methode de Monte Carlo par tirage blanc/noir appliquée à g a un mse de 0.04 et une moyenne de 2.55

#Monte Carlo par échantillonnage d'importance V1 (MC_I_1)
res_MC_I_1 = rep(NA, n_rep)
for (i in 1:n_rep) {
  res_MC_I_1[i] = MC_I_1(n)
}

truehist(xlim=range(I_True), res_MC_I_1, col = "yellow", main="Intégration de g estimée par Monte Carlo \n echantillonnage d'importance V1", xlab="variables aléatoires", ylab="Probabilités")
abline(v=I_True, lty=2, col="red")

mse_MC_I_1 = round(mean((res_MC_I_1-I_True)^2),2)
moy_MC_I_1 = round(mean(res_MC_I_1),2)
print(moy_MC_I_1)
print(mse_MC_I_1)
#La methode de Monte Carlo par tirage blanc/noir appliquée à g a un mse de 0 et une moyenne de 2.55

#Monte Carlo par échantillonnage d'importance V2 (MC_I_2)
res_MC_I_2 = rep(NA, n_rep)
for (i in 1:n_rep) {
  res_MC_I_2[i] = MC_I_2(n)
}

truehist(xlim=range(I_True), res_MC_I_2, col = "yellow", main="Intégration de g estimée par \n Monte Carlo echantillonnage d'importance V2", xlab="variables aléatoires", ylab="Probabilités")
abline(v=I_True, lty=2, col="red")

mse_MC_I_2 = round(mean((res_MC_I_2-I_True)^2),2)
moy_MC_I_2 = round(mean(res_MC_I_2),2)
print(moy_MC_I_2)
print(mse_MC_I_2)
#La methode de Monte Carlo par tirage blanc/noir appliquée à g a un mse de 0.03 et une moyenne de 2.55

par(mfrow=c(1,1))

#Exercice 3 : Simulation de variables aléatoires
#*********************************

# 1) Loi de poisson: Implémentation de la fonction my_rpois de simulation d'une variable
# de poisson en utilisant une méthode d'inversion de la fonction de repartition
# -------------------------------------------------------------------------------------

# Implémentation de my_rpois_one(lambda) : fonction qui renvoie une seule valeur tirée aléatoirement dans une loi de Poisson de paramètre Lambda
# Nous prendrons 2.4 comme valeur de lambde

lambda = 2.4

my_rpois_one = function(lambda){ # avec lambda >0
  u = runif(1) # une seule valeur issue de la loi uniforme
  x = 0
  while (ppois(x, lambda) < u){ # tant que la probabilité d'avoir x, avec une loi de poisson
    # de parametre lambda est inférieure à la valeur u tirée de la loi uniforme, on incrémente x de 1,
    # sinon, on renvoie x
    x = x + 1 
  }
  return(x)
}

#Exemple

cat("valeur suivant une loi de poisson de parametre",lambda,":", my_rpois_one(lambda))

# Impléméntation de la fonction my_rpois(n,lambda) , qui fait appelle n fois la fonction my_rpois_one et renvoie un vecteur de n valeurs tirées aléatoirement et indépendamment dans une loi de Poisson de paramètre lambda.

my_rpois = function(n, lambda) {
  n_pois = rep(NA, len=n)# initialision d'un vecteur de n valeurs
  for (i in 1:n) {
    n_pois[i] = my_rpois_one(lambda) # affectation progression des éléments du vecteur X, suivant une loi de poisson
  }
  return(n_pois)
}

#Exemple
n = 100000  # nombre de repetitions
n_pois = my_rpois(n, lambda) # vecteurs de n valeurs aléatoires tirées de la loi de posson my_r
head(n_pois) # pour visualiser que nous avons bien des entiers, et qui sont aléatoires
plot(table(n_pois)/n, col="blue", main="Distribution de variables aléatoires \n suivant une loi de poisson", xlab="variables aléatoires x", ylab="Probabilités")


# Verification de notre fonction de simulation my_rpois en utilisant la fonction rpois de R

val_x = seq(from=0, to=12) # entiers de 0 à 8
val_y = dpois(val_x, lambda) # valeurs suivant une loi de poisson
points(val_x, val_y, col="red") # representation des points sur notre courbe de my_rpois
# Notre simulation est bien en adéquation avec la loi de poisson de paramètre lambda (2.4). Les deux courbes (rouge et bleu) ont la meme allure

# 2) Loi discrète: Implémentation d'une fonction my_rdiscret de simulatio n de la variable aléatoire réelle discrète X définie par une loi de probabilité X
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

#Comme pour la fonction my_rpois_one, la fonction my_rdiscret_one() renvoie une seule valeur tirée aéatoitement dans la loi de X

my_rdiscret_one = function(){
  u = runif(1)
  x = c(-3,1.3,7,15.2) # valeur de x dans l'enoncé
  dx = c(0.1, 0.4, 0.3, 0.2) # probabilité d'avoir x
  px = cumsum(dx) # somme cumulée des probabilités dx
  i = 1
  while(px[i] < u){
    i=i+1
  }
  return(x[i])
}

cat("valeur suivant la loi discrète X:", my_rdiscret_one()) # une des valeurs X de l'enoncé

#La fonction my_rdiscret renvoie n  valeurs tirées de la loi X

my_rdiscret = function(n){
  n_discret = rep(NA, n) # initialisation d'un vecteur à n valeurs
  for(i in 1:n){
    n_discret[i] = my_rdiscret_one()
  }
  return(n_discret)
}

## Verification de notre implementation en faisant un tirage aléatoire avec remise des valeurs x, avec leurs probabilités associées
n = 10000
plot(table(my_rdiscret(n))/n, col="blue", main="Distribution de variables aléatoires X \n suivant une loi discrète", xlab="variables X", ylab="Probabilités")
valX = c(-3, 1.3, 7, 15.2) # valeurs à tirer
points(valX, table(sample(valX, n, replace=TRUE, prob=c(0.1, 0.4, 0.3, 0.2)))/n, col='red')
#Les deux courbes sont correlées


# 3) Loi de Laplace: Méthode de transformation
# Implémentation d'une fonction my_rlaplace de simulation d'une variable de Laplace en utilisant une méthode d'inversion de la fonction de répartition.
# --------------------------------------------

### Fonction de densité my_dlaplace()
#Si x < 0, my_dlaplace(x) = 1/2 exp(x)
#Si x > 0, my_dlaplace(x) = 1/2 exp(-x)

my_dlaplace = function(x) {
  return(1/2*exp(-abs(x)))
}
# Exemple
val_x = seq(from=-4, to=4, len=100) # 100 valeur de -4 à 4
val_y = my_dlaplace(val_x) 
plot(val_x, val_y, col="red", type="l", main="Fonction de densité de la loi de Laplace", xlab="valeurs de X", ylab="densité de laplace de X")

# Fonction de répartition my_plaplace()
# Si x < 0, my_plaplace(x) = 1/2 exp(x)
# Si x > 0, my_plaplace(x) = 1/2(2-exp(-x))
### Implémentation de my_plaplace_one() pour une seule valeur

my_plaplace_one = function(x) {
  if (x < 0) {
    p = 1/2*exp(x)
  } else {
    p = 1/2*(2-exp(-x))
  }
  return(p)
}


#### Impléméntation de my_plaplace() pour n valeurs

my_plaplace = function(x_vec) {
  n = length(x_vec)
  p_vec = rep(NA, n)
  for (i in 1:n) {
    p_vec[i] = my_plaplace_one(x_vec[i])
  }
  return(p_vec)
}

my_plaplace(c(-0.8, -0.6, -0.4, -0.2, 0, 2, 4, 6, 8, 10))

# Exemple
val_x = seq(from=-4, to=4, len=100)
val_y = my_plaplace(val_x)
plot(val_x, val_y, col="red", type="l", main="Fonction de repartion de la loi de Laplace", xlab="variables X", ylab="Probabilités")

### Fonction quantile qlaplace()
#Si p < 1/2, qlaplace(p) = log(2*p)
#Si p > 1/2, qlaplace(p) = -log(1-2*p)

my_qlaplace_one = function(p) {
  if (p < 1/2) {
    x = log(2*p)
  } else {
    x = -log(2-2*p)
  }
  return(x)
}


my_qlaplace = function(p_vec) {
  n = length(p_vec)
  x_vec = rep(NA, n)
  for (i in 1:n) {
    x_vec[i] = my_qlaplace_one(p_vec[i])
  }
  return(x_vec)
}

# Exemple

val_x = seq(from = 0, to = 1, len = 100)
val_y = my_qlaplace(val_x)
plot(val_x, val_y, col = "red", type = "l", main = "Fonction quantile de Laplace")

# La fonction quantile est bien la reciproque de la fonction de repartition

### Fonction rlaplace(): génération de nombres aléatoires suivant la loi de Laplace

my_rlaplace = function(n) {
  u = runif(n)
  x = my_qlaplace(u)
  return(x)
}

# Exemple
n = 10000
x = my_rlaplace(n)
head(x)

# Verification de notre implémentation de Laplace
truehist(x, col="yellow", main="Distribution des variables suivant une loi de Laplace", xlab="variables aléatoires", ylab="probabilités" )

val_x = seq(from = -7, to = 7, len = 1000)
val_y = my_dlaplace(val_x)
lines(val_x, val_y, col = "blue", lwt=1)

# les deux courbes se superposent













#Exercice 4 : Chaine de Markov et Inférence Bayésienne
#*****************************************************

#Expression de la loi à priori de theta

rtheta=runif(100,0,1)   #tirage uniforme des valeures aléatoires
dtheta=dunif(rtheta)    #Expression de la densité
#Vérification de la cohérence graphique
truehist(rtheta, main="Lois à priori de Theta entre 0 et 1")
lines(rtheta, dtheta, col="red")

#Densité de la loi à Posteriori d_target

# La varibale qualitative (statut mutant) suit une loi Binomiale
# La vraissemblance P[X=70] sachant theta s'exprime donc en dbinom en focntion de dtheta
P_X70_theta=dbinom(dtheta,100,0.7)

# Expression de la densité de la loi cible
d_target = function(x) {     #x ici est notre variable à estimer c'est à dire pi[theta]
  return(x*P_X70_theta/0.7)    # la densité est le produit de pi[theta] et de la vraissemblance à 70 divisé par la probabilité de tirer un mutant
}
#Visualisation graphique de d_target
val_x = seq(from=-10, to=8, len=100)
val_y = d_target(val_x)
plot(val_x, val_y, col="red", type="l", main = "d_target pour P[X=70]")
     
#Choix de la loi de proposition : loi normale
rprop <- function(from_x) {
  return(rnorm(1, mean=from_x, sd=1))
}

dprop <- function(to_x, from_x) {
  return(dnorm(to_x, mean=from_x, sd=1))
}

#Expression de la fonction cible selon un algo MCMC
r_target <- function(n){
  x = (rep(NA,n)) #Pas d'allocation de mémoire
  x[1] = 5
  for (i in 1:(n-1)){
    yi = rprop(x[i]) #Pour chaque xi on propose une valeur yi suivant la loi de proposition (démarche pas à pas)
    p = min (c(1, d_target(yi)/d_target(x[i])) * (dprop(x[i], yi)/dprop(yi, x[i])))
    u = runif(1)
    if(u <= p){
      x[i+1] = yi    #On accepte ou on rejette yi selon une règle d'acceptation/rejet.
    }else{
      x[i+1] = x[i]
    }
  }  
  return(x)         
}

n = 10000
x = r_target(n)

truehist(x, col="gray", main="r_target selon un algorithme MCMC")
val_x = seq(from=-10, to=10, len=1000)
val_y = d_target(val_x)
lines(val_x, val_y, col="red", lwd=2)
#d_target et r_target ne se supperposent pas. I doit y avoir un souci dans notre code.

#Malheureusement nous ne sommes pas allées plus loin pour cet exercice par manque de compréhension des notions bayésiennes.