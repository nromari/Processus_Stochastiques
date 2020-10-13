#Simualtions en Biologie - Méthodes stochastiques
#-------------------------------------------------
#------------------------------------------------

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
truehist(r_X(100000), main="Monte Carlo suivant l'importance - Lois beta de paramètres (2,1)")
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
truehist(r_X_V2(100000), main="Monte Carlo suivant l'importance - Lois beta de paramètres (2.7,1.5)")
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
n=100
n_rep=1000000
val_x = seq(from=0, to=4, len=1000)

#Monte Carlo par tirage blanc ou noir (MC_BN)
res_MC_BN = rep(NA, n_rep)
for (i in 1:n_rep) {
  res_MC_BN[i] = MC_BN(n)
}

truehist(xlim=range(I_True), res_MC_BN, col = "gray", main="Intégration de g estimée par Monte Carlo tirage blanc/noir")
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

truehist(xlim=range(I_True), res_MC_S, col = "gray", main="Intégration de g estimée par Monte Carlo echantillonnage simple")
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

truehist(xlim=range(I_True), res_MC_I_1, col = "gray", main="Intégration de g estimée par Monte Carlo echantillonnage d'importance V1")
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

truehist(xlim=range(I_True), res_MC_I_2, col = "gray", main="Intégration de g estimée par Monte Carlo echantillonnage d'importance V2")
abline(v=I_True, lty=2, col="red")

mse_MC_I_2 = round(mean((res_MC_I_2-I_True)^2),2)
moy_MC_I_2 = round(mean(res_MC_I_2),2)
print(moy_MC_I_2)
print(mse_MC_I_2)
#La methode de Monte Carlo par tirage blanc/noir appliquée à g a un mse de 0.03 et une moyenne de 2.55


