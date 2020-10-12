#Simualtions en Biologie - Méthodes stochastiques
#-------------------------------------------------
#------------------------------------------------

#Exercice 2 : Methodes Monte Carlo
#*********************************

#Il s'agit ici de comparer les methodes d'integration Monte Carlo.
#Les méthodes de Monte Carlo sont des modèles de simulation stochastiques.
#Trois methodes seront evaluees : MC "tirage noir ou blanc", MC simple, et MC suivant l'importance.

#1) La fonction à simuler g(x)
#-----------------------------

#Afin de mettre en place efficacement les modeles de simulation il faut bien connaître la fonction.
#En effet pour simuler une variable aletaoire continue il est necessaire de calculer son integrale.

g = function(x) {                 #g est une fonction exponentielle.
  return((exp(x)-1)/(exp(1)-1))
}

#Visualisation de g sur un interval [0-2]
val_x = seq(from=0, to=2, len=1000)
val_y = g(val_x)
plot(val_x, val_y, type="l", main="g(x) sur l'interval [0-2]")      #L'allure de la courbe est bien exponentielle.

#Verification de la valeur de l'integrale pour la valeur de max de valx (x=2).
I_True = (exp(2)-3)/(exp(1)-1)
print(I_True)
#La valeur réelle de g(2) est 2.554328
integrate(g, lower=0, upper=2)
#la valeur de l'integrale de g pour l'interval [0-2] est de 2.554328
#g(2) = integ(g)[0-2] --> On retrouve bien la même valeur.

#2) Methode de Monte Carlo par tirage noir ou blanc sur l'interval [0-2]
#-----------------------------------------------------------------------
#Le principe est basé sur le tirage au hasard des points.
#Il faut donc définir la surface et échantillonner l'espace de manière uniforme
#Si l'intégrale correspond à un rectangle les X et Y seront tirés de manière uniforme
#La délimitation du rectangle passe par le choix d'un majorant

#Choix du majorant : g(2) est le majorant optimale car sur l'interval choisit car g(2) est le plus petit majorant possible.

MC_BN = function(n) {
  m = g(2)    #choix du majorant
  u = runif(n, min = 0, max = 2)    #Tirage uniforme dans l'intervale de g
  v = runif(n, min = 0, max = m)    #Tirage uniforme dans le rectangle
  s = m*2
  ns = sum(v<g(u))    #nombre de sucés = nombre de tirage en dessous de g(u)
  return(ns/n*s)
}

MC_BN(1000) #A lancer plusieurs fois
# Nous lançons 3 fois : 2.535868 , 2.751529, 2.587924
#Les valeurs obtenues sont proches de la valeur attendues (pour rappel valeur vraie de l'intégrale = 2.554328)

