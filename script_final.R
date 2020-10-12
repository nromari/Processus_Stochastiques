#Simualtions en Biologie - Méthodes stochastiques
#-------------------------------------------------
#------------------------------------------------

#Exercice 2 : Methodes Monte Carlo
#---------------------------------

#Il s'agit ici de comparer les methodes d'integration Monte Carlo.
#Les méthodes de Monte Carlo sont des modèles de simulation stochastiques
#Trois methodes seront evaluees : MC "tirage noir ou blanc", MC simple, et MC suivant l'importance

#1) La fonction à simuler g(x)
#.............................

#Afin de mettre en place efficacement les modeles de simulation il faut bien connaître la fonction.
#En effet pour simuler une variable aletaoire continue il est necessaire de calculer son integrale.
g = function(x) {                 #g est une fonction exponenetielle
  return((exp(x)-1)/(exp(1)-1))
}

#Visualisation de g sur un interval [0-2]
val_x = seq(from=0, to=2, len=1000)
val_y = g(val_x)
plot(val_x, val_y, type="l")      #L'allure de la courbe est bien exponentielle

#Verification de la valeur de l'integrale pour la valeur de max de valx (x=2)
I_True = (exp(2)-3)/(exp(1)-1)
print(I_True)
#La valeur réelle de g(2) est 2.554328
integrate(g, lower=0, upper=2)
# la valeur de l'integrale de g pour l'interval [0-2] est de 2.554328
#g(2) = integ(g)[0-2] --> On retrouve bien la même valeur

#2)


