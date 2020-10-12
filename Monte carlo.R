# Execice2 : Efficacité des méthodes d'intégration Monte carlo
#-------------------------------------------------------------
#-------------------------------------------------------------

# Définition de la fonction
#--------------------------

g = function(x) {
  return((exp(x)-1)/(exp(1)-1))
}

val_x = seq(from=0, to=2, len=1000)
val_y = g(val_x)
plot(val_x, val_y, type="l")

I_True = (exp(2)-3)/(exp(1)-1) # oncalcul la valeur réelle
print(I_True)

integrate(g, lower=0, upper=2) # on retrouve la bonne valeur

# Méthode blanc et noir de 0 à 2
#-------------------------------

# Avant de faire la fonction, il faut choisir un majorant. Ici g(2) est le majorant optimale
# car sur l'interval choisit g(2) est le plus petit majorant possible

MC_BN = function(n) {
  m = g(2)
  u = runif(n, min = 0, max = 2)
  v = runif(n, min = 0, max = m)
  s = m*2
  ns = sum(v<g(u))
  return(ns/n*s)
}

MC_BN(1000) # valeurs prochent de l'attendu. A lancer plusieurs fois.

# Méthode simple (Espérance) de 0 à 2
#------------------------------------

MC_S = function(n) {
  x = runif(n, min=0, max=2)
  return(2*sum(g(x))/n) # le 2 correspond au dunif de la fonction. ici on peut le sortir
}

MC_S(1000) # résultat cohérent

# Méthode par l'importance
#-------------------------

# Créer une fonction qui permet de tirer et une fonction qui permet d'avoir la densité

# La fonction pour le tirage doit être du type bêta qui va échantilloner entre 0 et 1 (pas 0 2 car impossible avec loies de type bêta)
#Cf. wikipédia loi bêta bleue. pour avoir une bonne méthode MC par importance, il faut échantillonner avec une fonction qui ressemble à g 
#La solution pour faire du 0 2  c'est de tout multiplier par 2 or de la fonction

r_X = function(n) {
  return(2*rbeta(n, shape1 = 2, shape2 = 1))# On créé le Y=aX avec a=2
}
#A faire : un truehistr_x(10000)
#function de densité de probabilité de x
d_X = function(x) {
  return(1/2*dbeta(x/2, shape1 = 2, shape2 = 1))# ne pas oublier de diviser par 2 à cause de la transformation Y=aX car dbeta ne peut pas accepter des valeurs hor de 0:1
}
#A faire : un truehist de de r_X avec d_x pour vérifier que ça colle
#le but ici pour rappel est de construir l'échantillonage

I_importance_1 = function(n) {
  x = r_X(n)
  xx = g(x)/d_X(x)
  return(mean(xx))
}
I_importance_1(1000) #Les valeurs sont cohérentes

Truehist(r_X(1000000)) #!! Installer le bon package
val_x = seq(from=0, to=2, len=1000)
val_y = d_X(val(x))
lines(val_x, val_y, col=red)

#Evaluer les méthodes à l'aide du MSE
#------------------------------------

#Le but est de choisir le moins mauvais... ne pas prendre trop de point comme pour l'aiguille de buffon

mse = function(n,MC_BN,n_rep){
  vec = rep(MC_BN,n)
}






