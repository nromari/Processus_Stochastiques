# Exercice 3 Similation de variable aléatoires

require(MASS)

####### Loi de poisson #######

# 1-) Implémentation de la fonction my_rpois de simulation d'une variable
# de poisson en utilisant une méthode d'inversion de la fonction de repartition.

## Implémentation de my_rpois_one(lambda) : fonction qui renvoie une seule
## valeur tirée aléatoirement dans une loi de Poisson de paramètre Lambda
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

## Exemple

cat("valeur suivant une loi de poisson de parametre",lambda,":", my_rpois_one(lambda))

## Impléméntation de la fonction my_rpois(n,lambda) , qui fait appelle n fois la fonction my_rpois_one
## et renvoie un vecteur de n valeurs tirées aléatoirement et indépendamment dans une loi de Poisson de paramètre lambda.

my_rpois = function(n, lambda) {
  n_pois = rep(NA, len=n)# initialision d'un vecteur de n valeurs
  for (i in 1:n) {
    n_pois[i] = my_rpois_one(lambda) # affectation progression des éléments du vecteur X, suivant une loi de poisson
  }
  return(n_pois)
}

## Exemple


n = 100000
n_pois = my_rpois(n, lambda)
head(n_pois) # pour visualiser que nous avons bien des entiers, et qui sont aléatoires
plot(table(n_pois)/n, col="blue", main="Distribution de variables aléatoires \n suivant une loi de poisson", xlab="variables aléatoires x", ylab="Probabilités")
#require(MASS)
#truehist(table(n_pois)/n)
#truehist(n_pois)

# Verification de notre fonction de simulation my_rpois en utilisant la fonction rpois de R

val_x = seq(from=0, to=12) # entiers de 0 à 8
val_y = dpois(val_x, lambda) # valeurs suivant une loi de poisson
points(val_x, val_y, col="red") # representation des points sur notre courbe de my_rpois
# Notre simulation est bien en adéquation avec la loi de poisson de paramètre lambda. Les deux courbes (rouge et bleu ont la meme allure)


####### Loi discrète #######

# 2-) Implémentation d'une fonction my_rdiscret de simulatio n de la variable aléatoire réelle discrète X définie par une loi de probabilité X

#### Comme pour la fonction my_rpois_one, la fonction my_rdiscret_one() renvoie une seule valeur tirée aéatoitement dans la loi de X

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
### La fonction my_rdiscret renvoie n  valeurs tirées de la loi X

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

# Les deux courbes sont correlées


# 3-) Méthode de transformation

####### Loi de Laplace #######

# Implémentation d'une fonction my_rlaplace de simulation d'une variable de Laplace en utilisant une méthode d'inversion de la fonction de répartition.

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


####### Loi normale #######

# Utilisation de la fonction de distribution de laplace implementée plus haut
# pour implementer la fonction my_rnorm de similation (méthode de rejet)

# 4-) Méthode de rejet

## Plus petite valeur m telle f(x)<=m*g(x) où 
# f(x)= densité d'une variable normale centrée-reduite
# g(x) = densité d'une variable de laplace
# avec h(x)=f(x)/g(x)


h = function(x) {
  return(dnorm(x)/my_dlaplace(x))
}
# Exemple

val_x = seq(from=-10, to=10, len=1000)
val_y = h(val_x)
plot(val_x, val_y, col="red", type="l", main="Courbe de la fonction h(x)")

m = sqrt(2*exp(1)/pi)

abline(h = m, lty=2)

# Nous voyons que la droite y = m passe par le maximum de la courbe en rouge. ce qui nous permet de dire que 
# que m = sqrt(2*exp(1)/pi), est la plus petite valeur qui permet de verifier la condition f(x)<=m*g(x) quelque soit x réel.

# Verification de la condition avec la valeur de m
val_x = seq(from=-4, to=4, len=1000) # vecteur de 1000 valeurs entre -4 et 4
val_y = dnorm(val_x) # Distribution de la loi normale des 1000 valeurs du vecteurs (fonction f(x))
plot(val_x, val_y, col="blue", type="l", ylim=c(0, 0.7), main="Verification de la condition f(x) <= m * g(x)", xlab="variables alétoires X", ylab="Probabilité")
val_y = m*my_dlaplace(val_x) # fonction g(x)
lines(val_x, val_y, col="red")

# Nous voyons bien que pour x réel la fonction rouge (g(x)) est au dessus de celle bleue (f(x))

### Implémentation de la fonction my_rnorm_rejet(n)

my_rnorm_rejet = function(n) {
  m = sqrt(2*exp(1)/pi)
  x = my_rlaplace(n)
  u = runif(n)
  filter = which(u <= dnorm(x)/(m*my_dlaplace(x)))
  return(x[filter])
}

n = 100000
x = my_rnorm_rejet(n)
head(x)

truehist(x, col="yellow", main="Loi normale", xlab="variables aléatoires x", ylab="Probabilités")
val_x = seq(from = -4, to = 4, len = 1000)
val_y = dnorm(val_x)
lines(val_x, val_y, col = "red", lwd = 2)

rejet_1 = length(x)/n
rejet_2 = 1/m


# les deux valeurs sont très très proches. Nous pouvons donc dire notre implémentation a bien le taux de rejet attendu

# 5-) Algorithme MCMC

g1 = function(x){
  return(dnorm(x,-3,2))
}

g2 = function(x){
  return(dnorm(x, 0, 1))
}

g3 = function(x){
  return(dnorm(x,5,3))
}

my_dtarget = function(x){
  return(0.2*g1(x)+0.5*g2(x)+0.3*g3(x))
}


val_x <- seq(from=-10, to=8, len=100)
val_y <- my_dtarget(val_x)
plot(val_x, val_y, col="red", type="l", main="Représentation graphique de my_dtarget")

## Loi de proposition

# A-) Loi uniforme
my_rprop = function(x, delta){
  return(runif(1, x-delta, x+delta))
}

my_rtarget = function(n){
  valX = rep(NA, n)
  valX[1] = -10
  for(i in 2:n){
    valY = my_rprop(valX[i-1], 1)
    seuil = min(1, (my_dtarget(valY)/my_dtarget(valX[i-1])))
    U = runif(1)
    if( U <= seuil){
      valX[i] = valY
    }
    else{
      valX[i] = valX[i-1]
    }
  }
  return(valX)
}

n = 100000
res_target = my_rtarget(n)
truehist(res_target, col="yellow", main= "Distribution de variables aléatoires aléatoire continue \n par simulation MCMC: Loi uniforme", xlab="variables aléatoires X", ylab="Probabilités" )

valX = seq(-10, 12, len = 1000)
valY = my_dtarget(valX)
lines(valX, valY, type="l", col="red")
## B-) Loi normale

rprop <- function(from_x) {
  return(rnorm(1, mean=from_x, sd=1))
}

dprop <- function(to_x, from_x) {
  return(dnorm(to_x, mean=from_x, sd=1))
}

my_rtarget <- function(n){
  x = rep(NA,n)
  x[1]=5
  for (i in 1:(n-1)){
    yi = rprop(x[i]) 
    p = min (c(1, my_dtarget(yi)/my_dtarget(x[i])) * (dprop(x[i], yi)/dprop(yi, x[i])))
    u = runif(1)
    if(u <= p){
      x[i+1] = yi
    }else{
      x[i+1] = x[i]
    }
  }  
  return(x)         
}

n = 100000
x = my_rtarget(n)
head(x)
truehist(x, col="yellow",main= "Distribution de variables aléatoires aléatoire continue \n par simulation MCMC: Loi normale", xlab="variables aléatoires X", ylab="Probabilités" )
val_x <- seq(from=-10, to=15, len=1000)
val_y <- my_dtarget(val_x)
lines(val_x, val_y, col="red", lwd=2)

summary(x)
