#Transformation Loi exponenteielle (exos cours)

require(MASS)

n=10000
u=runif(n)# loi uniforme
truehist(u, col="gray")

a= 3.4
x=-log(1-u)/a
truehist(x, col="gray")# Loie en L

val_x = seq(from=0, to=3, len=1000)
val_y = dexp(val_x, rate = a) #fonction de densité d'une focntion exponentielle.
lines(val_x, val_y, col="red", lwd=2)

#refaire en partant d'une loie normale

u1=rnorm(n)# loi normale
truehist(u1, col="gray") # Pas d'allure en L ... on a pas la bonne loie !!
x1=-log(1-u1)/a
truehist(x1, col="gray")# Loie en L

# On peut aussi utiliser q exp pour l'expression de val_y mais attention il faut remplacer val_x par u.
val_y = qexp(u, rate = a)
truehist(x, col="gray")# Lois en L : on obtient bien la même chose.

# Exercice 3

#1) implémenter la fonction my_rpois = focntion de simulation d'une variable de Poisson de paramètre lamda

my_rpois_one <- function(lambda) { #Renvoie 1 valeure tirée au hasard dans une loie de poisson de param lamda
  u <- runif(1) # tire une valeure uniforme entre 0 et 1
  x <- 0 # définition de xinit
  while (ppois(x, lambda) < u) { # Appel de la loie de poisson en fonction de lambda
    x <- x+1 #itération
  }
  return(x)
}

my_rpois_one(10.5)

my_rpois <- function(n, lambda) { #répétition de my_rpois_one
  x <- rep(NA, len=n)
  for (i in 1:n) {
    x[i] <- my_rpois_one(lambda)
  }
  return(x)
}

lambda <- 6.4 #test
n <- 10000
x <- my_rpois(n, lambda)
head(x)

plot(table(x)/n) #représentation de la proportion
val_x <- seq(from=0, to=10) # vérification de la corrélation avec dpois
val_y <- dpois(val_x, lambda)
points(val_x, val_y, col="red")

#2) A FAIRE !!! cf vidéo

#3) Implémenter my_rlaplace = fonction de simulation d'une variable de Laplace par transformation (inversion de la fonction de répartition)

# Etape 1 : exprimer la fonction de Densité d’une Laplace "dlaplace" en sachant que :
#Si x < 0, f(x) = 1/2 exp(x)
#Si x > 0, f(x) = 1/2 exp(-x))

my_dlaplace <- function(x) { #On reprend ici la fonction g(x) donné dans l'énoncé
  return(1/2*exp(-abs(x)))
}

val_x <- seq(from=-4, to=4, len=100) #Représentation graphique de my_dlaplace
val_y <- my_dlaplace(val_x)
plot(val_x, val_y, col="red", type="l")


# Implémentatoin de la Fonction de répartition "plaplace" en sachant que :
#Si x < 0, F(x) = 1/2 exp(x)
#Si x > 0, F(x) = 1/2(2-exp(-x))

my_plaplace_one <- function(x) { #ici fonction qui renvoie un scalaire... puis faire une boucle pour avoir la fonction plapalce en vecteur, qui sera celle à inverser
  if (x < 0) {
    p <- 1/2*exp(x)
  } else {
    p <- 1/2*(2-exp(-x))
  }
  return(p)
}

my_plaplace <- function(x_vec) {
  n <- length(x_vec)
  p_vec <- rep(NA, n)
  for (i in 1:n) {
    p_vec[i] <- my_plaplace_one(x_vec[i])
  }
  return(p_vec)
}

my_plaplace(c(-0.2, -0.1, 0, 1, 3, 5)) #test
val_x <- seq(from=-4, to=4, len=100) # Visu graphique
val_y <- my_plaplace(val_x)
plot(val_x, val_y, col="red", type="l")

# Implémentation de la focntion Fonction quantile qlaplace
#Q(p) = … 
#Si p < 1/2
#Si p > 1/2

my_qlaplace_one <- function(p) {
  if (p < 1/2) {
    x <- log(2*p)
  } else {
    x <- -log(2-2*p)
  }
  return(x)
}

my_qlaplace <- function(p_vec) {
  n <- length(p_vec)
  x_vec <- rep(NA, n)
  for (i in 1:n) {
    x_vec[i] <- my_qlaplace_one(p_vec[i])
  }
  return(x_vec)
}

val_x <- seq(from=0, to=1, len=100)
val_y <- my_qlaplace(val_x)
plot(val_x, val_y, col="red", type="l")

# Implémentation de la focntion Fonction de simulation rlaplace

my_rlaplace <- function(n) {
  u <- runif(n)
  x <- my_qlaplace(u)
  return(x)
}

n <- 100000
x <- my_rlaplace(n)
head(x)

truehist(x)
val_x <- seq(from=-10, to=10, len=1000)
val_y <- my_dlaplace(val_x)
lines(val_x, val_y, col="red")


#4) Par méthode de rejet

h <- function(x) { #La fonction h va nous permettre de choisir un m optimale
  return(dnorm(x)/my_dlaplace(x))
}


val_x <- seq(from=-10, to=10, len=1000)
val_y <- h(val_x)
plot(val_x, val_y, col="red", type="l")

m <- sqrt(2*exp(1)/pi)
abline(h=m, lty=2)

val_x <- seq(from=-4, to=4, len=1000)
val_y <- dnorm(val_x)
plot(val_x, val_y, col="red", type="l", ylim=c(0, 0.6))
val_y <- m*my_dlaplace(val_x)
lines(val_x, val_y, col="blue")

my_rnorm_rejet <- function(n) {
  m <- sqrt(2*exp(1)/pi)
  x <- my_rlaplace(n)
  u <- runif(n)
  filter <- which(u <= dnorm(x)/(m*my_dlaplace(x)))
  return(x[filter])
}


n <- 100000
x <- my_rnorm_rejet(n)
head(x)

truehist(x, col="gray")
val_x <- seq(from=-4, to=4, len=1000)
val_y <- dnorm(val_x)
lines(val_x, val_y, col="red", lwd=2)

length(x)/n #On transforme 76% des données
1/m #donne le taux de transformation max atteignable : notre modele est donc optimal

#5) Implémentation d'une fonction de simulation par méthode MCMC

#Fonction de densité d_target
my_dtarget <- function(x) { #On reprend ici la fonction f(x) donné dans l'énoncé
  return(0.2*dnorm(x,-3,2)+0.5*dnorm(x,0,1)+0.3*dnorm(x,5,3))
}

val_x <- seq(from=-10, to=8, len=100) #Représentation graphique de my_dlaplace
val_y <- my_dtarget(val_x)
plot(val_x, val_y, col="red", type="l")

#
rprop <- function(from_x) {
  return(rnorm(1, mean=from_x, sd=1))
}

dprop <- function(to_x, from_x) {
  return(dnorm(to_x, mean=from_x, sd=1))
}

my_rtarget <- function(n){
  x = rep(NA,n) #Pas d'allocation de mémoire
  x[1] = 5
  for (i in 1:(n-1)){
    yi = rprop(x[i]) #Pour le point de départ soit on et yi
    p = min (c(1, my_dtarget(yi)/my_dtarget(x[i])) * (dprop(x[i], yi)/dprop(yi, x[i])))
            u = runif(1)
            if(u <= p){
              x[i+1] = yi
            }else{
              x[i+1] = xi
            }
  }  
  return(x)         
}

n = 10000
x = my_rtarget(n)

truehist(x, col="gray")
val_x <- seq(from=-10, to=10, len=1000)
val_y <- my_dtarget(val_x)
lines(val_x, val_y, col="red", lwd=2)



