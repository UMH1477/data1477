# librerías
library(tidyverse)
library(gridExtra)
library(kableExtra)
library(markovchain)
library(simmer)
library(simmer.plot)
library(simmer.bricks)
library(parallel)
library(queueing)


#####################################################################
## UNIDAD 1 CONCEPTOS BÁSICOS Y ESTIMACIÓN MC
#--------------------------------------------------------------------
# Estimación MC
#m=mean(iprob)
#error=sd(iprob)/sqrt(nsim)
#ic.low=round(m-qnorm(0.975)*error,3)
#ic.up=round(m+qnorm(0.975)*error,3)
#cat("Prob=",round(m,3),", IC95%=(",ic.low,",",ic.up,")")

#####################################################################
## UNIDAD 2 CMTD

#--------------------------------------------------------------------

# Definición de una CMTD
estados=as.character(0:1)
# Probabilidades de transición 
p1.matrix=matrix(c(0.51,0.49,0.9,0.1),byrow=TRUE,ncol=2,
                 dimnames=list(estados,estados))
# Procesos definidos por las probs.transición anteriores
idioma1 <- new("markovchain", states = estados,transitionMatrix = p1.matrix,byrow = TRUE)
# Simulación de nsim observaciones de dos textos en dos idiomas
nsim=1
texto1=rmarkovchain(nsim, idioma1)

#--------------------------------------------------------------------
# funciones en markovchain
# Probabilidades de transición**
#  - conditionalDistribution(proceso,t0)
#  - transitionProbability(proceso,t0,t1)

# Probabilidad de n pasos**
#   - (proceso)^n

# Tiempo esperado de primer paso**
#   - meanFirstPassageTime(proceso,destination)
#   - tiempo.pp(proceso,estado)

# Tiempo esperado de recurrencia**
#   - meanRecurrenceTime(proceso)

# Tiempo esperado de ocupación del proceso en n pasos**
#   - mocupa.proceso(proceso,n)

# Estado estacionario**
#   - verifyMarkovProperty(proceso)
#   - assessStationarity(proceso)
#   - period(proceso)
#   - steadyStates(proceso)

# Coste esperado total tras n pasos**
#   - expectedRewards(markovchain,n,rewards)

# Simulación**
#   - rmarkovchain(proceso,t0,include.t0=TRUE)

# Estimación matriz transición y compatibilidad**
#   - createSequenceMatrix(secuencia)
#   - markovchainFit(secuencia)
#   - verifyEmpiricalToTheoretical(secuencia,proceso)

#--------------------------------------------------------------------

# Función para calcular los tiempos de primer paso
tiempo.pp <- function(proceso, estado)
{
  # estados del proceso
  estados <- states(proceso)
  numestados <- length(estados)
  # posición de los estados deseados
  lestat <- length(estado)
  pos <- which(estados %in% estado)
  # matriz P_N
  P_N <- proceso[-pos,-pos]
  # vector de unos
  vector.1 <- matrix(rep(1, numestados-lestat), ncol=1)
  # sistema de ecuaciones
  sistema <- diag(numestados-lestat) - P_N
  # solución del sistema
  solucion <- solve(sistema, vector.1)
  return(solucion)
}

#--------------------------------------------------------------------
# Matriz de probs. transición de n pasos
ptran.n=function(ptran,n){
  # ptran es la matriz de transición de 1 paso
  # n son los pasos a dar
  i=1
  p=ptran
  while(i<n){
    p=p%*%ptran
    i=i+1
  }
  return(p)  
}
#--------------------------------------------------------------------
# Función para calcular los tiempos de ocupación
mocupa.proceso <- function(proceso, n)
{
  # Número de estados del proceso
  nestat <- dim(proceso)
  # Estados
  nombres<- names(proceso)
  # Generamos la matriz de ocupaciones
  mocupa <- diag(nestat)
  dimnames(mocupa) <- list(nombres, nombres)
  #  mocupa <- matrix(rep(0, nestat*nestat),
  #                 nrow = nestat, dimnames = list(nombres, nombres))
  # Bucle de cálculo de los tiempos de ocupación
  P=proceso[1:nestat,1:nestat]
  for (i in 1:n)
    mocupa <- mocupa + ptran.n(P,i)
  
  return(mocupa)
}
#--------------------------------------------------------------------
#####################################################################

## UNIDAD 3

#--------------------------------------------------------------------
# Superposición de PP
# Na ~ PP(lambda_a)
# Nb ~ PP(lambda_b)
# Nc ~ PP(lambda_a+lambda_c)
# Sin utilizar el resultado teórico, simulamos Na y Nb durante t minutos

simula.pp.superposicion = function(t,lambda_a,lambda_b,ta=NULL,tb=NULL){
  env=simmer()
  
  entrada_A =trajectory() %>%
    #log_("Entra a la autovía")%>%
    timeout(ta) #%>%
  #log_("Sale de la autovía")
  
  entrada_B =trajectory() %>%
    #log_("Entra a la autovía")%>%
    timeout(tb) #%>%
  #log_("Sale de la autovía")
  
  env %>%
    add_generator("entradaA",entrada_A,function() rexp(1,lambda_a)) %>%
    add_generator("entradaB",entrada_B,function() rexp(1,lambda_b)) %>%
    run(t)
}
#--------------------------------------------------------------------
# Adelgazamiento de PP
# N ~ PP(lambda)
# p y 1-p para optar por una bifurcación u otra

simula.pp.adelgazamiento = function(t,lambda,p){
  # p = probabilidad de ir a la izquierda
  env=simmer()
  
  entrada =trajectory() %>%
    branch(option=function() runif(1)<p, continue=FALSE, #si TRUE a la izqda
           trajectory() %>% log_("Tuerce a la izquierda") %>% set_attribute("giro",1)) %>% 
    # en otro caso a la derecha
    log_("Tuerce a la derecha") %>%
    set_attribute("giro",2)
  
  
  env %>%
    add_generator("entrada",entrada,function() rexp(1,lambda)) %>%
    run(t)
}
#t=10;lambda=2;p=0.5
#sim=simula.pp.adelgazamiento(t,lambda,p)
#get_mon_arrivals(sim)
#--------------------------------------------------------------------
#####################################################################

## UNIDAD 4 CMTC

# matriz.prob.trans(Rmat, ts, cal)
# tiempos.ocupacion(Rmat, Ts, cal)
# distr.lim.general(Rmat)
# tiempos.primer.paso(Rmat, A, estados) 

#--------------------------------------------------------------------
# Matriz de probabilidades de transición
matriz.prob.trans<- function(Rmat, ts, cal)
{
  # Algoritmo de uniformización para obtener P(t)
  
  
  # Parámetros de la función
  # Rmat: matriz de tasas 
  # ts: instante de tiempo
  # epsilon: error en la aproximación
  # cal: forma de obtener r, con dos valores 1 = máximo, 2 = suma finita
  epsilon <- 1e-05
  # Paso 2. Cálculo de r
  ris <- apply(Rmat, 1, sum)
  rlimit=ifelse(cal == 1, max(ris), sum(Rmat))
  # Paso 3. Cálculo de hat(P)
  hatP <- Rmat/rlimit
  diag(hatP) <- 1 - ris/rlimit
  # Paso 4. Cálculo de matrices y vectores accesorios
  rts <- rlimit*ts
  A <- hatP
  c <- exp(-rts)
  B <- c*diag(nrow(Rmat))
  suma <- c
  k <- 1
  # Bucle simulación
  while(suma < (1- epsilon))
  {
    c <- c*rts/k
    B <- B + c*A
    A <- A%*%hatP
    suma <- suma + c
    k <- k + 1
  }
  return(round(B, 4))
}
#--------------------------------------------------------------------
# Matriz de los tiempos de ocupación hasta Ts
tiempos.ocupacion<- function(Rmat, Ts, cal)
{
  # Algortimo de uniformización para obtener M(T)
  ################################################
  
  # Parámetros de la función
  # Rmat: matriz de tasas
  # ts: instante de tiempo
  # epsilon: error en la aproximación
  # cal: forms de obtener r con dos valores 1 = máximo, 2 = suma
  epsilon <- 1e-05
  # Paso 2. Calculo de r
  ris <- apply(Rmat, 1, sum)
  rlimit=ifelse(cal == 1, max(ris),sum(Rmat))
  # Paso 3. Calculo de hat(P)
  hatP <- Rmat/rlimit
  diag(hatP) <- 1 - ris/rlimit
  # Paso 4. 
  k <- 0
  A <- hatP
  # Paso 5.
  yek <- exp(-1*rlimit*Ts)
  ygk <- 1 - yek
  suma <- ygk
  # Paso 6.
  B <- ygk*diag(nrow(Rmat))
  # Bucle simulación
  cota <- Ts- epsilon
  while(suma/rlimit < cota)
  {
    k <- k + 1
    yek <- yek*(rlimit*Ts)/k
    ygk <- ygk - yek
    B <- B + ygk*A
    A <- A%*%hatP
    suma <- suma + ygk
  }
  return(round(B/rlimit, 4))
}
#--------------------------------------------------------------------
# Distribución límite CMTC
# Función para la resolución numérica de las ecuaciones de balance
distr.lim.general<-function(Rmat)
{
  # Parámetros de la función
  #=========================
  # Rmat: matriz de tasas del sistema
  
  # número de estados del sistema
  estados <- nrow(Rmat)
  # Calculamos r_i y lo colocamos en formato matriz
  sumarows <- diag(apply(Rmat, 1, sum), estados)
  # Matriz de coeficientes del sistema de ecuaciones de balance
  A <- t(R)-sumarows
  # Completamos la matriz añadiendo la restricción de suma de p`s igual a 1
  A <- rbind(A, rep(1, estados))
  # Vector de términos independientes del sistema
  CS <- c(rep(0, estados), 1)
  # Resolución del sistema
  ps <- qr.solve(A, CS)
  return(ps)
}
#--------------------------------------------------------------------
# Distribución límite para procesos de nacimiento-muerte
# Obtención de distribuciones límite 
distr.lim.nm <- function(estados, lambdas, mus)
{
  # Parámetros de la función
  # ========================
  # estados: número de estados del sistema (K)
  # lambdas: vector de tasas de nacimiento lambda0,..., lambda(K-1)
  # mus: vector de tasas de muerte mu1,...,muK
  
  # definimos vector de rho
  rhos <- rep(1, estados)
  # calculamos productos acumulados para lambda y mu
  prl <- cumprod(lambdas)
  prm <- cumprod(mus)
  # rellenamos rho con los productos acumulados
  rhos[2:estados] <- prl/prm
  # suma de rhos
  sumarhos <- sum(rhos)
  # vector de probabilidades
  ps <- rhos/sumarhos
  return(ps)
}
#--------------------------------------------------------------------
# Función para la resolución numérica de los tiempos de primer paso en CMTC
tiempos.primer.paso<-function(Rmat, A, estados)
{
  # Parámetros de la función
  #=========================
  # Rmat: matriz de tasas del sistema
  # A: vector de estados que queremos alcanzar
  # estados: conjunto de estados total (como carácter)
  
  # Estados como texto
  
  estados <- as.character(estados)
  # Tasas r
  rts <- diag(apply(Rmat, 1, sum), nrow(Rmat))
  # Seleccionamos el subconjunto de la matriz quitando fila y columna
  Rmod <- Rmat[-A, -A]
  rates <- rts[-A, -A]
  # Número de m´s a estimar
  lms <- nrow(Rmod)
  # Matriz de coeficientes del sistema de ecuaciones de balance
  B <- rates - Rmod
  # Vector de términos independientes del sistema
  CS <- rep(1, lms)
  # Resolución del sistema
  ps <- as.data.frame(qr.solve(B, CS))
  rownames(ps) <- paste("estado", estados)[-A]
  colnames(ps) <-"tiempo"
  return(ps)
}
#--------------------------------------------------------------------

# Simulación y Estimación con simmer
# nreplicas=2
# envs <- mclapply(1:nreplicas, function(i){
#  sistema(pars)%>%
#    wrap()},mc.set.seed=FALSE)

# almacenamos análisis de llegadas del sistema
# sim.arr<-as_tibble(get_mon_arrivals(envs))
# sim.arr %>% 
# group_by(replication) %>%
# summarise(m=activity_time)
#--------------------------------------------------------------------
#####################################################################

## UNIDAD 5

#--------------------------------------------------------------------
# Definición del entorno

# cola.in <- NewInput.MM1(lambda = 1, mu = 1, n = 4)
# NewInput.MM1K(lambda = lambda, mu = mu, k = k)
# NewInput.MMC(lambda = lambda, mu = muk, c = 2, n = 2)
# NewInput.MMCK(lambda = 6, mu = 4, c = 2, k = 5)

# cola=QueueingModel(cola.in)


# summary(s.MM1)
# RO, P0, Lq, Wq, L, W, Lqq, Wqq, Pn, Qn, FW, FWq
#--------------------------------------------------------------------


# FUNCIONES SIMULACIÓN COLAS
#################################################
# Sistema MM1
#################################################
cola.MM1 <- function(t, lambda, mu)
{
  # lambda: tasa de llegadas
  # mu: tasa de servicio
  
  # Funciones de tiempos
  tarrival <- function() rexp(1, lambda)
  tserver <- function() rexp(1, mu)
  
  # Trayectoria de servicio
  servicio <- trajectory() %>%
    visit("server", tserver)               
  
  # Entorno del sistema 
  #################################################
  simmer() %>%
    add_resource("server", capacity = 1, queue_size = Inf) %>%           
    add_generator("arrival", servicio, tarrival) %>% 
    run(until = t)     
}
#################################################
# Sistema M/M/1/K
#################################################
cola.MM1K <- function(t, lambda, mu, K)
{
  # lambda: tasa de llegadas
  # mu: tasa de servicio
  # K: capacidad del sistema
  
  # Funciones de tiempos
  tarrival <- function() rexp(1, lambda)
  tserver <- function() rexp(1, mu)
  
  # Tamaño de la cola
  qsize <- K - 1
  
  # Trayectoria de servicio
  servicio <- trajectory() %>%
    visit("server", tserver)               
  
  # Entorno del sistema 
  #################################################
  simmer() %>%
    add_resource("server", capacity = 1, queue_size = qsize) %>%           
    add_generator("arrival", servicio, tarrival) %>% 
    run(until = t)     
}
#################################################
# Sistema M/M/s
#################################################
cola.MMs <- function(t, lambda, mu, s)
{
  # lambda: tasa de llegadas
  # mu: tasa de servicio
  # s: servidores idénticos disponibles
  
  # Funciones de tiempos
  tarrival <- function() rexp(1, lambda)
  tserver <- function() rexp(1, mu)
  
  # Trayectoria de servicio
  servicio <- trajectory() %>%
    visit("server", tserver)               
  
  # Entorno del sistema 
  #################################################
  simmer() %>%
    add_resource("server", capacity = s, queue_size = Inf) %>%           
    add_generator("arrival", servicio, tarrival) %>% 
    run(until = t)     
}

#################################################
# Sistema M/M/s/K
#################################################
cola.MMsK <- function(t, lambda, mu, s, K)
{
  # lambda: tasa de llegadas
  # mu: tasa de servicio
  # s: servidores
  # K: capacidad del sistema
  
  # Funciones de tiempos
  tarrival <- function() rexp(1, lambda)
  tserver <- function() rexp(1, mu)
  
  # Tamaño de la cola
  qsize <- K - s
  
  # Trayectoria de servicio
  servicio <- trajectory() %>%
    visit("server", tserver)               
  
  # Entorno del sistema 
  #################################################
  simmer() %>%
    add_resource("server", capacity = s, queue_size = qsize) %>%           
    add_generator("arrival", servicio, tarrival) %>% 
    run(until = t)     
}
#####################################################################