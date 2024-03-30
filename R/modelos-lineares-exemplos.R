# Modelos Lineares - Exemplos R
# Material complementar: <https://gustavojy.github.io/mod-lin-r/>


# 1. Covariância, Correlação e Distância de Mahalanobis -------------------

## Vetores
y1 <- c(72,60,56,41,32,30,39,42,37,33,32,63,54,47,91,56,79,81,78,46,39,32,60,35,39,50,43,48)
y1

y2 <- c(66,53,57,29,32,35,39,43,40,29,30,45,46,51,79,68,65,80,55,38,35,30,50,37,36,34,37,54)
y2

y3 <- c(76,66,64,36,35,34,31,31,31,27,34,74,60,52,100,47,70,68,67,37,34,30,67,48,39,37,39,57)
y3

y4 <- c(77,63,58,38,36,26,27,25,25,36,28,63,52,43,75,50,61,58,60,38,37,32,54,39,31,40,50,43)
y4

## Matriz
Y <- cbind(y1, y2, y3, y4)
colnames(Y) <- c("North", "East", "South", "West")
Y

class(Y)

## Matriz de Variâncias e Covariâncias
### Opção 1
n <- nrow(Y)
p <- ncol(Y)
n; p

In <- diag(n)
In

jn <- matrix(data = 1, nrow = n, ncol = 1)
jn

Jnn <- matrix(data = 1, nrow = n, ncol = n)
Jnn

Sigma <- (1 / (n - 1)) * t(Y) %*% (In - (1/n) * Jnn) %*% Y
Sigma

### Opção 2
cov(Y)

## Matriz de correlação
### Opção 1
D <- sqrt(diag(Sigma))
D

corr <- solve(diag(D)) %*% Sigma %*% solve(diag(D))
corr

### Opção 2
cor(Y)

## Matriz de correlação retorna para Matriz de Variâncias e Covariâncias
Verifica <- diag(D) %*% corr %*% diag(D)
Verifica

## Distância de Mahalanobis (Distância padronizada)
### Opção 1
mi <- (1/n) * t(jn) %*% Y
mi

DM2 <- rep(0, n)
for (i in 1:n) {
  yi <- Y[i,]
  DM <- as.numeric((yi - mi) %*% solve(Sigma) %*% t(yi - mi))
  DM2[i] <- DM
}
DM2

### Opção 2
mahalanobis(x = Y, center = mi, cov = Sigma)

### Rank
rank <- rank(DM2)
data.frame(Y, DM2, rank)


# 2. Gráfico Normal Bivariada ------------------------------------------------

## Biblioteca
install.packages("plotly")
library(plotly)

## Vetores
y1 <- seq(from = -4, to = 4, by = 0.1)
y1

y2 <- seq(from = -4, to = 4, by = 0.1)
y2

## Coeficiente de correlação
r <- -0.75  # Modificar este valor para verificar as alterações gráficas

## Densidade bivariada da Normal

pi

z <- outer(y1, y2, function(y1,y2) { 
  phi <- 1/(2 * pi * sqrt(1 - r^2)) * exp(-(y1^2 - 2 * r * y1 * y2 + y2^2) / (2 * (1 - r^2)))
  return(phi)
})
z

## Gráfico Dinâmico 3D

plot_ly(x = y1, y = y2, z = z, type = "surface") |> 
  layout(title = paste("Densidade Normal Bivariada (r =", r, ")"))

## Link para aplicativo: <https://gustavojy.github.io/n-multi-graph-app/>.


# 3. Normal Multivariada - Propriedades --------------------------------------

## Exemplo 1 (4.4(a))----

mi <- c(3, 1, 2)
mi

Sigma <- matrix(c(4, 0, 2, 0, 1, -1, 2, -1, 3), nrow = 3, ncol = 3)
Sigma

### Teorema 1.1 (4.4A. i)
a <- c(1, -2, 1)
a

E_z <- t(a) %*% mi
E_z

var_z <- t(a) %*% Sigma %*% a
var_z

### Teorema 1.2 (4.4A. ii)
A <- matrix(c(1, -1, 1, 3, 1, -2), nrow = 2, ncol = 3, byrow = TRUE)
A

E_z1z2 <- A %*% mi
E_z1z2

cov_z1z2 <- A %*% Sigma %*% t(A)
cov_z1z2

### Teorema 2 (4.4B.)
#### Matriz A12
A12 <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)
A12

mi_12 <- A12 %*% mi
mi_12

Sigma_12 <- A12 %*% Sigma %*% t(A12)
Sigma_12

#### Matriz A13
A13 <- matrix(c(1, 0, 0, 0, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
A13

mi_13 <- A13 %*% mi; mi_13
mi_13

Sigma_13 <- A13 %*% Sigma %*% t(A13); Sigma_13
Sigma_13

### Teorema 3 (4.4C.)
a1 <- as.matrix(c(1, 0, 0))
a1

b2 <- as.matrix(c(0, 1, 0))
b2

cov_12 <- t(a1) %*% Sigma %*% b2
cov_12


## Exemplo 2 (4.4(b))----

### Teorema 4 (4.4D.) {#sec-T4}

mi <- c(2, 5, -2, 1)
mi

Sigma <- matrix(
  c(9, 0, 3, 3, 0, 1, -1, 2, 3, -1, 6, -3, 3, 2, -3, 7), 
  nrow = 4, byrow = TRUE
)
Sigma

Ay <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
Ay

mi_y <- Ay %*% mi
mi_y

Sigma_yy <- Ay %*% Sigma %*% t(Ay)
Sigma_yy

Ax <- matrix(c(0, 0, 1, 0, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
Ax

mi_x <- Ax %*% mi
mi_x

Sigma_xx <- Ax %*%Sigma %*% t(Ax)
Sigma_xx

Sigma_yx <- Ay %*% Sigma %*% t(Ax)
Sigma_yx

cov_ydx <- Sigma_yy - Sigma_yx %*% solve(Sigma_xx) %*% t(Sigma_yx)
cov_ydx

Sigma_yy
cov_ydx

## Exemplo 3 (4.4(c))----

### Teorema 4 (4.4D.) - Corolário 1

mi <- c(2, 5, -2, 1)
mi

Sigma <- matrix(
  c(9, 0, 3, 3, 0, 1, -1, 2, 3, -1, 6, -3, 3, 2, -3, 7), 
  nrow = 4, byrow = TRUE
)
Sigma

Ay <- matrix(c(1, 0, 0, 0), nrow = 1)
Ay

mi_y <- Ay %*% mi
mi_y

Sigma_yy <- Ay %*% Sigma %*% t(Ay)
Sigma_yy

Ax <- matrix(
  c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 
  nrow = 3, ncol = 4, byrow = TRUE
)
Ax

mi_x <- Ax %*% mi
mi_x

Sigma_xx <- Ax %*% Sigma %*% t(Ax)
Sigma_xx

sigma_yx <- Ay %*% Sigma %*% t(Ax)
sigma_yx

cov_ydx <- Sigma_yy - sigma_yx %*% solve(Sigma_xx) %*% t(sigma_yx)
cov_ydx

Sigma_yy
cov_ydx


# 4. Correlação Parcial ------------------------------------------------------

Sigma <- matrix(
  c(9, 0, 3, 3, 0, 1, -1, 2, 3, -1, 6, -3, 3, 2, -3, 7), 
  nrow = 4, byrow = TRUE
)
Sigma

D <- sqrt(diag(Sigma))
D

Ro <- solve(diag(D)) %*% Sigma %*% solve(diag(D))
Ro

Ay <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2, byrow = T)
Ay

Sigma_yy <- Ay %*% Sigma %*% t(Ay)
Sigma_yy

Dyy <- sqrt(diag(Sigma_yy))
Dyy

Ro_yy <- solve(diag(Dyy)) %*% Sigma_yy %*% solve(diag(Dyy))
Ro_yy

Ax <- matrix(c(0, 0, 1, 0, 0, 0, 0, 1), nrow = 2, byrow = T)
Ax

Sigma_xx <- Ax %*% Sigma %*% t(Ax)
Sigma_xx

Sigma_yx <- Ay %*% Sigma %*% t(Ax)
Sigma_yx

cov_ydx <- Sigma_yy - Sigma_yx %*% solve(Sigma_xx) %*% t(Sigma_yx)
cov_ydx

D <- sqrt(diag(cov_ydx))
D

Ro_ydx <- solve(diag(D)) %*% cov_ydx %*% solve(diag(D))
Ro_ydx

Ro_yy
Ro_ydx



# 5. Regressão Linear Simples -----------------------------------------------

## Vetores
y <- c(95,80,0,0,79,77,72,66,98,90,0,95,35,50,72,55,75,66)

x1 <- c(96,77,0,0,78,64,89,47,90,93,18,86,0,30,59,77,74,67)

## Análise exploratória
notas <- data.frame(y, x1)
notas

library(ggplot2)
ggplot(data = notas, aes(x = x1, y = y)) +
  geom_point()

## Estimação dos parâmetros beta
### Opção 1
n <- length(y)
n

jn <- matrix(data = 1, nrow = n, ncol = 1)
jn

Jnn <- jn %*% t(jn)
Jnn

In <- diag(n)
In

Beta1 <- as.numeric(t(x1) %*% (In - (1 / n) * Jnn) %*% y / (t(x1) %*% (In - (1 / n) * Jnn) %*% x1))
Beta1

Beta0 <- as.numeric((1 / n) * t(jn) %*% (y - Beta1 * x1))
Beta0

### Opção 2
lm(y ~ x1)

### Gráfico
#### Opção 1
ggplot(data = notas, aes(x = x1, y = y)) +
  geom_point() + 
  geom_abline(intercept = 10.7269, slope = 0.8726, color = "red")

#### Opção 2

ggplot(data = notas, aes(x = x1, y = y)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "red", se = FALSE)

## Estimação da variância

### SQRes
y_hat <- Beta0 + Beta1 * x1
y_hat

Res <- y - y_hat
Res

SQRes <- t(Res) %*% Res
SQRes

### Variância s2

s2 <- as.numeric(SQRes / (n - 2))
s2

### Resíduo padronizado
res_pad <- Res / sqrt(s2)
res_pad

### Variância dos dados originais

var_y <- (t(y) %*% (In - (1/n) * Jnn) %*% y) / (n - 1)
var_y

### Desvio padrão
s <- sqrt(s2)
s

## Teste de Hipóteses e Intervalo de Confiança para Beta1

### Variância e Erro Padrão
x_barra <- t(jn) %*% (x1 / n)
x_barra

var_Beta0 <- s2 * ((1 / n) + (x_barra^2 / (t(x1) %*% (In - (1 / n) * Jnn) %*% x1)))
var_Beta0

stderr_Beta0 <- sqrt(var_Beta0)
stderr_Beta0

var_Beta1 <- s2 / (t(x1) %*% (In - (1 / n) * Jnn) %*% x1)
var_Beta1

stderr_Beta1 <- sqrt(var_Beta1)
stderr_Beta1

### Valor da estatística t
t0 <- Beta0 / stderr_Beta0
t0

t1 <- Beta1 / stderr_Beta1
t1

### Intervalo de confiança
ttab <- qt(0.975, n - 2) # quantil t com probabilidade 0.975 e gl = n-2

liminf0 <- Beta0 - ttab * stderr_Beta0
liminf0

limsup0 <- Beta0 + ttab * stderr_Beta0
limsup0

liminf1 <- Beta1 - ttab * stderr_Beta1
liminf1

limsup1 <- Beta1 + ttab * stderr_Beta1
limsup1

### Coeficiente de Determinação

#### SQReg 
y_barra <- (1 / n) * Jnn %*% y
y_barra

Reg <- y_hat - y_barra
Reg

SQReg <- t(Reg) %*% Reg
SQReg

#### SQTotal
Tot <- y - y_barra
Tot

SQTotal <- t(Tot) %*% Tot
SQTotal

SQTotal <- SQReg + SQRes
SQTotal

#### Coef. determinação (r2)
R2 <- SQReg / SQTotal
R2

#### Coef de correlação (r)
r <- sqrt(R2)
r

#### Estatística t usando valor de r
tcalc1 <- r * sqrt(n - 2) / (sqrt(1 - r^2))
tcalc1

tcalc2 <- Beta1 / stderr_Beta1
tcalc2

### p-valor
p_valor <- 2 * (1 - pt(abs(tcalc1), n - 2))
p_valor

## Visão geral
modelo <- lm(y ~ x1)
summary(modelo)
