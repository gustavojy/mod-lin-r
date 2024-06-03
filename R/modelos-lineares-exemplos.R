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
# install.packages("plotly")
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

# 6. Regressão Linear Múltipla - Estimação ----------------------------------

## Dados
y <- c(2, 3, 2, 7, 6, 8, 10, 7, 8, 12, 11, 14)

x1 <- c(0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8)

x2 <- c(2, 6, 7, 5, 9, 8, 7, 10, 11, 9, 15, 13)

## Estimação dos parâmetros beta

n <- length(y)
n

jn <- matrix(data = 1, nrow = n, ncol = 1)
jn

Jnn <- jn %*% t(jn)
Jnn

In <- diag(n)
In

x0 <- jn
colnames(x0) <- ("x0")
x0

### x1 sozinho
X01 <- cbind(x0, x1)
X01

Beta01 <- solve(t(X01) %*% X01) %*% t(X01) %*% y
Beta01

### x2 sozinho
X02 <- cbind(x0, x2)
X02

Beta02 <- solve(t(X02) %*% X02) %*% t(X02) %*% y
Beta02

### x1 e x2 em conjunto
X012 <- cbind(x0, x1, x2)
X012

Beta012 <- solve(t(X012) %*% X012) %*% t(X012) %*% y
Beta012

## Estimação da variância

k <- 2

### Opção 1
s2 <- as.numeric((t(y) %*% y) - (t(Beta012) %*% t(X012) %*% y)) / (n - k - 1)
s2

### Opção 2
y_hat <- X012 %*% Beta012
colnames(y_hat) <- ("y_hat")

res <- y - y_hat

SQRes <- t(res) %*% res

s2 <- as.numeric(SQRes / (n - k - 1))

### cov(beta)
cov_Beta <- s2 * solve(t(X012) %*% X012)
cov_Beta

## Modelo de Regressão na Forma Centrada
x1x2 <- cbind(x1, x2)
x1x2

x1x2c <- (In - (1 / n) * Jnn) %*% x1x2
x1x2c

Xc <- cbind(x0, x1x2c)
Xc

Betac <- solve(t(Xc) %*% Xc) %*% t(Xc) %*% y
Betac

### Variância centrada
#### Opção 1
s2c <- as.numeric((t(y) %*% y) - (t(Betac) %*% t(Xc) %*% y)) / (n - k - 1)
s2c

#### Opção 2
y_hatc <- Xc %*% Betac
colnames(y_hatc) <- ("y_hatc")

res_c <- y - y_hatc

SQRes_c <- t(res_c) %*% res_c

s2c <- as.numeric(SQRes / (n - k - 1))

cov_Betac <- s2c * solve(t(X012) %*% X012)
cov_Betac

## Coeficiente de Determinação na Regressão com x-fixos
SQReg <- t(y) %*% (Xc %*% solve(t(Xc) %*% Xc) %*% t(Xc) - (1/n) * Jnn) %*% y
SQReg

SQTot <- t(y) %*% (In - (1/n) * Jnn) %*% y
SQTot

R2 <- SQReg / SQTot
R2

R2aj <- ((n - 1) * R2 - k) / (n - k - 1)
R2aj

## Resumo dos resultados
### Betas
data.frame(
  beta = Beta012,
  betac = Betac
) |> 
  kableExtra::kbl(
    align = "c",
    col.names = c("Beta", "Beta Centrado")
  ) |> 
  kableExtra::kable_styling(position = "center")

### s2
data.frame(
  s2 = s2,
  s2c = s2c
) |> 
  kableExtra::kbl(
    align = "c",
    col.names = c("s2", "s2 centrado")
  ) |> 
  kableExtra::kable_styling(position = "center")

### cov(Beta)
cat("Cov_Beta:\n")
print(cov_Beta, digits = 4)

cat("Cov_Beta Centrado:\n")
print(cov_Betac, digits = 4)

### y_hat
data.frame(
  y_hat = X012 %*% Beta012 |> round(3),
  y_hatc = Xc %*% Betac |> round(3)
) |> 
  kableExtra::kbl(
    align = "c",
    col.names = c("y_hat", "y_hat centrado")
  ) |> 
  kableExtra::kable_styling(position = "center")



# 7.1 Regressão Linear Múltipla - Teste Hip ------------------------------

## Dados

y <- c(2, 3, 2, 7, 6, 8, 10, 7, 8, 12, 11, 14)

x0 <- matrix(data = 1, nrow = length(y), ncol = 1)
colnames(x0) <- ("x0")

x1 <- c(0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8)

x2 <- c(2, 6, 7, 5, 9, 8, 7, 10, 11, 9, 15, 13)

X <- cbind(x0, x1, x2)
X

### Matrizes de interesse

k <- 2

n <- length(y)

jn <- matrix(data = 1, nrow = n, ncol = 1)

Jnn <- jn %*% t(jn)

In <- diag(n)

X12 <- cbind(x1, x2)
X12

### Matriz X12 centrada

Xc <- (In - (1 / n) * Jnn) %*% X12
Xc

## Teste de Regressão Global (beta1 = 0)

### Total

#### Soma de quadrados
SQTot <- t(y) %*% (In - (1 / n) * Jnn) %*% y
SQTot

#### Graus de Liberdade
gl_tot <- n - 1
gl_tot

### Regressão

#### Beta1
Beta1 <- solve(t(Xc) %*% Xc) %*% t(Xc) %*% y
Beta1

#### Soma de quadrados
SQReg <- t(Beta1) %*% t(Xc) %*% y
SQReg

#### Graus de Liberdade
gl_reg <- k
gl_reg

#### Quadrado Médio
QMReg <- SQReg / gl_reg
QMReg

### Resíduos

#### Soma de quadrados
SQRes <- SQTot - SQReg
SQRes

#### Graus de Liberdade
gl_res <- n - k - 1
gl_res

#### Quadrado Médio
QMRes <- SQRes / gl_res
QMRes

### Estatística F e p-valor
Fcalc1 <- QMReg / QMRes
Fcalc1

ftab1 <- qf(0.95, k, n - k - 1)
ftab1

p_valor1 <- 1 - pf(Fcalc1, k, n - k - 1)
p_valor1

## Teste sobre um subconjunto dos betas (beta2 = 0)

X01 <- cbind(x0, x1)  # Variáveis importantes em X1Beta1
X01

X2 <- as.matrix(x2)   # Variáveis desprezíveis em X2Beta2
X2

A1 <- X %*% solve(t(X) %*% X) %*% t(X) 
A1

A2 <- X01 %*% solve(t(X01) %*% X01) %*% t(X01)
A2

SQB2_B1 <- t(y) %*% (A1 - A2) %*% y
SQB2_B1

h <- ncol(X2)
h

QMB2_B1 <- SQB2_B1 / h
QMB2_B1

Fcalc2 <- QMB2_B1 / QMRes
Fcalc2

ftab2 <- qf(0.95, h, n - k - 1)
ftab2

p_valor2 <- 1 - pf(Fcalc2, h, n - k - 1)
p_valor2

## Hipótese Linear Geral (H0: C*Beta = 0 ou H0: beta1=beta2=0)
### H0: B1 = B2 = 0 usando HIPÓTESE LINEAR GERAL

C <- matrix(data = c(0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE)
C

Beta <- solve(t(X) %*% X) %*% t(X) %*% y
Beta

CBeta <- C %*% Beta
CBeta

SQHip <- t(CBeta) %*% solve(C %*% solve(t(X) %*% X) %*% t(C)) %*% CBeta
SQHip

gl_hip <- nrow(C)
gl_hip

QMHip <- SQHip / gl_hip
QMHip

Fcalc3 <- QMHip / QMRes
Fcalc3

ftab3 <- qf(0.95, gl_hip, n - k - 1)
ftab3

p_valor3 <- 1 - pf(Fcalc3, gl_hip, n - k - 1)
p_valor3



# 7.2 Regressão Linear Múltipla - Teste Hip -------------------------------

## Dados

y1 <- c(41.5,33.8,27.7,21.7,19.9,15.0,12.2,4.3,19.3,6.4,37.6,18.0,26.3,9.9,25.0,14.1,15.2,15.9,19.6)

y2 <- c(45.9,53.3,57.5,58.8,60.6,58.0,58.6,52.4,56.9,55.4,46.9,57.3,55.0,58.9,50.3,61.1,62.9,60.0,60.6)

n <- length(y2)
x0 <- rep(1, n)

x1 <- c(162,162,162,162,172,172,172,172,167,177,157,167,167,167,167,177,177,160,160)

x2 <- c(23,23,30,30,25,25,30,30,27.5,27.5,27.5,32.5,22.5,27.5,27.5,20,20,34,34)

x3 <- c(3,8,5,8,5,8,5,8,6.5,6.5,6.5,6.5,6.5,9.5,3.5,6.5,6.5,7.5,7.5)

x11 <- x1*x1

x22 <- x2*x2

x33 <- x3*x3

x12 <- x1*x2

x13 <- x1*x3

x23 <- x2*x3

X <- cbind(x0, x1, x2, x3, x11, x22, x33, x12, x13, x23)
X

Jnn <- matrix(data = 1, nrow = n, ncol = n)

In <- diag(n)

## Modelo completo

Beta <- solve(t(X) %*% X) %*% t(X) %*% y2
Beta |> round(4)

SQTotal <- t(y2) %*% (In - (1 / n) * Jnn) %*% y2
SQTotal

gltotal <- n - 1
gltotal

SQReg <- t(y2) %*% (X %*% solve(t(X) %*% X) %*% t(X) - (1 / n) * Jnn) %*% y2
SQReg

k <- ncol(X) - 1
gl_reg <- k
gl_reg

SQRes <- SQTotal - SQReg
SQRes

gl_res <- n - k - 1
gl_res

QMRes <- SQRes / gl_res
QMRes

## Hipótese 1 - H0: Beta4 = Beta5 = Beta6 = Beta7 = Beta8 = Beta9 = 0
X1 <- cbind(x0, x1, x2, x3)
X1

SQReg1 <- t(y2) %*% (X1 %*% solve(t(X1) %*% X1) %*% t(X1) - (1 / n) * Jnn) %*% y2
SQReg1

gl_regB1 <- ncol(X1) - 1
gl_regB1

SQB2_B1 <- SQReg - SQReg1
SQB2_B1

gl_regB2 <- gl_reg - gl_regB1
h <- gl_regB2
h

QMB2_B1 <- SQB2_B1 / h
QMB2_B1

Fcalc1 <- QMB2_B1 / QMRes
Fcalc1

p_valor1 <- 1 - pf(Fcalc1, h, gl_res)
p_valor1


## Hipótese 2 - H0: Beta1 = Beta2 = Beta3 = 0
X1 <- cbind(x0, x11, x22, x33, x12, x13, x23)
X1

SQReg1 <- t(y2) %*% (X1 %*% solve(t(X1) %*% X1) %*% t(X1) - (1 / n) * Jnn) %*% y2
SQReg1

gl_regB1 <- ncol(X1) - 1
gl_regB1

SQB2_B1 <- SQReg - SQReg1
SQB2_B1

gl_regB2 <- gl_reg - gl_regB1
h <- gl_regB2
h

QMB2_B1 <- SQB2_B1 / h
QMB2_B1

Fcalc2 <- QMB2_B1 / QMRes
Fcalc2

p_valor2 <- 1 - pf(Fcalc2, h, gl_res)
p_valor2


## Hipótese Linear Geral (Global)
### H0: Beta1 = Beta2 = Beta3 = 0
C <- matrix(
  c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0), 
  ncol = 10, nrow = 3, byrow = TRUE
)
C

CBeta <- C %*% Beta
CBeta

SQHip <- t(CBeta) %*% solve(C %*% solve(t(X) %*% X) %*% t(C)) %*% CBeta
SQHip

q <- nrow(C)

QMHip <- SQHip / q
QMHip

Fcalc3 <- QMHip / QMRes
Fcalc3

p_valor3 <- 1 - pf(Fcalc3, q, gl_res)
p_valor3

### H0: Beta1 = Beta2 = Beta3
C1 <- matrix(
  c(0, 1, -1,  0, 0, 0, 0, 0, 0, 0,
    0, 0,  1, -1, 0, 0, 0, 0, 0, 0), 
  ncol = 10, nrow = 2, byrow = TRUE
)
C1

C2 <- matrix(
  c(0, 2, -1, -1, 0, 0, 0, 0, 0, 0,
    0, 0,  1, -1, 0, 0, 0, 0, 0, 0), 
  ncol = 10, nrow = 2, byrow = TRUE
)
C2

C3 <- matrix(
  c(0, 1, 0,  -1, 0, 0, 0, 0, 0, 0,
    0, 0, 1, -1, 0, 0, 0, 0, 0, 0), 
  ncol = 10, nrow = 2, byrow = TRUE
)
C3

#### Função
hlg <- function(C) {
  
  # Calculando CBeta
  CBeta <- C %*% Beta |> round(3)
  
  # Calculando SQHip
  SQHip <- t(CBeta) %*% solve(C %*% solve(t(X) %*% X) %*% t(C)) %*% CBeta
  
  # Obtendo o número de linhas de C
  q <- nrow(C)
  
  # Calculando QMHip
  QMHip <- SQHip / q
  
  # Calculando Fcalc
  Fcalc3 <- QMHip / QMRes
  
  # Calculando p_valor
  p_valor3 <- 1 - pf(Fcalc3, q, gl_res)
  
  # Retornando resultados
  return(list(CBeta = CBeta, SQHip = SQHip, QMHip = QMHip, Fcalc = Fcalc3, p_valor = p_valor3))
}

#### Teste
# Definindo os valores para C
C_lista <- list(C1 = C1, C2 = C2, C3 = C3)

# Criando uma lista para armazenar os resultados
resultados <- list()

# Calculando os resultados para cada objeto C
for (i in 1:length(C_lista)) {
  resultados[[paste0("C", i)]] <- hlg(C_lista[[i]])
}

#### Resultado
tabela_resultados <- data.frame(
  SQHip = sapply(resultados, function(x) x$SQHip),
  QMHip = sapply(resultados, function(x) x$QMHip),
  Fcalc = sapply(resultados, function(x) x$Fcalc),
  p_valor = sapply(resultados, function(x) x$p_valor)
)

tabela_resultados



# 7.3 Regressão Linear Múltipla - Teste Hip -------------------------------

## Dados
y1 <- c(41.5,33.8,27.7,21.7,19.9,15.0,12.2,4.3,19.3,6.4,37.6,18.0,26.3,9.9,25.0,14.1,15.2,15.9,19.6)

n <- length(y1)
x0 <- rep(1, n)

x1 <- c(162,162,162,162,172,172,172,172,167,177,157,167,167,167,167,177,177,160,160)

x2 <- c(23,23,30,30,25,25,30,30,27.5,27.5,27.5,32.5,22.5,27.5,27.5,20,20,34,34)

x3 <- c(3,8,5,8,5,8,5,8,6.5,6.5,6.5,6.5,6.5,9.5,3.5,6.5,6.5,7.5,7.5)

Jnn <- matrix(data = 1, nrow = n, ncol = n)

In <- diag(n)

X <- cbind(x0, x1, x2, x3)
X

k <- ncol(X) - 1	     # número de variáveis regressoras
k

Beta <- solve(t(X) %*% X) %*% t(X) %*% y1
Beta |> round(4)

## Hipótese 1 - Teste de Regressão Global - H0: 2Beta1 = 2Beta2 = 2Beta3 = 0

C <- matrix(
  c(0, -2,  2,  0, 
    0,  0,  2, -1),
  nrow = 2, byrow = TRUE
)
C

CBeta <- C %*% Beta
CBeta

SQHip <- t(CBeta) %*% solve(C %*% (solve(t(X) %*% X)) %*% t(C)) %*% CBeta
SQHip

gl_Hip <- nrow(C)
gl_Hip

QMHip <- SQHip / gl_Hip
QMHip

SQRes <- t(y1) %*% (In - X %*% solve(t(X) %*% X) %*% t(X)) %*% y1
SQRes

gl_res <- n - k - 1
gl_res

QMRes <- SQRes / gl_res
QMRes

Fcalc1 <- QMHip / QMRes
Fcalc1

Ftab1 <- qf(0.95, gl_Hip, gl_res)
Ftab1

p_valor1 <- 1 - pf(Fcalc1, gl_Hip, gl_res) 
p_valor1

## Hipótese 2 - Teste de Regressão Hipótese Linear Geral para diversas matrizes C - H0: Beta1 = Beta2 = Beta3 = 0

### Matriz C1

C1 <- matrix(
  c(0, 1, -1,  0, 
    0, 0,  1, -1),
  nrow = 2, byrow = TRUE
)
C1

SQ_C1Beta <- t(C1 %*% Beta) %*% solve(C1 %*% (solve(t(X) %*% X)) %*% t(C1)) %*% C1 %*% Beta
SQ_C1Beta

### Matriz C2

C2 <- matrix(
  c(0, 1, -1,  0, 
    0, 1,  0, -1),
  nrow = 2, byrow = TRUE
)
C2

SQ_C2Beta <- t(C2 %*% Beta) %*% solve(C2 %*% (solve(t(X) %*% X)) %*% t(C2)) %*% C2 %*% Beta
SQ_C2Beta

### Matriz C3

C3 <- matrix(
  c(0, 2, -1, -1, 
    0, 0,  1, -1),
  nrow = 2, byrow = TRUE
)
C3

SQ_C3Beta <- t(C3 %*% Beta) %*% solve(C3 %*% (solve(t(X) %*% X)) %*% t(C3)) %*% C3 %*% Beta
SQ_C3Beta

## Hipótese 3 - H0: Beta1 = Beta2

### Hipótese Linear Geral

C <- matrix(c(0, 1, -1, 0), nrow = 1, byrow = TRUE)
C

CBeta <- C %*% Beta
CBeta

gl_Hip <- nrow(C)
gl_Hip

SQHip <- t(CBeta) %*% solve(C %*% (solve(t(X) %*% X)) %*% t(C)) %*% CBeta
SQHip

QMHip <- SQHip / gl_Hip
QMHip

Fcalc2 <- QMHip / QMRes
Fcalc2

p_valor2 <- 1 - pf(Fcalc2, gl_Hip, gl_res)
p_valor2


### Modelo completo e modelo reduzido

x12 <- x1 + x2           # Modelo reduzido: y = B0 + B12(x1+x2) + B3x3 + e
x12

Xr <- cbind(x0, x12, x3) # Matriz Xr do modelo reduzido
Xr

SQHip2 <- t(y1) %*% (X %*% solve(t(X) %*% X) %*% t(X) - Xr %*% solve(t(Xr) %*% Xr) %*% t(Xr)) %*% y1
SQHip2



# 7.4 Regressão Linear Múltipla - Métodos Bonferroni e Scheffé ------------

y <- c(2,3,2,7,6,8,10,7,8,12,11,14)

n <- length(y)

x0 <- rep(1, n)

x1 <- c(0,2,2,2,4,4,4,6,6,6,8,8)

x2 <- c(2,6,7,5,9,8,7,10,11,9,15,13)

X <- cbind(x0, x1, x2)
X

k <- ncol(X) - 1

In <- diag(n)

Beta <- solve(t(X) %*% X) %*% t(X) %*% y
Beta

a0 <- matrix(data = c(1, 0, 0), ncol = 1, byrow = TRUE)
a0

Beta0 <- t(a0) %*% Beta
Beta0

a1 <- matrix(data = c(0, 1, 0), ncol = 1, byrow = TRUE)
a1

Beta1 <- t(a1) %*% Beta
Beta1

a2 <- matrix(data = c(0,0,1), ncol = 1, byrow = TRUE)
a2

Beta2 <- t(a2) %*% Beta
Beta2

SQRes <- t(y) %*% (In - X %*% solve(t(X) %*% X) %*% t(X)) %*% y
SQRes

gl_res <- n - k - 1
gl_res

QMRes <- SQRes / gl_res
QMRes

s <- sqrt(QMRes)
s

G <- solve(t(X) %*% X)
G

t1 <- Beta1 / (s %*% sqrt(t(a1) %*% G %*% a1))
t1

t2 <- Beta2 / (s %*% sqrt(t(a2) %*% G %*% a2))
t2

p_valor_t1 <- 2 * (1 - pt(abs(t1), gl_res))
p_valor_t1

p_valor_t2 <- 2 * (1 - pt(abs(t2), gl_res))
p_valor_t2

p <- 1 - 0.05 / (2 * k)
p

t_tab <- qt(0.975, n - k - 1)
t_tab

t_Bon <- qt(p, n - k - 1) 		# calcula t-tabelado para Método de Bonferroni
t_Bon

t_Scheffe <- sqrt((k + 1) %*% qf(0.95, k + 1, n - k - 1)) |> as.numeric() # calcula t-tabelado para Método de Scheffé
t_Scheffe




# 8. Modelos de Análise de Variância (ANOVA) ------------------------------

## 8.1 ANOVA balanceada com um fator - Estimação --------------------------

y <- as.vector(c(14,16,15,18,19,17))
y

X <- matrix(c(
  1, 1, 0, 
  1, 1, 0,
  1, 1, 0, 
  1, 0, 1, 
  1, 0, 1, 
  1, 0, 1
),
ncol = 3, byrow = TRUE
)
X

XLX <- t(X) %*% X
XLX

XLy <- t(X) %*% y
XLy

install.packages("MASS")
library(MASS)     # função ginv()

rank_X <- sum(diag(X %*% ginv(X)))
rank_X

rank_XLX <- sum(diag(XLX %*% ginv(XLX)))
rank_XLX

p <- ncol(X)         # Número de parâmetros
p

k <- rank_X          # Posto(X)
k


defRank <- p - k     # Deficiência de rank - mostra que é de posto incompleto
defRank

### Inversa generalizada

Beta0 <- ginv(XLX) %*% t(X) %*% y
Beta0

### Solução com restrição $\tau_1 = 0$ (R)

y_hat0 <- X %*% Beta0
y_hat0

TR <- c(0, 1, 0)
TR

W <- rbind(X, TR)
W

z <- matrix(c(y, 0), ncol = 1)
z

Beta1 <- round(solve(t(W) %*% W) %*% t(W) %*% z, 2)
Beta1

y_hat1 <- X %*% Beta1
y_hat1

### Solução com restrição $\tau_2 = 0$ (SAS)

TS <- c(0, 0, 1)
TS

W <- rbind(X, TS)
W

z <- matrix(c(y, 0), ncol = 1)
z

Beta2 <- round(solve(t(W) %*% W) %*% t(W) %*% z, 2)
Beta2

y_hat2 <- X %*% Beta2
y_hat2

### Solução com restrição $\tau_1 + \tau_2 = 0$ (estatística experimental)

TE <- c(0, 1, 1)
TE

W <- rbind(X, TE)
W

z <- matrix(c(y, 0), ncol = 1)
z

Beta3 <- round(solve(t(W) %*% W) %*% t(W) %*% z, 2)
Beta3

y_hat3 <- X %*% Beta3
y_hat3

### Funções estimáveis

#### $\beta = \mu + \tau_1$

L1 <- c(1, 1, 0)
L1

L1Beta0 <- L1 %*% Beta0
L1Beta0

L1Beta1 <- L1 %*% Beta1
L1Beta1

L1Beta2 <- L1 %*% Beta2
L1Beta2

L1Beta3 <- L1 %*% Beta3
L1Beta3

#### $\beta = \tau_1 + \tau_2$

L2 <- c(0, 1, 1)
L2

L2Beta0 <- L2 %*% Beta0
L2Beta0

L2Beta1 <- L2 %*% Beta1
L2Beta1

L2Beta2 <- L2 %*% Beta2
L2Beta2

L2Beta3 <- L2 %*% Beta3
L2Beta3


## 8.2 ANOVA balanceada com um fator - Teste de Hipótese ------------------

# install.packages("MASS")
library(MASS)     # função ginv()

y <- as.vector(
  c(14.29, 19.10, 19.09, 16.25, 15.09, 16.61, 19.63,
    20.06, 20.64, 18.00, 19.56, 19.47, 19.07, 18.38,
    20.04, 26.23, 22.74, 24.04, 23.37, 25.02, 23.27)
)
y

X <- matrix(c(
  rep(1,21),
  rep(1,7),rep(0,14),
  rep(0,7),rep(1,7),rep(0,7),
  rep(0,14),rep(1,7)
),
ncol = 4, byrow = FALSE)
X

kn <- length(y)
kn

p <- ncol(X)
p

k <- sum(diag(X %*% ginv(X)))
k

XLX <- t(X) %*% X
XLX

XLy <- t(X) %*% y
XLy

### Vetor com as estimativas a partir da inversa generalizada

Beta <- ginv(XLX) %*% XLy
Beta

### Vetor com as estimativas a partir da inversa generalizada utilizando a submatriz diagonal

igXLX <- (1/7) * matrix(
  c(0, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1),
  ncol = 4, byrow = TRUE
)
igXLX

Beta2 <- igXLX %*% XLy
Beta2

## Funções estimáveis

### $\beta = \alpha_1 - \alpha_2$ (estimável)

L1 <- c(0, 1, -1, 0)
L1

L1Beta <- t(L1) %*% Beta
L1Beta

L1Beta2 <- t(L1) %*% Beta2
L1Beta2

### $\beta = \alpha_1 + \alpha_2 + \alpha_3$ (não estimável)

L2 <- c(0, 1, 1, 1)
L2

L2Beta <- t(L2) %*% Beta
L2Beta

L2Beta2 <- t(L2) %*% Beta2
L2Beta2

## Teste de hipótese: $H_0: \mu_1 = \mu_2 = \mu_3$

### Modelo Completo x Modelo Reduzido

In <- diag(kn)
Jn <- matrix(1, nrow = kn, ncol = kn)

# Total
Tot <- In - (1 / kn) * Jn

SQTotal <- t(y) %*% Tot %*% y
SQTotal

gl_total <- round(sum(diag(Tot %*% ginv(Tot))))
gl_total

# Tratamentos
A <- X %*% ginv(t(X) %*% X) %*% t(X) - (1 / kn) * Jn

SQTrat <- t(y) %*% A %*% y
SQTrat

gl_trat <- round(sum(diag(A %*% ginv(A))))
gl_trat

QMTrat <- SQTrat / gl_trat
QMTrat

# Resíduo
B <- In - X %*% ginv(t(X) %*% X) %*% t(X)

SQRes <- t(y) %*% B %*% y
SQRes

gl_res <- round(sum(diag(B %*% ginv(B))))
gl_res

QMRes <- SQRes / gl_res
QMRes

Fcalc <- QMTrat / QMRes
Fcalc

ftab <- qf(0.95, gl_trat, gl_res)
ftab

p_valor <- 1 - pf(Fcalc, gl_trat, gl_res)
p_valor

### Contrastes

#### Caso 1: As linhas são linearmente independentes (l.i.) e ortogonais

C1 <- c(0, 2, -1, -1)
C1Beta <- C1 %*% Beta
C1Beta

C2 <- c(0, 0, 1, -1)
C2Beta <- C2 %*% Beta
C2Beta

# Partição da SQTrat
SQC1Beta <- t(C1Beta) %*% solve(t(C1) %*% ginv(t(X) %*% X) %*% C1) %*% (C1Beta)
SQC1Beta

SQC2Beta <- t(C2Beta) %*% solve(t(C2) %*% ginv(t(X) %*% X) %*% C2) %*% (C2Beta)
SQC2Beta

SomaSQ <- SQC1Beta + SQC2Beta
SomaSQ

SQTrat

# QM
QMC1Beta <- SQC1Beta / 1
QMC1Beta

QMC2Beta <- SQC2Beta / 1
QMC2Beta

# F e p-valor
FC1 <- QMC1Beta / QMRes
FC1

FC2 <- QMC2Beta / QMRes
FC2

ftab <- qf(0.95, 1, 18)
ftab

p_valorC1 <- 1 - pf(FC1, 1, gl_res)
p_valorC1

p_valorC2 <- 1 - pf(FC2, 1, gl_res)
p_valorC2

#### Caso 2: As linhas são l.i. e não ortogonais

C1 <- c(0, 2, -1, -1)
C1Beta <- C1 %*% Beta
C1Beta

C3 <- c(0, 1, 0, -1)
C3Beta <- C3 %*% Beta
C3Beta

# Partição da SQTrat
SQC1Beta <- t(C1Beta) %*% solve(t(C1) %*% ginv(t(X) %*% X) %*% C1) %*% (C1Beta)
SQC1Beta

SQC3Beta <- t(C3Beta) %*% solve(t(C3)%*% ginv(t(X) %*% X) %*% C3) %*% (C3Beta)
SQC3Beta

# QM
QMC1Beta <- SQC1Beta / 1
QMC1Beta

QMC3Beta <- SQC3Beta / 1
QMC3Beta

# F e p-valor
FC1 <- QMC1Beta / QMRes
FC1

FC3 <- QMC3Beta / QMRes
FC3

ftab <- qf(0.95, 1, 18)
ftab

p_valorC1 <- 1 - pf(FC1, 1, gl_res)
p_valorC1

p_valorC3 <- 1 - pf(FC3, 1, gl_res)
p_valorC3

#### A soma das somas de quadrados dos contrastes não ortogonais não resulta na soma de quadrados dos tratamentos:

SomaSQ <- SQC1Beta + SQC3Beta
SomaSQ

SQTrat



## 8.3 ANOVA balanceada com dois fatores ----------------------------------

# install.packages("MASS")
library(MASS)

y <- c(39.02, 38.79, 35.74, 35.41, 37.02, 36.00, 38.96, 39.01, 35.58, 35.52, 35.70, 36.04)
y

X <- matrix(c(
  1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
  1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
  1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
  1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
  1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
  1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
  1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
  1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
  1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
  1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
  1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1
), 
nrow = 12, byrow = TRUE)
X

rank_X <- sum(diag(ginv(X) %*% X))
rank_X

npar <- ncol(X)    # número de parâmetros
npar

a <- 2    # níveis do fator A
b <- 3    # níveis do fator B
n <- 2    # número de repeticões

abn <- a * b * n    # número total de observacões
abn

X0 <- X[, 1]        # constante
X0

XA <- X[, 2:3]      # fator A
XA

XB <- X[, 4:6]      # fator B
XB

XAB <- X[, 7:12]    # combinacão dos níveis dos dois fatores
XAB

### Estimação dos parâmetros

#### Estimação de beta usando inversa generalizada de Moore-Penrose

Beta <- ginv(t(X) %*% X) %*% t(X) %*% y
Beta

#### Estimação de beta usando inversa generalizada (Searle) de X'X

XLX <- t(X) %*% X
XLX

iXLX <- (1 / 2) * kronecker(
  matrix(c(0,0,0,1), byrow = TRUE, ncol = 2), diag(6)
)
fractions(iXLX)

Betag <- iXLX %*% t(X) %*% y
Betag

### Teste de hipótese

In <- diag(abn)
Jn <- matrix(1, nrow = abn, ncol = abn)

# Total
Tot <- In - (1 / abn) * Jn

SQTotal <- t(y) %*% Tot %*% y
SQTotal

gl_total <- round(sum(diag(Tot %*% ginv(Tot))))
gl_total

# Resíduo
PR <- diag(abn) - X %*% ginv(t(X) %*% X) %*% t(X)

SQRes <- t(y) %*% PR %*% y
SQRes

gl_res <- round(sum(diag(PR %*% ginv(PR))))
gl_res

QMRes <- SQRes / gl_res
QMRes

# Cálculo SQAxB (interação) - forma quadrática
X1 <- cbind(X0, XA, XB)
X1

PAB <- X %*% ginv(t(X) %*% X) %*% t(X) - X1 %*% ginv(t(X1) %*% X1) %*% t(X1)
PAB

SQAB <- t(y) %*% PAB %*% y
SQAB

glAB <- round(sum(diag(ginv(PAB) %*% PAB)))
glAB

QMAB <- SQAB / glAB
QMAB

FAB <- QMAB / QMRes
FAB

ftabAB <- qf(0.95, glAB, gl_res)
ftabAB

p_valorAB <- 1 - pf(FAB, glAB, gl_res)
p_valorAB

# fator A - SQ(A)
PA <- XA %*% ginv(t(XA) %*% XA) %*% t(XA) - Jn/abn
PA

SQA <- t(y) %*% PA %*% y
SQA

glA <- round(sum(diag(ginv(PA) %*% PA)))
glA

QMA <- SQA / glA
QMA

FA <- QMA / QMRes
FA

ftabA <- qf(0.95, glA, gl_res)
ftabA

p_valorA <- 1 - pf(FA, glA, gl_res)
p_valorA

# fator B - SQ(A)
PB <- XB %*% ginv(t(XB) %*% XB) %*% t(XB) - Jn/abn
PB

SQB <- t(y) %*% PB %*% y
SQB

glB <- round(sum(diag(ginv(PB) %*% PB)))
glB

QMB <- SQB / glB
QMB

FB <- QMB / QMRes
FB

ftabB <- qf(0.95, glB, gl_res)
ftabB

p_valorB <- 1 - pf(FB, glB, gl_res)
p_valorB

### Somas de quadrados usando Hipótese Linear Geral

CA <- 
  (1/3) * matrix(
    c(0, 3, -3, 0, 0, 0, 1, 1, 1, -1, -1, -1), nrow = 1, byrow = TRUE
  )
CA

CB <- 
  (1/2) * matrix(
    c(0, 0, 0, 2, -2, 0, 1, -1, 0, 1, -1, 0,
      0, 0, 0, 2, 0, -2, 1, 0, -1, 1, 0, -1), nrow = 2, byrow = TRUE
  )
CB

CAxB <- matrix(
  c(0, 0, 0, 0, 0, 0, 1, -1, 0, -1, 1, 0,
    0, 0, 0, 0, 0, 0, 1, 0, -1, -1, 0, 1), nrow = 2, byrow = TRUE
)
CAxB

SQ_CA <- t(CA %*% Beta) %*% solve(CA %*% ginv(t(X) %*% X) %*% t(CA)) %*% (CA %*% Beta)
SQ_CA

SQ_CB <- t(CB %*% Beta) %*% solve(CB %*% ginv(t(X) %*% X) %*% t(CB)) %*% (CB %*% Beta)
SQ_CB

SQ_CAxB <- t(CAxB %*% Beta) %*% solve(CAxB %*% ginv(t(X) %*% X) %*% t(CAxB)) %*% (CAxB %*% Beta)
SQ_CAxB

### Estimabilidade no Modelo Superparametrizado Sem Restrições

#Verifica se LBeta = a1-a2 e LBeta = a1-a2 +(1/3(g11+g12+g13)-(1/3)(g21=g22+g23) são estimáveis - não é estimável

L1 <- t(matrix(c(0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1))
L1

ver <- t(X) %*% X %*% ginv(t(X) %*% X)

verL1 <- ver %*% L1 |> round(2)
verL1

L1Beta <- t(L1) %*% Beta
L1Beta

#Estimabilidade de L2*Beta = a1-a2 + (1/3(g11+g12+g13)-(1/3)(g21+g22+g23) é estimavel no modelo SEM restricao nos parametros - é estimavel!

L2 <- (1/3) * t(matrix(c(0, 3, -3, 0, 0, 0, 1, 1, 1, -1, -1, -1), nrow = 1))
L2

ver <- t(X) %*% X %*% ginv(t(X) %*% X)

verL2 <- ver %*% L2 |> round(2)
verL2

L2Beta <- t(L2) %*% Beta
L2Beta

# Matriz T de condições marginais: T*Beta = 0

t <- matrix(
  c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1), 
  nrow = 7, ncol = 12, byrow = TRUE)

rank_T <- sum(diag(ginv(t) %*% t))
rank_T

W <- rbind(X, t)

rank_W <- sum(diag(ginv(W) %*% W))
rank_W

yr <- c(y, rep(0, 7))
yr

Beta_R <- solve(t(W) %*% W) %*% t(W) %*% yr   # Beta sujeito às condições marginais
Beta_R

### Estimabilidade no Modelo Superparametrizado Com Restrições

# LiBeta é estimável se Li = verLi. Verifica que L1Beta = a1-a2 é ESTIMÁVEL no modelo COM RESTRIÇÃO nos parâmetros

L1 <- t(matrix(c(0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1))
L1

ver <- t(W) %*% W %*% solve(t(W) %*% W)

verL1 <- ver %*% L1 |> round()
verL1

L1Beta_r <- t(L1) %*% Beta_R
L1Beta_r # Mostra que L1Beta = a1-a2 É estimável no modelo', 'COM restrição nos parâmetros

# Verifica que L2Beta = a1-a2 a1-a2 + (1/3(g11+g12+g13) - (1/3)(g21+g22+g23) é ESTIMÁVEL no modelo', 'COM RESTRIÇÃO nos parâmetros

L22 <- (1/3) * t(matrix(c(0, 3, -3, 0, 0, 0, 1, 1, 1, -1, -1, -1), nrow = 1))
L22 |> round(3)

verL22 <- ver %*% L22 |> round(3)
verL22

L2Beta_R <- t(L22) %*% Beta_R
L2Beta_R      # Mostra que L2Beta = a1-a2 +(1/3(g11+g12+g13)-(1/3)(g21=g22+g23) É estimável no modelo', 'COM restrição nos parâmetros

### Hipótese linear geral

# ESTIMABILIDADE NO MODELO SUPERPARAMETRIZADO SEM RESTRIÇÃO
CA <- (1 / 3) * c(0,  3, -3,  0,  0,  0,  1,  1,  1, -1, -1, -1)
CA |> round(2)

CABeta <- CA %*% Beta_R
CABeta

# No modelo COM restrição nos parâmetros:

CA

CABeta_R <- CA %*% Beta_R
CABeta_R

SQ_A <- t(CABeta) %*% solve(t(CA) %*% ginv(t(X) %*% X) %*% CA) %*% CABeta
SQ_A

SQ_A_R = t(CABeta_R) %*% solve(t(CA) %*% solve(t(W) %*% W) %*% CA) %*% CABeta_R
SQ_A_R

