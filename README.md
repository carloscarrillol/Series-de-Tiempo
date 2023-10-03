# Series-de-Tiempo
Puedes ver el código ejecutado en
#### https://rpubs.com/CarlosCarrillol/1091190



---
title: "Series de Tiempo"
output: html_notebook
---




<p style="text-align:center;">
*Quiero que quede claro que no tengo idea de lo que estoy apunto de escribir, únicamente es una compilación de documentos que he encontrado sobre modelado de series de tiempo. Por lo tanto, descubriremos juntos y al mismo tiempo (eso no es verdad) que hace cada una de las siguientes lineas. Por favor no juzgues, tuve la decencia de advertirte.*
</p>

<p style="text-align:center;">
-C
</p>

<br><br><br><br><br><br>

Vamos a utilizar la paqueteria zoo y las librerias forecast, para hacer pronosticos (supongo) y tseries.


```{r}
install.packages("zoo")
library(forecast)
library(tseries)
```

# Procesos Autorregresivos $AR(p)$

$\textbf{Definición 1.1}$ Un $\textit{proceso autorregresivo}$ es un proceso en el que los datos futuros dependen linealmente de los datos historicos. El modelo se define como

$$Y_{t}=C+\sum_{i=1}^{p}\phi_{i}Y_{t-i}+\varepsilon_{t} \tag{1.1}$$
De manera matricial podemos reescribirlo como
$$Y_{t}= Y_{t-i}\cdot\Phi+\varepsilon_{t}$$
donde $c$ es constante y $\Phi$ es un vector columna de $p$ entradas,  $Y_{t-i}$ es un vector de $p$ entradas donde $i=1,...,p$, tal que cada entrada es la i-ésima observación antes de $t$, $\varepsilon_{t}$ es un término de error que se distribuyue como Ruido Blanco (RB de ahora en adelante) con media 0 y varianza $\sigma_{\varepsilon}^{2}$. De manera extensiva:

$$Y_{t}=C+\begin{pmatrix}Y_{t-1} & Y_{t-2} & \cdots & Y_{t-p}\end{pmatrix} \begin{pmatrix} \phi_{1}\\ \phi_{2}\\ \vdots \\ \phi_{p}\end{pmatrix} + \varepsilon_{t} \tag{2.1}$$
Note que el modelo AR es un modelo de regresión lineal de la serie $Y_{t}$ y sus valores rezagados.

Para la definición no es necesario suponer que el ruido blanco tiene una distribución particular, pero más adelante supondremos que el ruido blanco tiene distribución normal.
Una notación muy útil para considerar este tipo de modelos se basa en el uso del operador de retardo $B$. El modelo $AR(p)$ es


$$(1-\phi_{1} B - \phi_{2} B^{2}-\phi_{3}B^{3}-\cdots -\phi_{p}B^{p})Y_{t}=\varepsilon_{t} \tag{3.1}$$
La ecuación (3.1) se puede resumir usando la siguiente notación:
$$\phi(B)X_{t}=\varepsilon_{t}$$
Donde $\phi(B)$ es un polinomio de grado $p$ llamado operador autorregresivo, definido por

$$\phi(B)=(1-\phi_{1} B - \phi_{2} B^{2}-\phi_{3}B^{3}-\cdots -\phi_{p}B^{p})$$

### Proceso Autoregresivo de Orden 1 $AR(1)$

En este proceso la variable $Y_{t}$ depende del valor rezagado $Y_{t-1}$
$$Y_{t}=\phi Y_{t-1}+\varepsilon_{t} \tag{4.1}$$
$\textbf{Definición 2.1}$ Se dice que un proceso estocástico $Y_{t}$ es $\textit{estríctamente estacionario}$ si todas las distribuciones conjuntas de las observaciones en diferentes puntos temporales son las mismas.


$\textbf{Definición 3.1}$ Se dice que un proceso estocástico $Y_{t}$ es $\textit{débilmente estacionario}$ si su media $\mathbb{E}(Y_{t})$ y su varianza $\mathbb{V}(Y_{t})$ son constantes a través del tiempo y la autocovarianza $\mathbb{C}(Y_{t},Y_{t+h})$ entre cualesquiera dos observaciones no depende de $t$.


$\textbf{Teorema  1}$ Si un proceso estocastico $Y_{t}$ es estríctamente estacionario, entonces es débilmente estacionario.

Si iteramos la ecuación (4.1) hacia atrás, tenemos que,

\begin{align}
Y_{t}&=\left(\phi(\phi Y_{t-2}+\varepsilon_{t-1})\right)+\varepsilon_{t}\\ \\
    &= \phi^{2}Y_{t-2}+\phi\varepsilon_{t-1}+\varepsilon_{t}\\ \\
\Leftrightarrow Y_{t}&=\phi^{k}Y_{t-k}+\sum_{i=0}^{k-1}\phi^{i}\varepsilon_{t-i}
\end{align}

Si continuamos este proceso de iteración al pasado (siempre que $|\phi|<1)$ y $Y_{t}$ sea estacionario, podemos representar el modelo $AR(1)$ como un proceso lineal:

$$Y_{t}=\lim_{k\rightarrow\infty}\phi^{k}Y_{t-k}+\sum_{i=0}^{k-1}\phi^{i}\varepsilon_{t-i}$$
Aplicando el límite y usando el hecho de que $|\phi|<1$, tenemos que,
$$Y_{t}=\sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t-i}\tag{5.1}\sim MA(\infty)$$
Por lo tanto, la ecuación (5.1) define un proceso $AR(1)$ estacionario con media,
\begin{align}
\mathbb{E}(Y_{t})&=\mathbb{E}\left(\sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t-i}\right)\\ \\
          &=\sum_{i=0}^{\infty}\phi^{i}\mathbb{E}(\varepsilon_{t-i})=0
\end{align}

varianza,

\begin{align}
\mathbb{V}(Y_{t})&=\mathbb{V}\left(\sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t-i}\right)\\ \\
        &=\sum_{i=0}^{\infty}\phi^{2i}\mathbb{V}(\varepsilon_{t-i})=\sigma_{\varepsilon}^{2}\sum_{i=0}^{\infty}\phi^{2i}
\end{align}

y autocovarianza,

\begin{align}
\mathbb{C}(Y_{t},Y_{t+h})&=\mathbb{E}\left(\sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t-i}\cdot \sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t+h-i}\right) - \mathbb{E}\left(\sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t-i}\right)\cdot\mathbb{E}\left(\sum_{i=0}^{\infty}\phi^{i}\varepsilon_{t+h-i}\right)\\ \\
            &=\mathbb{E}\left((\varepsilon_{t}+\phi\varepsilon_{t-1}+\phi^{2}\varepsilon_{t-2}+\cdots)\cdot(\varepsilon_{t+h}+\phi\varepsilon_{t+h-1}+\phi^{2}\varepsilon_{t+h-2}+\cdots+\phi^{h}\varepsilon_{t}+\phi^{h+1}\varepsilon_{t-1}+\cdots)\right)
\end{align}

dado que $\varepsilon_{t}\sim RB(0,\sigma_{\varepsilon}^{2})$, tenemos que los momentos cruzados $\mathbb{E}(\varepsilon_{t}\varepsilon_{s})=0$ sii $t\not=s$, luego

\begin{align}
\mathbb{C}(Y_{t},Y_{t+h})&=\mathbb{E}\left((\varepsilon_{t}+\phi\varepsilon_{t-1}+\phi^{2}\varepsilon_{t-2}+\cdots)\cdot(\varepsilon_{t+h}+\phi\varepsilon_{t+h-1}+\phi^{2}\varepsilon_{t+h-2}+\cdots+\phi^{h}\varepsilon_{t}+\phi^{h+1}\varepsilon_{t-1}+\cdots)\right)\\ \\
&=\sigma_{\varepsilon}^{2}\phi^{h}\sum_{i=0}^{\infty}\phi{2i}=\frac{\sigma_{\varepsilon}^{2}\phi^{h}}{1-\phi^{2}},\quad h>0
\end{align}

### Simulación $AR(1)$


```{r Simulación de un proceso AR(1)}
set.seed(1)             #Vamos a fijar la semilla en 1 para generar un mismo vector de números aleatorios

X <- matrix(rep(0,50),ncol=1)       #Generamos este vector con 0 para guardar los datos de una serie
inicio <- 5                         #Valor de la primera observación de la serie
phi <- 0.8                          #Parámetro phi de la ecuación (5.1)
epsilon <- rnorm(49,0,1.5)          #Vamos a modelar el error de la regresión con una distribución normal con media 0 y varianza 1.5
X[1]<- inicio                       #Inicializamos la serie con el valor inicio.

for (i in 2:nrow(X)){
  X[i]=X[i-1]*phi+epsilon[i-1]
}                                   #Generamos el resto de las observaciones

X <- ts(X,frequency = 1, 
        start = 1971)               #Formato de serie temporal con observaciones por año
```

```{r}
library(ggplot2)
X_df <- data.frame(fecha = time(X),
                       valor = as.vector(X)) #ggplot2 utiliza data frames, por lo tanto modificaré la serie original.

ggplot(X_df, aes(x = fecha, y = valor)) +
  geom_line() +
  labs(title = "Mi Serie de Tiempo", 
       x = "Fecha", y = "Valor")            #Grafica la serie
```

### Pronostico $AR(1)$

Sea la serie de tiempo $Y_{t}$ dada por
$$Y_{t}=\phi_{0}+\phi_{1}Y_{t-1}+\varepsilon_{t}$$
$$\varepsilon_{t}\sim N(0,\sigma^{2})$$
por lo tanto $Y_{t}\sim N(\phi_{0}+\phi_{1}Y_{t-1},\sigma^{2})$ ya que la serie se conoce hasta $t-1$, y $Y_{t-1}\in\mathbb{R}$. Para realizar un pronostico en $t+1$ se puede escribir como
$$Y_{t+1}=\phi_{0}+\phi_{1}Y_{t}+\varepsilon_{t+1}$$
$$\varepsilon_{t+1}\sim N(0,\sigma^{2})$$
Entonces,

\begin{align}
Y_{t+1}&=\phi_{0}+\phi_{1}(\phi_{0}+\phi_{1}Y_{t-1}+\varepsilon_{t})+\varepsilon_{t+1}\\ \\
      &=\phi_{0}+\phi_{0}\phi_{1}+\phi_{1}^{2}Y_{t-1}+\phi_{1}\varepsilon_{t}+\varepsilon_{t+1}
\end{align}

Por lo tanto, 
$$Y_{t+1}\sim N(\phi_{0}+\phi_{0}\phi_{1}+\phi_{1}^{2}Y_{t-1},\sigma^{2}(1+\phi^{2}))$$
```{r}
epsilon.for<-rnorm(1,0,1.5)
X.for<-X[length(X)]*phi+epsilon.for
X.for
```

### Serie Filtrada $AR$

El filtrado de la serie se expresa como 
$$Y_{t}-Y_{t-1}=\varepsilon_{t}$$
Donde la serie filtrada resulta ser el error del proceso.
```{r}
Xf <- c(rep(0,length(X)-1))          #Vector donde se almacenan serie filtrada

for (i in 2:length(Xf)){
  Xf[i]=X[i]-X[i-1]*phi
}
Xf <- ts(Xf, frequency = 1, start = 1071)
```

```{r}
Xf_df <- data.frame(fecha = time(Xf),
                    valor = as.vector(Xf))

ggplot(Xf_df,aes(x = fecha, y = valor))+
  geom_line()+
  labs(title = "Gráfico de la Serie Filtrada", 
       x = "Fecha", y = "Valor") 
```

### Estimación de procesos $AR$

ar.ols estima un modelo autorregresivo por el metodo de minimos cuadrados ordinarios (ols por sus siglas en ingles)

```{r}
ar <- ar.ols(
  X,                    #Datos a modelar
  aic = TRUE,           #Criterio AIC
  order.max =1,         #Proceso maximo AR(1)
  intercept = FALSE,    #Intercepto falso
  demean = FALSE        #Omite un valor medio que automaticamente calcula ar.ols
)
ar
```

La estimación del parámetro $\phi$ del modelo $AR(1)$ por medio de la función ar.ols fue de $\hat{\phi}=0.7289$ y el parámetro con el que se genero la serie era de $\phi=0.8$.

Para calcular los valores ajustados del modelo autorregresivo que se obtuvo anteriormente, se puede realizar de manera manual o por medio de la función fitted como se presenta a continuación:

```{r}
X.fit <- fitted(ar)
```
Grafica de los valores reales contra los estimados:
```{r}
X.fit_df <- data.frame(fecha = time(X),
                    valor = as.vector(X.fit))
ggplot() +
  geom_line(data = X.fit_df, aes(x = fecha, y = valor, color = "Serie Estimada")) +
  geom_line(data = X_df, aes(x = fecha, y = valor, color = "Serie Real")) +
  labs(title = "Serie Real vs. Serie Estimada",
       x = "Fecha",
       y = "Valor") +
  scale_color_manual(values = c("Serie Estimada" = "blue", "Serie Real" = "red"))
```

#### Cálculo de los residuales:
```{r}
residuals <- ar$resid #es la diferencia entre los valores observados y los estimados
residuals
```
Gráfico de los residuos
```{r}
residuals_df <- data.frame(residuals,
                           fecha = time(X),
                           valor = as.vector(residuals))
ggplot(residuals_df, aes(x = fecha, y = valor))+
  geom_line()+
  labs(title = "Gráfico de los Residuos", 
       x = "Fecha", y = "Valor") 
```

#### Diagnostico de los residuos:

La distribución normal de los residuos es una suposición fundamental detrás de muchos métodos de estimación de parámetros, como la estimación de máxima verosimilitud (MLE) y los mínimos cuadrados ordinarios (OLS). Cuando los residuos se distribuyen normalmente, las estimaciones de los parámetros del modelo son eficientes y no sesgadas.

La normalidad de los residuos es importante para hacer pronósticos y predicciones precisas. Los intervalos de predicción y las predicciones puntuales se basan en la suposición de que los residuos son normalmente distribuidos.

Para hacer el análisis de normalidad de los residuos, es necesario calcular su densidad primero:
```{r}
# Distribución de los residuos
residuals <- na.omit(residuals)
densidad_res <-density(residuals)
densidad_res_df <- data.frame(x = densidad_res$x,
                              y = densidad_res$y)

# Distribución de la normal con 1000 puntos
normal <- seq(min(residuals)-2, max(residuals)+2, length.out = 1000)
media <- mean(residuals)
desviacion <- sd(residuals)
densidad_normal <- dnorm(normal, mean = media, sd = desviacion)
normal_df <- data.frame(x = normal, densidad = densidad_normal)


#Gráfico
ggplot(densidad_res_df, aes(x = x, y = y)) +
  geom_line(color = "blue", linetype = "solid", size = 1) +
  geom_line(data = normal_df, aes(x = x, y = densidad), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Densidad Observaciones vs. Densidad Normal",
       x = "Valor",
       y = "Densidad") +
  scale_color_manual(values = c("blue", "red"))
```
#### Prueba ACF y Correlograma de los residuos:

La ACF es una función estadística que se utiliza para medir y visualizar la autocorrelación entre una serie temporal y sus propios rezagos. En otras palabras, la ACF muestra cómo se relaciona cada observación con las observaciones anteriores en diferentes períodos de tiempo. La ACF se calcula mediante la correlación entre la serie temporal original y versiones desplazadas de sí misma a lo largo del tiempo. La ACF es una herramienta importante para identificar patrones de autocorrelación en una serie temporal, como ciclos, estacionalidad y tendencias.

```{r}
par(mar = c(5, 5, 4, 2))
residuals_acf <- acf(residuals, lag.max = 20,
                     main = "Correlograma de los Residuos", cex.main = 1.2)

```

La falta de autocorrelación después del primer período sugiere que las observaciones de la serie temporal son esencialmente independientes entre sí. Esto significa que no hay una relación sistemática y significativa entre las observaciones en momentos sucesivos, lo que es una suposición importante en muchas técnicas estadísticas y de modelado.

Además, la ausencia de autocorrelación significativa en los rezagos posteriores puede indicar que el comportamiento de la serie es estable y no muestra patrones de estacionalidad, tendencia o ciclos que se repiten a lo largo del tiempo.

#### Prueba PACF
La PACF se utiliza comúnmente para determinar el orden adecuado de un modelo AR. Cuando se observa un decaimiento significativo en la PACF después de un rezago específico, indica que ese rezago es el orden adecuado del modelo AR. Por ejemplo, si la PACF muestra un decaimiento significativo después del rezago 2, sugiere un modelo AR(2). 

La PACF puede ayudar a diferenciar entre modelos AR y modelos MA (que veremos en el siguiente apartado). Si la PACF muestra una caída abrupta después de un rezago específico y luego permanece cerca de cero, sugiere un modelo AR, mientras que una PACF que decae lentamente sugiere un modelo MA. 

```{r}
par(mar = c(5, 5, 4, 2))
residuals_pacf <- pacf(residuals, lag.max = 20,
                     main = "PACF de los Residuos", cex.main = 1.2)
```


En el contexto de modelado de series temporales, si la PACF de los residuos está dentro del intervalo de confianza en todos los rezagos, indica que el modelo elegido es adecuado para describir la estructura de dependencia en los datos. Los residuos parecen no contener información adicional que no se haya capturado en el modelo.


# Procesos de Medias Moviles $MA(q)$

$\textbf{Definición 2.1}$ Se dice que una serie temporal $Y_{t}$ admite una representación de medias móviles (MA) de orden $q$ , y se denota por $MA(q)$ , si es susceptible de ser modelizada a través de la ecuación:
\begin{align}
Y_{t}&=\theta_{1}\varepsilon_{t-1}+\theta_{2}\varepsilon_{t-2}+\cdots+\theta_{q}\varepsilon_{t-q}+\varepsilon_{t}\\ \\
    &=\sum_{i=1}^{q}\theta_{i}\varepsilon_{t-i}+\varepsilon_{t},  \tag{1.2}
\end{align}

Donde $Y_{t}$ es una variable aleatoria concebida como realizaciones de un proceso estocástico en los momentos del tiempo $t$ , que se caracterizan por $\mathbb{E}(Y_{t})=\mathbb{E}(Y_{t−1})=\mathbb{E}(Y_{t−2})=\cdots$, $\theta_{1},...,\theta_{q}$ son parámetros del modelo que, junto con la varianza deben ser estimados, y $\varepsilon_{t}\sim RB(0,\sigma_{\varepsilon}^{2})$.

$\textbf{Definición 2.1}$ Se dice que un proceso MA(1) es $\textit{estacionario en media}$ si
$$\mathbb{E}(Y_{t})=\mathbb{E}(\varepsilon_{t}+\theta\varepsilon_{t-1})=0$$
$\textbf{Definición 2.3}$ Se dice que un proceso MA(1) es $\textit{estacionario en covarianza}$ si

\begin{align}
\gamma_0&=\mathbb{E}(Y_{t}-\mathbb{E}(Y_{t}))^{2}=\mathbb{E}(Y_{t})^{2}=\mathbb{E}(\varepsilon_{t}+\theta\varepsilon_{t-1})^{2}\\\\
      &=\mathbb{E}(\varepsilon_{t}^{2}+2\theta^{2}\varepsilon_{t}\varepsilon_{t-1}+\theta^{2}\varepsilon_{t-1}^{2})\\\\
      &=\mathbb{E}(\varepsilon_{t}^{2})+2\theta^{2}\mathbb{E}(\varepsilon_{t}\varepsilon_{t-1})+\mathbb{E}(\varepsilon_{t-1}^{2})\\ \\
      &=\sigma_{\varepsilon}^{2}+\sigma_{\varepsilon}^{2}\theta^{2}=\sigma_{\varepsilon}^{2}(1+\theta^{2})<\infty
\end{align}

### Simulación $MA(1)$

```{r}
set.seed(1)

Y <- matrix(rep(0,50), ncol=1)
epsilon <- rnorm(50,0,1.5)
theta <-0.9
cte <- 0

Y[1] <- epsilon[1]

for (i in 2:nrow(Y)){
  Y[i] <- cte + theta*epsilon[i-1] + epsilon[i]
}

Y <- ts(Y, frequency = 1, start = 1970)
```


```{r}
Y_df <- data.frame(fecha = time(Y),
                   valor = as.vector(Y))

ggplot(Y_df,aes(x=fecha, y=valor))+
  geom_line()+
  geom_point()+
  labs(title = "Gráfico de MA", xlab("Fecha"), ylab("Valor"))
```

### Pronóstico con Proceso $MA$

Para hacer la estimación de la serie vamos a suponer que,
$$\hat{Y_{t}}=\phi_{0}+\phi_{1}e_{t-1}+e_{t}$$
$$e_{t}\sim\sigma_{\varepsilon}\varepsilon_{t}$$

$$\varepsilon\sim N(0,1)$$
donde $Y_{t}$ es el último valor conocido de la serie y $e_{t}$ es el último valor del ruido conocido. 

Suponga que $\phi_{0}=0$ $Y_{0}=e_{t}$, $\hat{Y_{1}}=\phi_{1}e_{0}+e_{1}$, $\hat{Y_{2}}=\phi_{1}e_{1}+e_{2}$, por lo tanto
$$\hat{Y}_{t+1}=\phi_{1}e_{t}+e_{t+1}$$
### Serie Filtrada $MA(1)$
```{r}
Yf <- matrix(rep(0,50),ncol=1)

for (i in 2:nrow(Yf)){
  Yf[i-1] <- Y[i]+epsilon[i-1]
}

Yf <- ts(Yf, frequency = 1, start = 1970)
```

```{r}
Yf_df <- data.frame(fecha = time(Yf),
                    valor = as.vector(Yf))

ggplot(Yf_df, aes(x = fecha, y = valor))+
  geom_line()+
  geom_point()+
  labs(title = "Gráfico Serie Filtrada MA", xlab("Fehca"),
       ylab("Valor"))
```

### Estimación de Proceso $MA$
```{r}
ma <- arima(Y,
          order = c(0,0,2),   #(p=0,i=0,q=1) es un proceso MA(1)
          include.mean = F)   
ma
```
En la presentación de los resultados de la estimación del modelo, el termino s.e significa error estandar y se refiere a una estimación de la desviación estandar.
$$SE_{\hat{Y}}=\frac{s}{\sqrt{n}}$$
Donde $s$ es la desviación estandar  y $n$ el número de observaciones.



```{r}
Y.fit <- fitted(ma)

Y.fit_df <- data.frame(fecha = time(X),
                    valor = as.vector(Y.fit))
ggplot() +
  geom_line(data = Y.fit_df, aes(x = fecha, y = valor, color = "Serie Estimada")) +
  geom_line(data = Y_df, aes(x = fecha, y = valor, color = "Serie Real")) +
  labs(title = "Serie Real vs. Serie Estimada",
       x = "Fecha",
       y = "Valor") +
  scale_color_manual(values = c("Serie Estimada" = "blue", "Serie Real" = "red"))
```

#### Calculo de los residuales del modelo MA:

El modelo ajustado por medio de la función arima y asignado a la variable "ma" tienen varios parametros internos. Para poder ver cual es la estructura interna y su composición utilice str( ). Los residuales del modelo ya estan calculados internamente y se calculan:
```{r}
residuals <- ma$residuals
residuals_df <- data.frame(fecha = time(residuals),
                           valor = as.vector(residuals))

ggplot(residuals_df, aes(x = fecha, y = valor))+
         geom_line()+
         geom_point()+
         labs(title = "Gráfico de los residuos del proceso MA",
              xlab("Fecha"), ylab("Valor"))
```

#### Diagnostico de los residuales MA:

```{r}
# Distribución de los residuos
residuals <- na.omit(residuals)
densidad_res <-density(residuals)
densidad_res_df <- data.frame(x = densidad_res$x,
                              y = densidad_res$y)

# Distribución de la normal con 1000 puntos
normal <- seq(min(residuals)-2, max(residuals)+2, length.out = 1000)
media <- mean(residuals)
desviacion <- sd(residuals)
densidad_normal <- dnorm(normal, mean = media, sd = desviacion)
normal_df <- data.frame(x = normal, densidad = densidad_normal)


#Gráfico
ggplot(densidad_res_df, aes(x = x, y = y)) +
  geom_line(color = "blue", linetype = "solid", size = 1) +
  geom_line(data = normal_df, aes(x = x, y = densidad), color = "red", linetype = "dashed", size = 1) +
  labs(title = "Densidad Observaciones vs. Densidad Normal",
       x = "Valor",
       y = "Densidad") +
  scale_color_manual(values = c("blue", "red"))
```

#### Prueba ACF y Correlograma de los residuos:
```{r}
par(mar = c(5, 5, 4, 2))
acf(residuals, 
    na.action = na.pass,
    ylim=c(-1,1),
    main = "Corelograma de los Residuos")
```

#### Prueba PACF:

```{r}
par(mar = c(5, 5, 4, 2))
pacf1 <- pacf(residuals,
             na.action = na.pass,
             ylim=c(-1,1),
             main = "PACF de los residuos")
```

## Procesos autorregresivos de media movil $ARMA$

Los modelos $ARMA(p,q)$ son modelos autorregresivos de media movil se componen de una parte AR y la otra MA, siguen la siguiente ecuación:

\begin{align}
Y_{t}&=\phi_{0}+\phi_{1}Y_{t-1}+\cdots+\phi_{p}Y_{t-p}+\varepsilon_{t}+\theta_{1}\varepsilon_{t-1}+\cdots+\theta_{q}\varepsilon_{t-q} \\ \\
    &=\phi_{0}+\sum_{i=1}^{p}\phi_{i}Y_{t-i}+\sum_{j=1}^{q}\theta_{j}\varepsilon_{t-j} + \varepsilon_{t}\tag{3.1}
\end{align}

### Simulación del proceso $ARMA$

```{r}
set.seed(1)

Z <- arima.sim(n = 50,                         #número de observaciones simuladas
               list(ar = c(0.82,-0.4),         #parámetro phi del modelo AR
                    ma = c(-0.2,0.24)),        #parámetro theta del modelo MA
               innov=rnorm(50,0,1))            #error asosciado distribuido normalmente                  
```

```{r}
Z_df<-data.frame(time = c(1970:2019),
                 value = as.vector(Z))

ggplot(Z_df,aes(x = time,
               y = value))+
  geom_line()+
  geom_point()+
  labs(title = "Simulación modelo ARMA")
```


```{r}
#Vamos a ajustar con un modelo ARMA(2,2)

par(mar = c(5, 5, 4, 2))
acf(Z, 
    na.action = na.pass,
    ylim=c(-1,1))
```
```{r}
Z.fit <- arima(Z, order = c(2,0,2),
               include.mean = FALSE)
Z.fit
```
```{r}
Z.fit_df <- data.frame(time = c(1970:2019),
                       value = as.vector(fitted(Z.fit)))

ggplot() +
  geom_line(data = Z.fit_df, aes(x = time, y = value, color = "Serie Estimada")) +
  geom_line(data = Z_df, aes(x = time, y = value, color = "Serie Real")) +
  labs(title = "Serie Real vs. Serie Estimada",
       x = "Fecha",
       y = "Valor") +
  scale_color_manual(values = c("Serie Estimada" = "blue", "Serie Real" = "red"))
```
```{r}
Z.for <- predict(Z.fit,n.ahead=2)
Z.for
```
```{r}
#Actualización de la serie:
pronostico <- data.frame(c(2020,2021),as.vector(Z.for$pred))
colnames(pronostico)<-colnames(Z.fit_df)
Z.fit_df<-rbind(Z_df,pronostico)

ggplot() +
  geom_line(data = Z.fit_df, aes(x = time, y = value, color = "Pronostico")) +
  geom_line(data = Z_df, aes(x = time, y = value, color = "Serie Real")) +
  labs(title = "Serie Real y Pronóstico",
       x = "Fecha",
       y = "Valor") +
  scale_color_manual(values = c("Pronostico" = "blue", "Serie Real" = "red"))

```

### Diferenciación de la serie de tiempo

Para utilizar esta metodología para el modelamiento y posterior pronóstico de series de tiempo se debe verificar la estacionaridad de la serie. Esta se puede diagnosticar mediante el test de Dickey-Fuller, si la serie resulta no estacionaria se debe diferenciar la serie hasta alcanzar la estacionariedad.

Sea $Y_{t}$ una serie de tiempo, su diferenciación simple esta dada por:
$$\Delta Y_{t}=Y_{t}-Y_{t-1}$$
```{r}
# Prueba Dickey-Fuller de la serie orignial
adf.test(Y)

## Ejemplo de diferenciación simple
Y.d <- diff(Y)  # Serie X con presencia de raiz Unit.
adf.test(Y.d)
```
La serie era estacionaria desde el pricipio por eso el valor $p$ en ambas pruebas no cambia.






## Proceso autorregresivo integrado de media móvil $ARIMA$

El modelo autorregresivo integrado de media movil utiliza diferenciación de la serie y regresiones para modelar la serie de tiempo, se considera que es un modelo dinamico ya que los datos de predicciones dependen de comportamientos pasados y no de variables externas a la serie.

El proceso $ARIMA(p,d,q)$ donde $p,d,q$ son números enteros no negativos, se puede representar matematicamente como:



$$Y_{t}=-(\Delta^{d}Y_{t}-Y_{t})+\phi_{0}+\sum_{i=1}^{p}\phi_{i}\Delta^{d}Y_{t-i}+\sum_{j=1}^{q}\theta_{j}\varepsilon_{t-j}+\varepsilon_{t}$$
Donde:

* $\phi_{0}$ es una constante
* $d$ corresponde a las diferencias necesarias para generar que la serie sea estacionaria
* $\phi_{1},...,\phi_{p}$ son los parametros asociados a la parte autorregresiva del modelo
* $\theta_{1},...,\theta_{q}$ son los parametros asociados a la parte de medias moviles del modelos
* $\varepsilon_{t}$ es el termino de error.


### Simulación proceso $ARIMA$


```{r}
set.seed(1)

W <- arima.sim(n = 50,                        #50 observaciones
               list(order = c(1,1,1),         #ARIMA(1,1,1) AR(1), Dif(1), MA(1)
                    ar = 0.812,               #Phi-AR
                    ma = 0.214))              #Theta-MA
head(W)

W_df <- data.frame(time = c(1970:2020),
                   value = as.vector(W))

ggplot(W_df, aes(x = time, y = value))+
  geom_line()+
  geom_point()+
  labs(title = "Simulación proceso ARIMA(1,1,1)", 
       xlab("Año"), ylab("Valor"))
```

### Estimacion modelo ARIMA
```{r}
W.fit <- arima(W,                    # Datos de la serie
               order = c(1,1,1),     # Proceso tipo ARIMA(1,1,1) ar(1), diff(1),ma(1)
               include.mean = T)  
W.fit

W.fit_df <- data.frame(time = W_df$time,
                       value = as.vector(fitted(W.fit)))

ggplot() +
  geom_line(data = W.fit_df, aes(x = time, y = value, color = "Serie Estimada")) +
  geom_line(data = W_df, aes(x = time, y = value, color = "Serie Real")) +
  labs(title = "Serie Real vs. Serie Estimada",
       x = "Fecha",
       y = "Valor") +
  scale_color_manual(values = c("Serie Estimada" = "blue", "Serie Real" = "red"))
```
### Predicción del modelo $Arima(1,1,1)$

```{r}
W.for<-predict(W.fit,n.ahead=3)
W.for

pronostico <- data.frame(c(2021:2023),
                         as.vector(W.for$pred))

colnames(pronostico) <- colnames(W.fit_df)
W.fit_df <- rbind(W.fit_df,pronostico)

ggplot() +
  geom_line(data = W.fit_df, aes(x = time, y = value, color = "Pronostico")) +
  geom_line(data = W_df, aes(x = time, y = value, color = "Serie Real")) +
  labs(title = "Serie Real y Pronóstico",
       x = "Fecha",
       y = "Valor") +
  scale_color_manual(values = c("Pronostico" = "blue", "Serie Real" = "red"))
```


<br><br><br><br><br><br><br><br><br><br><br><br>
<p style="text-align:center;">
*Cuidate bb*
</p>


