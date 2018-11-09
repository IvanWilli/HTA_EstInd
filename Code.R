############################## Hipertensión Arterial en Argentina ######################################
############################## Transición y Mortalidad Diferencial #####################################
################# Una estimación indirecta a partir de datos transversales (2009-2013) #################

### Iván Williams (Jun2017)
# 1) Se preparan los datos
# 2) Supuestos y definiciones previas
# 3) Optimización
# 4) Resultados

#1)#############################################################Preparación de datos

#levanto data 
source("Prev+Mort+Canada.R")


#####Prevalencia observada
prevObs <- data.frame(prev_09_18a88$prev,prev_13_18a88$prev)

#####Prev ajustadas
#Defino edad de ajuste
x_aj = c(30,70)
prev3070 <- data.frame(x=30:70,
                     prev09=subset(prev_09_18a88, x>=x_aj[1] & x<=x_aj[2])$prev,
                     prev13=subset(prev_13_18a88, x>=x_aj[1] & x<=x_aj[2])$prev)

#Exponencial por cada año
#09
fitExp_09 <- lm(data= prev3070, log(prev09)~x+1)
    summary(fitExp_09)$r.squared
#13
fitExp_13 <- lm(data= prev3070, log(prev13)~x+1)
    summary(fitExp_13)$r.squared
#Graf0913
    plot(prev3070$x,prev3070$prev09,col=2,main="Prevalencias suavizadas",xlab="Edad",ylab="Prevalencia")
      points(prev3070$x,prev3070$prev13,col=3)
      lines(prev3070$x,exp(predict(fitExp_09)),col=2) 
      lines(prev3070$x,exp(predict(fitExp_13)),col=3)
      
###prevalencia final suavizada
x_suav = 30:70
prev09_ic <- predict(fitExp_09,newdata=data.frame(x=x_suav),interval='confidence',level=0.95,se.fit = TRUE)
prev13_ic <- predict(fitExp_13,newdata=data.frame(x=x_suav),interval='confidence',level=0.95,se.fit = TRUE)
prev_ic <- data.frame(x=x_suav,"09"=prev09_ic$fit,"13"=prev13_ic$fit)

###mortalidad
#la interpolación de la tabla de mortalidad nacional a 2011 no es realizada aquí (no estimado aquí, pronto a incluir)

#spline para desagregar por edad simple
library("splines")
spl11 <- interpSpline(mort_11$x,mort_11$l)
edades <- data.frame(x = 30:70)
mort11 <- data.frame(x=edades,q=predict(spl11,edades$x))

#función para construir q a partir de l
qs<-function(s)
{   newq<-data.frame(x = 30:66,y=seq(1,37))
    for (i in 1:(nrow(s)-4))
    {
        newq[i,2]<-(s[i,3]-s[i+4,3])/s[i,3]
    }
    newq}
q11<-qs(mort11)

# 2)#####################################supuestos y definiciones previas
###set de bandas para rr y q(HTA)
#rr
    rr_li<-c(2,1)
    rr_ls<-c(6,2.5)
    x_l<-c(30,66)
    dat_li=data.frame(rr_li,x_l)
    dat_ls=data.frame(rr_ls,x_l)
    model_li <- lm(log(rr_li) ~ x_l + 1) 
    model_ls <- lm(log(rr_ls) ~ x_l + 1)
    x <- data.frame(x_l=seq(30,66,4))
    pred_li_rr<-data.frame(x=seq(30,66,4),p=exp(predict(model_li,x)))
    pred_ls_rr<-data.frame(x=seq(30,66,4),p=exp(predict(model_ls,x)))
    a1.li<-exp(model_li$coefficients[1])
    b1.li<-model_li$coefficients[2]
    a1.ls<-exp(model_ls$coefficients[1])
    b1.ls<-model_ls$coefficients[2]
    Can_df_<-Can_df[3:10,]
    model_<- lm(data = Can_df_,log(rr) ~ x + 1)
    x_ <-data.frame(x=seq(30,66,4))
    pred_rr<-data.frame(x=seq(30,66,4),
                      p=exp(predict(model_,x_)))
#q(HTA)
    inc_li<-c(1,10)/100
    inc_ls<-c(10,60)/100
    x_l<-c(32.5,67.5)
    model_li <- lm(log(inc_li) ~ x_l + 1) 
    model_ls <- lm(log(inc_ls) ~ x_l + 1)
    x <- data.frame(x_l=seq(30,66,4))
    pred_li_inc<-data.frame(x=seq(30,66,4),p=exp(predict(model_li,x)))
    pred_ls_inc<-data.frame(x=seq(30,66,4),p=exp(predict(model_ls,x)))
    a2.li<-exp(model_li$coefficients[1])
    b2.li<-model_li$coefficients[2]
    a2.ls<-exp(model_ls$coefficients[1])
    b2.ls<-model_ls$coefficients[2]
    Can_df_<-Can_df[3:10,]
    model_<- lm(data = Can_df_,log(inc*4/100) ~ x + 1)
    x_ <-data.frame(x=seq(30,66,4))
    pred_inc<-data.frame(x=seq(30,66,4),
                      p=exp(predict(model_,x_)))
#gráfico de bandas
    library("ggplot2")
    incgraf<-ggplot() +
        geom_point(data = Can_df, aes(x=x, y=inc/100)) +
        geom_line(data = pred_inc, aes(x=x, y=p/4))+
        scale_y_continuous(name = "i(x)", limits=c(0,0.1))+
        scale_x_continuous(name = "Edad", breaks = seq(20,90,10))+
        ggtitle("Tasa de Incidencia")+
        theme(plot.title = element_text(size = 10),
              axis.text=element_text(size=9),
              axis.title=element_text(size=10))
    rrgraf<-ggplot() +
        geom_point(data = Can_df, aes(x=x, y=rr)) +
        geom_line(data = pred_rr, aes(x=x, y=p)) +
        scale_y_continuous(name = "rr(x)", limits=c(1,5))+
        scale_x_continuous(name = "Edad", breaks = seq(20,90,10))+
        ggtitle("Sobremortalidad")+
        theme(plot.title = element_text(size = 10),
              axis.text=element_text(size=9),
              axis.title=element_text(size=10))
    library("gridExtra")
    grid.arrange(incgraf,rrgraf,ncol=2, 
                 top="Valores observados y suavización exponencial de edades 30-70. Canadá. Años 2007/8")

###Formas funcionales por edad (se le llama inc a q(HTA) erróneamente)
    rr<-function(x) {a1*exp(x*b1)}
    der_rr<-function(x) {a1*b1*exp(x*b1)}
    inc<-function(x) {a2*exp(x*b2)}
    der_inc<-function(x) {a2*b2*exp(x*b2)}


#3)##########################################OPtimización
#función de optimización que toma como argumentos
 #prev al 09 y su desvío standar por edad
 #prev al 13 y su desvío standar por edad

auglag_funcion<-function(prev1, prev2, prev1.se.fit, prev2.se.fit, q){

        #sorteo prevalencias de cada año
        prev1.rnd=data.frame(x=30:70,prev=rnorm(n=length(prev1),mean = prev1, sd = prev1.se.fit))
        prev1.rnd=lm(data=prev1.rnd,log(prev)~x+1)
        prev1.rnd=exp(predict(prev1.rnd,data.frame(x=30:70)))
        prev2.rnd=data.frame(x=30:70,prev=rnorm(n=length(prev2),mean = prev2, sd = prev2.se.fit))
        prev2.rnd=lm(data=prev2.rnd,log(prev)~x+1)
        prev2.rnd=exp(predict(prev2.rnd,data.frame(x=30:70)))
        
        #complemento (sanos)
        Noprev1<-1-prev1.rnd
        Noprev2<-1-prev2.rnd
        
        #agregación en letras según GY(2009)
        Y<-c();B<-c();C<-c();D<-c()
        for (i in 1:length(q)){
            Y[i]<-Noprev2[i+4]-Noprev1[i]/(1-q[i])
            B[i]<-Noprev1[i]*q[i]/(1-q[i])
            C[i]<-Noprev1[i]
            D[i]<-Noprev1[i]/(1-q[i])}
        
        #edades de cohortes al 2009
        x<-30:66
        
        #función objetivo
        fn<-function(p){
            sum((Y+B/(C+(1-C)*p[1]*exp(x*p[2]))+D*p[3]*exp(x*p[4]))^2)
            }
        
        #sistema de restricciones de desigualdad (bandas)
        hin<-function(p){
            h <- rep(NA, 68)
            cont<-1
            for (k in seq(30,66,4))
            {
                #rr
                h[cont] <- p[1]*exp(k*p[2])-a1.li*exp(k*b1.li)
                cont<-cont+1
                h[cont] <- a1.ls*exp(k*b1.ls)-p[1]*exp(k*p[2])
                cont<-cont+1
                h[cont] <- -p[1]*p[2]*exp(k*p[2])
                cont<-cont+1
                
                #inc
                h[cont] <- p[3]*exp(k*p[4])-a2.li*exp(k*b2.li)
                cont<-cont+1
                h[cont] <- a2.ls*exp(k*b2.ls)-p[3]*exp(k*p[4])
                cont<-cont+1
                h[cont] <- p[3]*p[4]*exp(k*p[4])
                cont<-cont+1
                
                #parámetros
                h[61]<-p[1]-a1.li
                h[62]<-a1.ls-p[1]
                h[63]<-p[2]-b1.li
                h[64]<-b1.ls-p[2]
                h[65]<-p[3]-a2.li
                h[66]<-a2.ls-p[3]
                h[67]<-p[4]-b2.li
                h[68]<-b2.ls-p[4]
            }
            return(h)
        }
        
        #restricciones de igualdad (rr=1.1 en edades finales)
        heq<-function(p){
            h<-rep(NA,1)
            h[1]<-p[1]*exp(100*p[2])-1.1
            h
        }
        
        #loopeo soluciones
        soluciones<-data.frame(a1=NA,b1=NA,a2=NA,b2=NA,c2=NA,con=NA, val=NA, err=NA)
        
        m<-30:66 #edades de chekeo de restr
        a1<-c();b1<-c();a2<-c();b2<-c()
        for (i in 1:1){ #si se quiere se puede randomizar la solución inicial. Converge
            randoma1<-runif(1)
            randomb1<-runif(1)
            randoma2<-runif(1)
            randomb2<-runif(1)
            a1_i<-randoma1*a1.li+(1-randoma1)*a1.ls
            b1_i<-randomb1*b1.li+(1-randomb1)*b1.ls
            a2_i<-randoma2*a2.li+(1-randoma2)*a2.ls
            b2_i<-randomb2*b2.li+(1-randomb2)*b2.ls
            iniciales<-c(a1_i,b1_i,a2_i,b2_i)
        
            #optim
        solucion <- auglag(par=iniciales, fn=fn, hin=hin, heq=heq,
                               control.outer = list(itmax= 300, method="nlminb"))
            #control de soluciones
            soluciones[i,1]<-solucion$par[1]
            soluciones[i,2]<-solucion$par[2]
            soluciones[i,3]<-solucion$par[3]
            soluciones[i,4]<-solucion$par[4]
            soluciones[i,5]<-0
            soluciones[i,6]<-solucion$convergence
            soluciones[i,7]<-solucion$value
            soluciones[i,8]<-0
            a1<<-soluciones[i,1];b1<<-soluciones[i,2]
            a2<<-soluciones[i,3];b2<<-soluciones[i,4]
            for (j in m) {
                if (inc(j)<0 | inc(j)>1 | der_inc(j)<0 | 
                    rr(j)<1 | rr(j)>10 | der_rr(j)>0)
                {soluciones[i,8]<-1}
            }
        }
        solOK<-subset(soluciones, soluciones$err==0)
        return(solOK)
}

# 4)###########################################Resultados: estimación 0913

library("alabama")
set.seed(2378)
solOk_0913<-data.frame(a1=NA,b1=NA,a2=NA,b2=NA,c2=NA,
                       con=NA, val=NA, err=NA)
n_sim=20
for (i in 1:n_sim){
    sol0913<-auglag_funcion(exp(prev09_ic$fit[,1]), exp(prev13_ic$fit[,1]), 
                            prev09_ic$se.fit, prev13_ic$se.fit, q11[,2])
    solOk_0913[i,]<-head(sol0913[order(sol0913$val),],2)[1,]
}

#chek
sol0913=subset(solOk_0913,!is.na(err))
nrow(sol0913)

###las funciones resultantes de inc y rr son
inc_0913<-data.frame();rr_0913<-data.frame()
    for (i in 1:nrow(sol0913)){
        a1<-sol0913[i,1]; b1<-sol0913[i,2]
        a2<-sol0913[i,3]; b2<-sol0913[i,4]; c2<-sol0913[i,5]
        ord<-1
        for (j in 30:66) {
            inc_0913[i,ord]<-inc(j)
            rr_0913[i,ord]<-rr(j)
            ord<-ord+1}
    }

##Media, Mediana y percentiles de ambas funciones
quantfun97.5 <- function(x) (quantile(x, 0.975))
quantfun2.5 <- function(x) (quantile(x, 0.025))
inc_0913_CI<-data.frame(x=30:66,
                        media=apply(X=inc_0913, MARGIN = 2, FUN = mean)
                        , mediana=apply(X=inc_0913, MARGIN = 2, FUN = median)
                        , lwr=apply(X=inc_0913, MARGIN = 2, FUN = quantfun2.5)
                        , up=apply(X=inc_0913, MARGIN = 2, FUN = quantfun97.5)
)
rr_0913_CI<-data.frame(x=30:66,
                       media=apply(X=rr_0913, MARGIN = 2, FUN = mean)
                       , mediana=apply(X=rr_0913, MARGIN = 2, FUN = median)
                       , lwr=apply(X=rr_0913, MARGIN = 2, FUN = quantfun2.5)
                       , up=apply(X=rr_0913, MARGIN = 2, FUN = quantfun97.5)
)

###Gráfico de ajuste final
    prev09_g <- data.frame(x=30:70, exp(prev09_ic$fit[,1:3]))
    prev13_g <- data.frame(x=30:70, exp(prev13_ic$fit[,1:3]))
    Noprev09 <- 1-prev09_g$fit
    Noprev13 <- 1-prev13_g$fit
    Noprev13_estim <- data.frame()
    x = 34:70 
    library(ggplot2)
    graf_aj_0913 <- function(Noprev09,inc_0913,rr_0913,q11){
        graf <- ggplot()
        # j las simulaciones
        # i las edades
        for (j in 1:n_sim){
            #j=1
            Noprev13_estim[j,1]<-0
            Noprev13_estim[j,2]<-0
            Noprev13_estim[j,3]<-0
            Noprev13_estim[j,4]<-0
            for (i in 1:37){
            #i=3    
                Noprev13_estim[j,i+4]<-Noprev09[i]/(1-q11[i,2])-Noprev09[i]/(1-q11[i,2])*(q11[i,2]/(Noprev09[i]+(1-Noprev09[i])*rr_0913[j,i])+inc_0913[j,i])
            }
            data1=data.frame(edad=x,
                             prev=(rep(1,37)-t(Noprev13_estim[j,5:41]))*100,row.names = NULL)
            names(data1)<-c("edad","prev")
            graf<-graf+geom_line(data=data1,
                                    aes(x=edad,y=prev),color="red",alpha=0.1)
            graf<-graf+geom_point(data=data1,
                                  aes(x=edad,y=prev),color="red",alpha=0.1)
        }
        graf <- graf+geom_line(data=data.frame(edad=x[1:37],prev=exp(prev13_ic$fit[,1][5:41])*100),
                             aes(x=edad,y=prev),color="black",linetype="dashed")    
        graf <- graf+geom_point(data=prev_13_18a88[17:53,],
                              aes(x=x,y=prev*100),color="grey")
        return(graf)
    }
    grafaj0913 <- graf_aj_0913(Noprev09,inc_0913,rr_0913,q11)
    grafaj0913 <- grafaj0913+
        ggtitle("Prevalencia estimada por edad como reultado del modelo. Total del país. Año 2013")+labs(x="Edad",y="%")+
        theme(plot.title = element_text(size = 10),
          axis.text=element_text(size=9),
          axis.title=element_text(size=10))
    

###Cuadro y gráfico ppales de cada función
    x=30:66
    library(ggplot2)
    gInc0913 <- ggplot()+
        geom_line(data=as.data.frame(inc_0913_CI), aes(y=log(media*100),x=x), color="red")+
        geom_line(data=as.data.frame(inc_0913_CI), aes(y=log(mediana*100),x=x), color="red", alpha=0.8, linetype="12345678")+
        geom_line(data=as.data.frame(inc_0913_CI), aes(y=log(lwr*100),x=x), color="red", alpha=0.5, linetype="dashed")+
        geom_line(data=as.data.frame(inc_0913_CI), aes(y=log(up*100),x=x), color="red", alpha=0.5, linetype="dashed")+
        #scale_y_continuous(limits=c(0,30))+
        scale_x_continuous(breaks = seq(30,70,5))+ xlab("Edad")+ ylab("log(qH(x,4))")+
        ggtitle("Pobabilidad de transición a HTA en 4 años")+
        theme(plot.title = element_text(size = 10),
          axis.text=element_text(size=9),
          axis.title=element_text(size=10),
          legend.text=element_text(size=rel(0.8)),
          legend.title=element_text(size=rel(0.8)))
    gRR0913<-ggplot()+
        geom_line(data=as.data.frame(rr_0913_CI), aes(y=media,x=x), color="green")+
        geom_line(data=as.data.frame(rr_0913_CI), aes(y=mediana,x=x), color="green",alpha=0.8, linetype="12345678")+
        geom_line(data=as.data.frame(rr_0913_CI), aes(y=lwr,x=x), color="green", alpha=0.5, linetype="dashed")+
        geom_line(data=as.data.frame(rr_0913_CI), aes(y=up,x=x), color="green", alpha=0.5, linetype="dashed")+
        scale_y_continuous(limits=c(1,6))+
        scale_x_continuous(breaks = seq(30,70,5))+ xlab("Edad")+ ylab("rr(x,4)")+
        ggtitle("Sobremortalidad de la HTA  en 4 años")+
        theme(plot.title = element_text(size = 10),
              axis.text=element_text(size=9),
              axis.title=element_text(size=10),
              legend.text=element_text(size=rel(0.8)),
              legend.title=element_text(size=rel(0.8)))
    grid.arrange(gInc0913,gRR0913,ncol=2, 
                 top="Soluciones por edad. Total del país. Período 2009/13")
    
####esperanza de vida de HTA 09-11

#mortalidad resultante: HTA y no HTA
q11_NoHTA <- q11[,2]/(Noprev09[1:37]+(1-Noprev09[1:37])*rr_0913_CI$mediana)
q11_HTA <- q11_NoHTA*rr_0913_CI$mediana
    plot(x[1:37],q11_NoHTA, ylim=c(0,0.15))
    points(x[1:37],q11_HTA,col=2)
#obtener lx a partir de qx
tabla_l <- function(q)
{
    l<-rep(0,11)
    l[1]<-100000
    p<-1-q
    for (i in 1:10){
        l[i+1]<-p[i]*l[i]}
    l}

#funciones de sobrevivencia de cada subpoblación
l11_HTA<-tabla_l(q11_HTA[seq(1,37,4)])
l11_NoHTA<-tabla_l(q11_NoHTA[seq(1,37,4)])
l11<-tabla_l(q11[seq(1,37,4),2])
l<-data.frame(l11,l11_NoHTA,l11_HTA)
plot(seq(30,70,4),l11_NoHTA,ylim = c(0,100000))
  points(seq(30,70,4),l11_HTA,col=2)
  points(seq(30,70,4),l11,col=3)

#desagrego por edad simple con spline de lx
library("splines")
l11_HTA_spl <- interpSpline(seq(30,70,4),l11_HTA)
l11_HTA_splp <- data.frame(x=30:70,l=predict(l11_HTA_spl,30:70)[[2]])
l11_NoHTA_spl <- interpSpline(seq(30,70,4),l11_NoHTA)
l11_NoHTA_splp <- data.frame(x=30:70,l=predict(l11_NoHTA_spl,30:70)[[2]])
l11_spl <- interpSpline(seq(30,70,4),l11)
l11_splp <- data.frame(x=30:70,l=predict(l11_spl,30:70)[[2]])
#gráf    
    plot(30:70,l11_HTA_splp$l);points(seq(30, 70, 4),l11_HTA,col=2)
    lsplines<-data.frame(x=30:70,l11_splp[,2],l11_HTA_splp[,2],l11_NoHTA_splp[,2])  ###OK
    plot(l11_HTA_splp$x,l11_HTA_splp$l, col=2, type="l", ylim = c(50000,100000));
        points(l11_NoHTA_splp$x,l11_NoHTA_splp$l, col=1, type="l");
        points(l11_splp$x,l11_splp$l, col=3, type="l")


#construyo q simples
tabla_q<-function(l)
{
    q<-rep(1,40)
    for (i in 1:(length(l)-1)){
        q[i]<-(l[i]-l[i+1])/l[i]}
    q}
q_NoHTA<-tabla_q(l11_NoHTA_splp$l)
q_HTA<-tabla_q(l11_HTA_splp$l)
q_<-tabla_q(l11_splp$l)
    plot(30:66,q11_NoHTA,ylim = c(0,0.05));points(30:69,q_NoHTA,t="l",col=2)
    plot(30:69,q_HTA/q_NoHTA);points(rr_0913_CI$x,rr_0913_CI$mediana,col=2)

#Sigo con el resto de la tabla lx a Lx y Tx
L11_HTA_spl<-data.frame(x=30:69,L=predict(l11_HTA_spl,30.5:69.5)[[2]])
L11_NoHTA_spl<-data.frame(x=30:69,L=predict(l11_NoHTA_spl,30.5:69.5)[[2]])
L11_spl<-data.frame(x=30:69,L=predict(l11_spl,30.5:69.5)[[2]])
    plot(L11_HTA_spl$x+0.5,L11_HTA_spl$L,col=2);lines(l11_HTA_splp)
T11_HTA_spl<-c();T11_NoHTA_spl<-c();T11_spl<-c()
    for (i in 1:40) {T11_HTA_spl[i]<-sum(L11_HTA_spl$L[i:40])}
    for (i in 1:40) {T11_NoHTA_spl[i]<-sum(L11_NoHTA_spl$L[i:40])}
    for (i in 1:40) {T11_spl[i]<-sum(L11_spl$L[i:40])}

#ev condicionada a tener/no tener/ HTA hasta los 70
ev11_HTA<-data.frame(edad=30:69, ev=T11_HTA_spl/l11_HTA_splp$l[1:40])
ev11_NoHTA<-data.frame(edad=30:69, ev=T11_NoHTA_spl/l11_NoHTA_splp$l[1:40])
ev11<-data.frame(edad=30:69, ev=T11_spl/l11_splp$l[1:40])

###cuadro y gráf mortalidad
df.ev<-data.frame(Edad=30:69,round(q_NoHTA,5), round(q_HTA,5), 
                  round(ev11_NoHTA$ev,1),round(ev11_HTA$ev,1),
                  Diferencia=round(ev11_NoHTA$ev-ev11_HTA$ev,1));
                  View(df.ev)
df.evgr=rbind(data.frame(x=df.ev[,1],ev=df.ev[,2],cond="NoHTA"),
              data.frame(x=df.ev[,1],ev=df.ev[,3],cond="HTA"))

#graf1    
ggplot()+geom_line(data=ev11_NoHTA,aes(x=edad,y=ev),color=1)+
        geom_point(data=ev11_NoHTA,aes(x=edad,y=ev),color=1)+
        geom_line(data=ev11_HTA,aes(x=edad,y=ev),color=2)+
        geom_point(data=ev11_HTA,aes(x=edad,y=ev),color=2)+
        scale_y_continuous(name = "EV",
                           breaks = seq(0,40,5))+
        scale_x_continuous(name = "Edad")+
        ggtitle("Esperanza de vida condicionada a HTA y temporaria por edad.\n Total del País. Años 2009/11")+
        theme(plot.title = element_text(size = 14),
                axis.text=element_text(size=10),
                axis.title=element_text(size=12))

###q(HTA)
#q(HTA) por edad simple
library(splines)
linc<-tabla_l(inc_0913_CI$mediana[seq(1,37,4)])
linc_spl<-interpSpline(seq(30,70,4),linc)
linc_splp<-data.frame(x=30:70,l=predict(linc_spl,30:70)[[2]])
linc_splp$q<-0
for (i in 1:nrow(linc_splp)-1){
    linc_splp$q[i]<-(linc_splp$l[i]-linc_splp$l[i+1])/linc_splp$l[i]}

###tabla conjunta multiestado (con prevalentes al inicio)
#prevalencia inicial a los 30: promedio movil de 5 en ambos años
prev_x30=0.2 
Noinc2<-data.frame(x=30:70,
                  l=c(100000*(1-prev_x30),rep(0,40)),
                  qinc=c(1-0.5*c(q_NoHTA,0))*linc_splp$q,
                  qmort=c((1-0.5*linc_splp$q)*c(q_NoHTA,0)),
                  qmortinc=c(q_HTA,0),
                  dinc=rep(0,41),
                  linc=c(100000*prev_x30,rep(0,40)),
                  dmort=rep(0,41),
                  dincmort=rep(0,41))
for (i in 1:(nrow(Noinc2)-1)){
    Noinc2$dinc[i]<-Noinc2$qinc[i]*Noinc2$l[i]
    Noinc2$linc[i+1]<-Noinc2$qinc[i]*Noinc2$l[i]+Noinc2$linc[i]*(1-q_HTA[i])
    Noinc2$dincmort[i]<-Noinc2$linc[i]*q_HTA[i]
    Noinc2$dmort[i]<-Noinc2$qmort[i]*Noinc2$l[i]
    Noinc2$l[i+1]<-Noinc2$l[i]-Noinc2$dmort[i]-Noinc2$dinc[i]}
Noinc2[,6:9]<-round(Noinc2[,6:9],0)

#Lxs
spNoinc2<-interpSpline(30:70,Noinc2$l)
spNoinc2<-predict(spNoinc2,30.5:69.5)[[2]]
Noinc2$L<-c(spNoinc2,0)
spNoinci2<-interpSpline(30:70,Noinc2$linc)
spNoinci2<-predict(spNoinci2,30.5:69.5)[[2]]
Noinc2$Li<-c(spNoinci2,0)

#Ts
    for (i in 1:40) {Noinc2$TL[i]<-sum(Noinc2$L[i:40])}
    for (i in 1:40) {Noinc2$TLi[i]<-sum(Noinc2$Li[i:40])}
Noinc2$ev11_libre<-c(round(Noinc2$TL[1:40]/Noinc2$l[1:40],1),0)
Noinc2$ev11_NoCond<-c(round((Noinc2$TL[1:40]+Noinc2$TLi[1:40])/
                               (Noinc2$l[1:40]+Noinc2$linc[1:40]),1),0)
#tabla result
TablaEVinc2<-Noinc2[Noinc2$x %in% c(30,40,50,60),]

#gráfico de decrementos
plot(Noinc2$x[1:40],Noinc2$dinc[1:40],type="l",col=2,ylim=c(0,2000), 
    xlab="Edad",ylab=expression(~d[x]), lwd=3,cex.lab=0.7)
    lines(Noinc2$x[1:40],Noinc2$dmort[1:40],col=1,ylim=c(0,2000), 
     xlab="Edad",ylab=expression(~d[x]), lwd=3)
    title("Salidas desde el estado sin HTA. Cohorte hipotética.\nTotal del País. Período 2009/13",cex.main = 0.8)
    legend("topleft", c(expression(~d[x]~~""~a~""~HTA),
                    expression(~d[x]~""~a~""~Muerte)),
                    col = c(2,1),text.col = "black", lty=c(1,1), cex = 0.6, bty = "n")
