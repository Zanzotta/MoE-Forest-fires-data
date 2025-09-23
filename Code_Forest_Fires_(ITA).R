#########################################################################################################
#  Appendice: codice R incendi Montesinho Park (ITA)
#########################################################################################################

# Autori: Boiocchi Federico, Zanzottera Andrea 

rm(list=ls())

# carichiamo le libraries

library(mclust)
library(tidyverse)
library(mixtools)
library(flexmix)
library(GGally)
library(plotly)
library(grid)
library(jpeg)
library(mclust)
library(Rmixmod)
library(caret)


# carichiamo il dataset

forest<-read.csv("C:\\Users\\andre\\Desktop\\HW stat comp\\forestfires.csv") 

forest<-forest[-c(1:8,12)] # selezioniamo le variabili utili (Temp, Wind, RH)
# sono state selezionate le piu interpretabili, tranne Rain (zero inflated)

# istogramma risposta
g1<-ggplot(forest, aes(x=area)) +  
  geom_histogram(aes(y=after_stat(ncount)),fill="green3", color="white", binwidth=10) +       
  labs(title = "Andamento Area", 
       x = "Area bruciata (ettari)", 
       y = "densità di freq.")
g1
# la risposta non ci sembra adatta a una regressione di misture poichè  è troppo asimmetrica
# modifichiamo la risposta con log(area+1) --> risposta diventerà simile a mistura di 3 norm
# le correlazioni sono molto basse, dovremo controllare se nei sottogruppi sono migliorate


cor(forest)
#pairs(forest,pch=19,col=season)

# trasformiamo la risposta
forest[,4]<-log(forest$area+1) 

# rinominiamo le variabili
names<-colnames(forest)
colnames(forest)<-c(names[-4],"log_area")


# centriamo e standardizziamo i dati
forest<-as.data.frame(scale(forest)) 


# istogramma risposta trasformata

g2<-ggplot(forest, aes(x=log_area))+  
  geom_histogram(aes(y=after_stat(ncount)),fill="green3", color="white", binwidth=0.3,alpha=0.5) +       
  geom_density(color = "black" ,
               size = 0.6)+
  labs(title = "Andamento trasformata Area", 
       x = "ln(Area bruciata+1)", 
       y = "densità di freq.") 
g2

anyNA(forest)

pairs<-ggpairs(forest)


##################################################################################################
#    Mixture of Experts
##################################################################################################


set.seed(123) # settiamo il seed per avere risultati sempre uguali

# esplicative : Temp, RH,Wind
# risposta : Log_area

# le covariate entrano sia a modellare le medie di gruppo che le mixing proportions

# stima del il miglior MoE 


fit<-stepFlexmix(log_area~temp+RH+wind,data=forest, k=2:4, 
                 model=FLXMRglm(~temp+RH+wind,family=c("gaussian")),
                 concomitant = FLXPmultinom(~temp+RH+wind),
                 nrep = 5,verbose = TRUE, drop = F, unique = FALSE)

# scegliamo il modello a 3 cluster invece che 2 perchè quello a due non mostra 
# un andamento marcatamente diverso tra cluster dell'effetto di wind e RH su log-area (lo stesso vale per altre coppie di esplicative)
labs<-fit@models$'3'@cluster  # etichette previste (modello a 3 cluster)
labs<-factor(labs)
best3fit<-fit@models$'3'  # miglior modello

# matrice delle posterior probabilities legata alla clusterizzazione fatta con MoE
pp<-best3fit@posterior$scaled
pp<-matrix(pp,nrow=517,ncol=3)


# parametri gamma e sigma (delle 3 regressioni LM in questo caso)

parameters(best3fit)  # questo per mostrare gli effetti diversificati
summary(best3fit)
ICL(best3fit)

# coefficienti della regressione multinomiale logistica (beta)
parameters(best3fit,which="concomitant")


# misure di bontà:
# incertezza media
uncertainty<-matrix(apply(pp,1,function(t) 1-max(t)),ncol=1)
unc<-mean(uncertainty)
Avg_unc<-unc

ind1 <- which(labs == 1)
ind2 <- which(labs == 2)
ind3 <- which(labs == 3)


# correlazioni nei sottogruppi

forest1<-forest[ind1,]
forest2<-forest[ind2,]
forest3<-forest[ind3,]

cor1<-cor(forest1)
cor2<-cor(forest2)
cor3<-cor(forest3) # non completamente calcolabile per via degli zeri

  
# plot 3D con le 2 esplicative RH e Wind (scelte arbitrariamente) 
# e log_area come risposta. Nuovle di punti differenziate per gruppo
# e piani di regressione



# Estende i limiti dei dati per i piani di regressione
x_min <- min(forest$wind)
x_max <- max(forest$wind)
y_min <- min(forest$RH)
y_max <- max(forest$RH)

# Griglie per le previsioni (includendo i limiti)
grid_x1 <- seq(x_min, x_max, length.out = 200)
grid_y1 <- seq(y_min, y_max, length.out = 200)
grid1 <- expand.grid(wind = grid_x1, RH = grid_y1)

grid_x2 <- seq(x_min, x_max, length.out = 200)
grid_y2 <- seq(y_min, y_max, length.out = 200)
grid2 <- expand.grid(wind = grid_x2, RH = grid_y2)

grid_x3 <- seq(x_min, x_max, length.out = 200)
grid_y3 <- seq(y_min, y_max, length.out = 200)
grid3 <- expand.grid(wind = grid_x3, RH = grid_y3)

# Modelli
model1 <- lm(log_area ~ wind + RH, data = forest[ind1, ])
model2 <- lm(log_area ~ wind + RH, data = forest[ind2, ])
model3 <- lm(log_area ~ wind + RH, data = forest[ind3, ])

# Previsioni per le griglie
grid1$log_area <- predict(model1, newdata = grid1)
grid2$log_area <- predict(model2, newdata = grid2)
grid3$log_area <- predict(model3, newdata = grid3)


plot <- plot_ly() %>%
  # Prima nuvola di punti
  add_markers(
    x = forest$wind[ind1], y = forest$RH[ind1], z = forest$log_area[ind1],
    marker = list(color = "red", size = 4, line = list(color = "grey", width = 0.5)),
    name = 'Group 1'
  ) %>%
  # Seconda nuvola di punti
  add_markers(
    x = forest$wind[ind2], y = forest$RH[ind2], z = forest$log_area[ind2],
    marker = list(color = "black", size = 4, line = list(color = "grey", width = 0.5)),
    name = 'Group 2'
  )%>%
  # Terza nuvola di punti
  add_markers(
    x = forest$wind[ind3], y = forest$RH[ind3], z = forest$log_area[ind3],
    marker = list(color = "gold", size = 4, line = list(color = "grey", width = 0.5)),
    name = 'Group 3'
  ) %>%
  # Piano di regressione per il primo gruppo
  add_surface(
    x = grid_x1, y = grid_y1, z = matrix(grid1$log_area, nrow = 200, ncol = 200),
    showscale = FALSE, opacity = 0.4, name = 'Regression Plane 1'
  ) %>%
  # Piano di regressione per il secondo gruppo
  add_surface(
    x = grid_x2, y = grid_y2, z = matrix(grid2$log_area, nrow = 200, ncol = 200),
    showscale = FALSE, opacity = 0.4, name = 'Regression Plane 2'
  ) %>%
  # Piano di regressione per il terzo gruppo
  add_surface(
    x = grid_x3, y = grid_y3, z = matrix(grid3$log_area, nrow = 200, ncol = 200),
    showscale = FALSE, opacity = 0.4, name = 'Regression Plane 3'
  ) %>%
  layout(scene = list(
    xaxis = list(title = 'Wind'),
    yaxis = list(title = 'RH'),
    zaxis = list(title = 'Log Area')
  ))
plot




###################################################################################################################
#  Classification
###################################################################################################################

# estraiamo casualmente il training set (70 percento delle osservazioni) sul quale verrà allenato e selezionato il classifier
# (di fatto estraiamo il suo complementare (test set) e lo mettiamo da parte per la fase di classification)
set.seed(123)
(n<-nrow(forest)) # numero dei pazienti (u.s)
test.set.labels<-sample(1:n,100)   # test set: 155 osservazioni sono circa il 30% dei dati

# Consideriamo le etichette stimate dal FULL MEM 
# come ground truth e facciamo classification sul test set

# fase di learning e validation del classifier

# nel stimare il classifier consideriamo come features le esplicative
# cerchiamo il migliore tra tutti i classifier in termini di stima CV del Mer
# i valori per la cross-validation sono lasciati di default
# ripetiamo la mixmodlearn piu volte e prendiamo il modello migliore

forest<-read.csv("C:\\Users\\andre\\Desktop\\HW stat comp\\forestfires.csv") 
forest<-forest[-c(1:8,12,13)]
forest<-as.data.frame(scale(forest))
#View(forest)
Mod<-list()
mod<-list()
CV<-numeric()
for (i in 1:300){
  res = mixmodLearn(forest[-test.set.labels,], labs[-test.set.labels], 
                    models=mixmodGaussianModel(family="all",equal.proportions=FALSE),
                    criterion=c('CV','BIC'))
  Mod[[i]]<-res
  mod[[i]]<-res@bestResult@model
  CV[i]<-res@bestResult@criterionValue
}


CV_min<-min(CV)
k<-which.min(CV)
mod_opt<-mod[[k]]
Mod_opt<-Mod[[k]]
res<-Mod_opt

BIC = CV = rep(NA ,length(res@models@listModels) )   # faccio un vettore da riempire di lunghezza 14



for (i in 1: length(res@models@listModels)){
  ind = which(res@results [[i]] @model == res@models@listModels) # faccio l'associazione tra risultati e nome del modello. Di fatto vado a cercare chi in results occupa la posizione i-esima e mi chiedo qual è l'indice di quel modello nella lista semplice-complesso
  CV[ind] = res@results [[i]] @criterionValue [1] 
  BIC[ind] = res@results [[i]] @criterionValue [2]
}

# gli andamenti di stima CV del MER
par(mfrow=c(2,1))
plot(BIC ,type='b',pch=19,cex=0.7,xlab='',xaxt='n',col =2,lwd=0.5); axis(1,at=1: length(
  res@results),labels=substr(res@models@listModels ,10 ,30),cex.axis =0.6
  ,las=2)
abline(v=which.min(BIC), col=1, lty =2)

par(mfrow=c(1,1))
plot(CV ,type='b',pch=19,cex=0.9,xlab='',xaxt='n',col =3); axis(1,at=1: length(
  res@results),labels=substr(res@models@listModels,10 ,30),cex.axis =0.6
  ,las =2)
abline(v=which.min(CV), col=1, lty =2)
par(mfrow=c(1,1))
min(CV)

# fase di classfication vera e propria sul test set

PREDICTION<- mixmodPredict(data=forest[test.set.labels,], 
                           classificationRule=Mod_opt@bestResult)

PREDICTION@partition  # etichette assegnate dal classifier
PREDICTION@proba # matrice delle posterior probabilities


# vera confusion matrix da guardare

pred_lab <- as.factor(PREDICTION@partition)
ground_lab <-as.factor(labs[test.set.labels])
levels(pred_lab)<-c("1","2","3")
str(pred_lab)
str(ground_lab)
confusionMatrix(data=pred_lab,reference=ground_lab)



#################################################################################################
# grafici

forest<-read.csv("C:\\Users\\andre\\Desktop\\HW stat comp\\forestfires.csv") 

graf_day<-ggplot(forest, aes(x=day,          
                             y= log(area+1),      
                             fill=labs))+
  scale_fill_manual(values=c("red","black","gold"))+# sottogruppi
  geom_boxplot(position=position_dodge(preserve="single"),alpha=0.6)

library(ggridges)
ggplot(forest, aes(x=log(area+1),     
                   y=month,     # mettere anche day  
                   fill=month,alpha=0.7)) + # metter anche day
  geom_density_ridges(show.legend = F)+
  theme_ridges()

forest<-forest[-c(1:8,12)] 
forest[,4]<-log(forest$area+1)


# ggpairs con densità non parametriche
graf_pairs <- ggpairs(forest,
                      columns = 1:4,
                      mapping = aes(alpha = 0.6),  # Impostazione dell'alpha per la trasparenza
                      lower = list(continuous = wrap("points", size = 0.5)),  # Ridurre la grandezza dei punti negli scatter plots
                      diag = list(continuous = function(data, mapping, ...) {
                        ggplot(data = data, mapping = mapping) +
                          geom_density(aes(y = ..density..), fill = "gray", color = "black", alpha = 0.5, ...)  # Densità riempita di grigio
                      })) +
theme_minimal()


# codice per la mappa

mappa <- readJPEG("C:\\Users\\andre\\Desktop\\HW stat comp\\Immagine 2025-01-25 152208.jpg")

dim_mappa <- dim(mappa) # Ottiene altezza, larghezza e numero di canali (RGB)

# Dimensioni della mappa in pixel
altezza <- dim_mappa[1]
larghezza <- dim_mappa[2]

#creo una griglia 9x9
n_col <- 9
n_row <- 9

# Calcolare i limiti della griglia in base alle dimensioni dell'immagine
x_grid <- seq(0, larghezza, length.out = n_col + 1)
y_grid <- seq(0, altezza, length.out = n_row + 1)


# recupero dal dataset iniziale le variabili X e Y
ff <- read.csv("C:\\Users\\andre\\Desktop\\HW stat comp\\forestfires.csv") 

# traduco le variabili X e Y del dataset iniziale
# in coordinate in pixel disegnabili sull'immagine

coord_pixel_x <- numeric()
coord_pixel_y <-numeric()
for (i in 1:517){
  coord_pixel_x[i] <- x_grid[ff$X[i]+1]
  coord_pixel_y[i] <- altezza-y_grid[ff$Y[i]]
}

# aggiungo del rumore se no avrei osservazioni tutte sovrapposte

coord_pixel_x_jitter <-jitter(coord_pixel_x,amount=50)
coord_pixel_y_jitter <-jitter(coord_pixel_y,amount=20)

coord<-data.frame(cbind(coord_pixel_x_jitter,coord_pixel_y_jitter))
str(coord)

gruppi <- factor(labs)

# Creare il plot con la mappa, la griglia e i punti
map<-ggplot() +
  # per  Visualizzare l'immagine come sfondo
  annotation_custom(
    rasterGrob(mappa, 
               width = unit(1, "npc"), 
               height = unit(1, "npc"),
               interpolate = TRUE),
    xmin = 0, xmax = larghezza,
    ymin = 0, ymax = altezza
  ) +
  # Aggiungere la griglia
  geom_vline(xintercept = x_grid, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = y_grid, color = "grey", linetype = "dashed") +
  # Aggiungere i punti colorati in base alle etichette
  geom_point(data = coord, aes(x = coord_pixel_x_jitter, y =coord_pixel_y_jitter , fill =gruppi), 
             color = "black",                 # Bordo nero
             size = 2,                        # Dimensione dei punti
             shape = 21) +                    # Cerchio pieno con bordo
  scale_fill_manual(values = c( "red","black","gold")) + # Colori per etichetta   
  # Adattare il rapporto della figura
  coord_fixed(ratio = 1, xlim = c(0, larghezza), ylim = c(0, altezza)) +
  # Personalizzo gli assi
  scale_x_continuous(breaks = seq(0, larghezza, by = larghezza / n_col),
                     labels = 0:9) +
  scale_y_continuous(breaks = seq(0, altezza, by = altezza / n_row),
                     labels = 9:0) +
  labs(
    title = "Incendi Montesinho park",
    x = "X jittered (Griglia 0-9)",
    y = "Y jittered (Griglia 0-9)"
  ) +
  theme_minimal()





