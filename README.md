# code-B121


### PACKAGES UTILISES

library("vegan")
library("ade4")
library("adegraphics")
library("berryFunctions")
library("dplyr")
library("tidyverse")
library("ggpubr")
library("rstatix")
library("RColorBrewer")
library("betapart")
library("ggpubr")
library("ggsci")
library("ggplot2")
library("factoextra")
library("FactoMineR")
library("ape")
library("iNEXT")
library("hclust")




### IMPORTATION DES JEUX DE DONNEES ###

load("Biigletot.RData")
load("sp2_8.RData")
load("Fusion2_8.RData")
load("MI_8.RData")
load("BI_8.RData")
load("FH_8.RData")
load("GR_8.RData")
load("HI_8.RData")
load("NH_8.RData")
load("SK_8.RData")
load("UI_8.RData")
load("Faune.RData")



### MISE EN FORME DU JEU DE DONNEES ###
   ## Changer les NA du jeu de données en 0:
Biigletot[is.na(Biigletot)] = 0
  
   ## Enlever les taxons non identifiés:
sp2_8<-sp2_8[,-16]
sp2_8<-sp2_8[,-12]
sp2_8<-sp2_8[,-21]

## Regrouper les données par site:
LeaGroup<-spe2_8 %>% group_by(station) %>% summarise(across(everything(),sum,na.rm=TRUE))


### CALCUL RICHESSE SPECIFIQUE PAR STATION ###
  ## Richesse spécifique (nombre d'espèces par station) grâce à la fonction specnumber:
Richessespe_8<-specnumber(sp2_8[,-1])

  ## Plot de la richesse spécifique:
barplot(Richessespe_8, names.arg = sp2_8$station, col=terrain.colors(8),ylab = "Nombre d'espèces par station", xlab="Stations", 
        main="Richesse spécifique des 8 stations", ylim=c(0,20))



### COURBES D'ACCUMULATION DES ESPECES PAR STATION ###
  ## Séparation des données par station grâce à la fonction subset:
MI<-subset(Biigletot, station == "MI")
BI<-subset(Biigletot, station=="BI")
FH<-subset(Biigletot, station=="FH")
GR<-subset(Biigletot, station=="GR")
HI<-subset(Biigletot, station=="HI")
NH<-subset(Biigletot, station=="NH")
SK<-subset(Biigletot, station=="SK")
UI<-subset(Biigletot, station=="UI")

  ## Courbe d'accumulation des espèces par station avec la fonction specaccum:
AccMI_8<-specaccum(MI_8[,-1])
AccBI_8<-specaccum(BI_8[,-1])
AccFH_8<-specaccum(FH_8[,-1])
AccGR_8<-specaccum(GR_8[,-1])
AccHI_8<-specaccum(HI_8[,-1])
AccNH_8<-specaccum(NH_8[,-1])
AccSK_8<-specaccum(SK_8[,-1])
AccUI_8<-specaccum(UI_8[,-1])


  ## Plot des courbes d'accumulation
par(mfcol=c(4,2)) #permet de disposer tous les plots sur une même page divisée en 4 lignes et 2 colonnes

plotMI_8<-plot(AccMI, col="blue", ci.type="poly", ci.col="lightblue", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
     le site MI", ylim=c(0,23))

plotBI_8<-plot(AccBI, col="blue", ci.type="poly", ci.col="yellow", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site BI", ylim=c(0,17))

plotFH<-plot(AccFH, col="blue", ci.type="poly", ci.col="purple", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site FH", ylim=c(0,17))

plotGR<-plot(AccGR, col="blue", ci.type="poly", ci.col="red", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site GR", ylim=c(0,18))

plotHI<-plot(AccHI, col="blue", ci.type="poly", ci.col="green", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site HI",ylim=c(0,17))

plotNH<-plot(AccNH, col="blue", ci.type="poly", ci.col="orange", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site NH", ylim=c(0,17))

plotSK<-plot(AccSK, col="blue", ci.type="poly", ci.col="brown", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site SK", ylim=c(0,19))

plotUI<-plot(AccUI, col="blue", ci.type="poly", ci.col="darkgreen", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site UI", ylim=c(0,15))


### CALCUL DES ESTIMATEURS (CHAO, JACKKNIFE 1 ET 2,  BOOTSTRAP) avec les fonctions specpool et poolaccum ###

specpool(MI[,-1])
poolMI<-poolaccum(MI[,-1])
plot(poolMI, main="Calcul de diversité de la station MI par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(BI[,-1])
poolBI<-poolaccum(BI[,-1])
plot(poolBI, main="Calcul de diversité de la station BI par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(FH[,-1])
poolFH<-poolaccum(FH[,-1])
plot(poolFH, main="Calcul de diversité de la station FH par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(GR[,-1])
poolGR<-poolaccum(GR[,-1])
plot(poolGR, main="Calcul de diversité de la station GR par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(HI[,-1])
poolHI<-poolaccum(HI[,-1])
plot(poolHI, main="Calcul de diversité de la station HI par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(NH[,-1])
poolNH<-poolaccum(NH[,-1])
plot(poolNH, main="Calcul de diversité de la station NH par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(SK[,-1])
poolSK<-poolaccum(SK[,-1])
plot(poolSK, main="Calcul de diversité de la station SK par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")

specpool(UI[,-1])
poolUI<-poolaccum(UI[,-1])
plot(poolUI, main="Calcul de diversité de la station UI par extrapolation 
     selon 4 méthodes (Chao, Jackknife 1, Jackknife 2, Boostrap)")




### CALCUL DES INDICES DE DIVERSITE ###
  ## Tri des données pour ne garder que la faune hors porifères:
Faune<-sp2_8[, -(2:8)]
Faune<-Faune[, -3]
Faune<-Faune[, -(20:23)]

  ## Shannon, fonction diversity: 
Shannon<-diversity(Faune[-1])
barplot(Shannon, names.arg = Faune$station, col=brewer.pal(n = 8, name = "Blues"),ylab = "Indice de Shannon", xlab="Stations", 
        main="Indice de Shannon des 8 stations étudiées", ylim=c(0,2))
abline(0,0) # permet d'ajouter une ligne sur le graphique à l'origine

  ## Simpson :
diversity(Faune[-1], index = "simpson" )
Simpson<-diversity(Faune[-1], index = "simpson" )
barplot(Simpson, names.arg = Faune$station, col=brewer.pal(n = 8, name = "Purples"),ylab = "Indice de Simpson", xlab="Stations", 
        main="Indice de Simpson des 8 stations étudiées", ylim=c(0,1))
abline(0,0)

  ## Piélou:
S<-specnumber(Faune[,-1])
J<-Shannon/log(S)
barplot(J, names.arg = sp2_8$station, col=brewer.pal(n = 8, name = "YlGn"),ylab = "Indice de Pielou", xlab="Stations", 
        main="Indice de Pielou des 8 stations étudiées", ylim=c(0,1))
abline(0,0)

  

  ## Indice de Bray-Curtis, fonction vegdist:
PrasAbs<-ifelse(Spe2_8>0,1,0)  # Permet de changer les données en présence-absence

Bray<-vegdist(PrasAbs[,-1], method = "bray", binary = TRUE, diag=TRUE)

     # Clustering des données de l'indice de Bray-Curtis :
Clust<-hclust(Bray) 
DendroBray<-as.dendrogram(Clust)
plot(DendroBray)




###ABONDANCE RELATIVE DES ESPECES ANIMALES HORS PORIFERES
load(file="Pourcentage_sp.RData")

PlotPourcen<-ggplot(Pourcentage_sp, aes(Station, Taux, fill=Taxons), na.rm=TRUE)+
  geom_bar(stat="identity", alpha=0.9)+
  scale_y_continuous(limits=c(0, 100))+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  labs(x="Stations",
       y="Abondance proportionnelle des taxons")

PlotPourcen




###RECOUVREMENT PAR LES ALGUES###
setwd(dir="C:/Users/Utilisateur/Desktop/Memoire/Biigle/Recouvrement")
load(file = "TotSpVeg.RData")


	## Division des données en fichier par espèce par site:
HimantoBIJ<-subset(BI,Espece=="Himantothallus grandifolius")
MonomastraBIJ<-subset(BI, Espece == "Monostroma hariotii")
DesmasteriaBIJ<-subset(BI, Espece=="Desmarestia")
IridiaBIJ<-subset(BI, Espece=="Iridaea cordata")
PlacomiumBIJ<-subset(BI, Espece=="Plocamium hookeri")
TrematocarpusBIJ<-subset(BI,Espece=="Trematocarpus antarcticus")
LeptoBIJ<-subset(BI,Espece=="Leptophytum coulmanicum")
UnidBIJ<-subset(BI, Espece=="unidalgae")


	##BI:
#49 photos donc 49 * 8 294 192 = 406 415 408  pixels

sommeiridiaBIJ<-colSums(IridiaBIJ [,-(1:2)])
PourIridiaBIJ<-c(sommeiridiaBIJ*100/406415408)
sommetrematoBIJ<-colSums(TrematocarpusBIJ[,-(1:2)])
PourTrematoBIJ<-c(sommetrematoBIJ*100/406415408)
sommeleptoBIJ<-colSums(LeptoBIJ[,-(1:2)])
PourLeptoIJ<-c(sommeleptoBIJ*100/406415408)
SUIJ<-colSums(UnidBIJ[,-(1:2)])
PUIJ<-c(SUIJ*100/406415408)
SUIJ
PUIJ



	## FH:
# 44 photos analysées donc 44* 8 294 192 = 364944448 pixels

HimantoFHJ<-subset(FH,Espece=="Himantothallus grandifolius")
MonomastraFHJ<-subset(FH, Espece == "Monostroma hariotii")
DesmasteriaFHJ<-subset(FH, Espece=="Desmarestia")
IridiaFHJ<-subset(FH, Espece=="Iridaea cordata")
PlacomiumFHJ<-subset(FH, Espece=="Plocamium hookeri")
TrematocarpusFHJ<-subset(FH,Espece=="Trematocarpus antarcticus")
LeptoFHJ<-subset(FH,Espece=="Leptophytum coulmanicum")
UnidFHJ<-subset(FH, Espece=="unidalgae")

sommeHimantoFHJ<-colSums(HimantoFHJ [,-(1:2)])
PourHimantoFHJ<-c(sommeHimantoFHJ*100/364944448)
sommedesmasFHJ<-colSums(DesmasteriaFHJ[,-(1:2)])
PourdesmasJ<-c(sommedesmasFHJ*100/364944448)
sommeDesmasteriaFHJ<-colSums(DesmasteriaFHJ[,-(1:2)])
PourDesmasFHJ<-c(sommeDesmasteriaFHJ*100/364944448)
sommePlacoFHJ<-colSums(PlacomiumFHJ[,-(1:2)])
PourPlacoFHJ<-c(sommePlacoFHJ*100/364944448)
sommeTremaFHJ<-colSums(TrematocarpusFHJ[,-(1:2)])
PourTremaFHJ<-c(sommeTremaFHJ*100/364944448)
sommeLeptoFHJ<-colSums(LeptoFHJ[,-(1:2)])
PourLeptoFHJ<-c(sommeLeptoFHJ*100/364944448)
sommeIriFHJ<-colSums(IridiaFHJ[,-(1:2)])
PourIriFHJ<-c(sommeIriFHJ*100/364944448)
sommeMonoFHJ<-colSums(MonomastraFHJ[,-(1:2)])
PourMonoFHJ<-c(sommeMonoFHJ*100/364944448)
SUIJ<-colSums(UnidFHJ[,-(1:2)])
PUIJ<-c(SUIJ*100/364944448)
SUIJ
PUIJ


	## GR
# 79 photos analysées donc 79* 8 294 192 = 655241168 pixels

Hi_GRJ<-subset(GR,Espece=="Himantothallus grandifolius")
M_GRJ<-subset(GR, Espece == "Monostroma hariotii")
D_GRJ<-subset(GR, Espece=="Desmarestia")
I_GRJ<-subset(GR, Espece=="Iridaea cordata")
P_GRJ<-subset(GR, Espece=="Plocamium hookeri")
T_GRJ<-subset(GR,Espece=="Trematocarpus antarcticus")
L_GRJ<-subset(GR,Espece=="Leptophytum coulmanicum")
A_GRJ<-subset(GR, Espece=="Ascoseira mirabilis")
UnidGRJ<-subset(GR, Espece=="unidalgae")


S__Hi_GRJ<-colSums(Hi_GRJ [,-(1:2)])
P_Hi_GRJ<-c(S__Hi_GRJ*100/655241168)
S__Hi_GRJ
P_Hi_GRJ
S_M_GRJ<-colSums(D_GRJ[,-(1:2)])
P_M_GRJ<-c(S_M_GRJ*100/655241168)
S_M_GRJ
S_D_GRJ<-colSums(D_GRJ[,-(1:2)])
P_D_GRJ<-c(S_D_GRJ*100/655241168)
S_D_GRJ
P_D_GRJ
S_P_GRJ<-colSums(P_GRJ[,-(1:2)])
P_P_GRJ<-c(S_P_GRJ*100/655241168)
S_P_GRJ
P_P_GRJ
S_T_GRJ<-colSums(T_GRJ[,-(1:2)])
P_T_GRJ<-c(S_T_GRJ*100/655241168)
S_T_GRJ
P_T_GRJ
S_L_GRJ<-colSums(L_GRJ[,-(1:2)])
P_L_GRJ<-c(S_L_GRJ*100/655241168)
S_L_GRJ
S_I_GRJ<-colSums(I_GRJ[,-(1:2)])
P_I_GRJ<-c(S_I_GRJ*100/655241168)
S_I_GRJ
P_I_GRJ
S_A_GRJ<-colSums(A_GRJ [,-(1:2)])
P_A_GRJ<-c(S_A_GRJ*100/655241168)
S_A_GRJ
P_A_GRJ
SUIJ<-colSums(UnidGRJ[,-(1:2)])
PUIJ<-c(SUIJ*100/655241168)
SUIJ
PUIJ

	## HI
# 88 photos analysées donc 88* 8 294 192 = 729888896 pixels

Hi_HIJ<-subset(HI,Espece=="Himantothallus grandifolius")
M_HIJ<-subset(HI, Espece == "Monostroma hariotii")
D_HIJ<-subset(HI, Espece=="Desmarestia")
I_HIJ<-subset(HI, Espece=="Iridaea cordata")
P_HIJ<-subset(HI, Espece=="Plocamium hookeri")
T_HIJ<-subset(HI,Espece=="Trematocarpus antarcticus")
L_HIJ<-subset(HI,Espece=="Leptophytum coulmanicum")
UnidHIJ<-subset(HI, Espece=="unidalgae")

S_Hi_HIJ<-colSums(Hi_HIJ [,-(1:2)])
P_Hi_HIJ<-c(S_Hi_HIJ*100/729888896)
S__Hi_HIJ
P_Hi_HIJ
S_M_HIJ<-colSums(D_HIJ[,-(1:2)])
P_M_HIJ<-c(S_M_HIJ*100/729888896)
S_M_HIJ
S_D_HIJ<-colSums(D_HIJ[,-(1:2)])
P_D_HIJ<-c(S_D_HIJ*100/729888896)
S_D_HIJ
P_D_HIJ
S_I_HIJ<-colSums(I_HIJ[,-(1:2)])
P_I_HIJ<-c(S_I_HIJ*100/729888896)
S_I_HIJ
P_I_HIJ
S_P_HIJ<-colSums(P_HIJ[,-(1:2)])
P_P_HIJ<-c(S_P_HIJ*100/729888896)
S_P_HIJ
P_P_HIJ
S_T_HIJ<-colSums(T_HIJ[,-(1:2)])
P_T_HIJ<-c(S_T_HIJ*100/729888896)
S_T_HIJ
P_T_HIJ
S_L_HIJ<-colSums(L_HIJ[,-(1:2)])
P_L_HIJ<-c(S_L_HIJ*100/729888896)
S_L_HIJ
P_L_HIJ
SUIJ<-colSums(UnidHIJ[,-(1:2)])
PUIJ<-c(SUIJ*100/729888896)
SUIJ
PUIJ


	## MI
# 192 photos analysées donc 192 * 8 294 192 = 1592484860 pixels

Hi_MIJ<-subset(MI,Espece=="Himantothallus grandifolius")
M_MIJ<-subset(MI, Espece == "Monostroma hariotii")
D_MIJ<-subset(MI, Espece=="Desmarestia")
I_MIJ<-subset(MI, Espece=="Iridaea cordata")
P_MIJ<-subset(MI, Espece=="Plocamium hookeri")
T_MIJ<-subset(MI,Espece=="Trematocarpus antarcticus")
L_MIJ<-subset(MI,Espece=="Leptophytum coulmanicum")
UnidMIJ<-subset(MI, Espece=="unidalgae")

S_Hi_MIJ<-colSums(Hi_MIJ [,-(1:2)])
P_Hi_MIJ<-c(S_Hi_MIJ*100/1592484860)
S_Hi_MIJ
P_Hi_MIJ
S_M_MIJ<-colSums(D_MIJ[,-(1:2)])
P_M_MIJ<-c(S_M_MIJ*100/1592484860)
S_M_MIJ
P_M_MIJ
S_D_MIJ<-colSums(D_MIJ[,-(1:2)])
P_D_MIJ<-c(S_D_MIJ*100/1592484860)
S_D_MIJ
P_D_MIJ
S_I_MIJ<-colSums(I_MIJ[,-(1:2)])
P_I_MIJ<-c(S_I_MIJ*100/1592484860)
S_I_MIJ
P_I_MIJ
S_P_MIJ<-colSums(P_MIJ[,-(1:2)])
P_P_MIJ<-c(S_P_MIJ*100/1592484860)
S_P_MIJ
P_P_MIJ
S_T_MIJ<-colSums(T_MIJ[,-(1:2)])
P_T_MIJ<-c(S_T_MIJ*100/1592484860)
S_T_MIJ
P_T_MIJ
S_L_MIJ<-colSums(L_MIJ[,-(1:2)])
P_L_MIJ<-c(S_L_MIJ*100/1592484860)
S_L_MIJ
P_L_MIJ
SUIJ<-colSums(UnidMIJ[,-(1:2)])
PUIJ<-c(SUIJ*100/592484860)
SUIJ
PUIJ


	## NH
# 166 photos analysées donc 166* 8 294 192 = 1376835870 pixels

Hi_NHJ<-subset(NH,Espece=="Himantothallus grandifolius")
M_NHJ<-subset(NH, Espece == "Monostroma hariotii")
D_NHJ<-subset(NH, Espece=="Desmarestia")
I_NHJ<-subset(NH, Espece=="Iridaea cordata")
P_NHJ<-subset(NH, Espece=="Plocamium hookeri")
T_NHJ<-subset(NH,Espece=="Trematocarpus antarcticus")
L_NHJ<-subset(NH,Espece=="Leptophytum coulmanicum")
UnidNHJ<-subset(NH, Espece=="unidalgae")

S_Hi_NHJ<-colSums(Hi_NHJ [,-(1:2)])
P_Hi_NHJ<-c(S_Hi_NHJ*100/1376835870)
S_Hi_NHJ
P_Hi_NHJ
S_M_NHJ<-colSums(D_NHJ[,-(1:2)])
P_M_NHJ<-c(S_M_NHJ*100/1376835870)
S_M_NHJ
P_M_NHJ
S_D_NHJ<-colSums(D_NHJ[,-(1:2)])
P_D_NHJ<-c(S_D_NHJ*100/1376835870)
S_D_NHJ
P_D_NHJ
S_I_NHJ<-colSums(I_NHJ[,-(1:2)])
P_I_NHJ<-c(S_I_NHJ*100/1376835870)
S_I_NHJ
P_I_NHJ
S_P_NHJ<-colSums(P_NHJ[,-(1:2)])
P_P_NHJ<-c(S_P_NHJ*100/1376835870)
S_P_NHJ
P_P_NHJ
S_T_NHJ<-colSums(T_NHJ[,-(1:2)])
P_T_NHJ<-c(S_T_NHJ*100/1376835870)
S_T_NHJ
P_T_NHJ
S_L_NHJ<-colSums(L_NHJ[,-(1:2)])
P_L_NHJ<-c(S_L_NHJ*100/1376835870)
S_L_NHJ
P_L_NHJ
SUIJ<-colSums(UnidNHJ[,-(1:2)])
PUIJ<-c(SUIJ*100/1376835870)
SUIJ
PUIJ



	## SK
# 133 photos analysées donc 133* 8 294 192 = 1103127540 pixels

Hi_SKJ<-subset(SK,Espece=="Himantothallus grandifolius")
M_SKJ<-subset(SK, Espece == "Monostroma hariotii")
D_SKJ<-subset(SK, Espece=="Desmarestia")
I_SKJ<-subset(SK, Espece=="Iridaea cordata")
P_SKJ<-subset(SK, Espece=="Plocamium hookeri")
T_SKJ<-subset(SK,Espece=="Trematocarpus antarcticus")
L_SKJ<-subset(SK,Espece=="Leptophytum coulmanicum")
UnidSKJ<-subset(SK, Espece=="unidalgae")

S_Hi_SKJ<-colSums(Hi_SKJ [,-(1:2)])
P_Hi_SKJ<-c(S_Hi_SKJ*100/1103127540)
S_Hi_SKJ
P_Hi_SKJ
S_M_SKJ<-colSums(D_SKJ[,-(1:2)])
P_M_SKJ<-c(S_M_SKJ*100/1103127540)
S_M_SKJ
P_M_SKJ
S_D_SKJ<-colSums(D_SKJ[,-(1:2)])
P_D_SKJ<-c(S_D_SKJ*100/1103127540)
S_D_SKJ
P_D_SKJ
S_I_SKJ<-colSums(I_SKJ[,-(1:2)])
P_I_SKJ<-c(S_I_SKJ*100/1103127540)
S_I_SKJ
P_I_SKJ
S_P_SKJ<-colSums(P_SKJ[,-(1:2)])
P_P_SKJ<-c(S_P_SKJ*100/1103127540)
S_P_SKJ
P_P_SKJ
S_T_SKJ<-colSums(T_SKJ[,-(1:2)])
P_T_SKJ<-c(S_T_SKJ*100/1103127540)
S_T_SKJ
P_T_SKJ
S_L_SKJ<-colSums(L_SKJ[,-(1:2)])
P_L_SKJ<-c(S_L_SKJ*100/1103127540)
S_L_SKJ
P_L_SKJ
SUIJ<-colSums(UnidSKJ[,-(1:2)])
PUIJ<-c(SUIJ*100/1103127540)
SUIJ
PUIJ


	## UI
# 103 photos analysées donc 103* 8 294 192 = 854301776 pixels

Hi_UIJ<-subset(UI,Espece=="Himantothallus grandifolius")
M_UIJ<-subset(UI, Espece == "Monostroma hariotii")
D_UIJ<-subset(UI, Espece=="Desmarestia")
I_UIJ<-subset(UI, Espece=="Iridaea cordata")
P_UIJ<-subset(UI, Espece=="Plocamium hookeri")
T_UIJ<-subset(UI,Espece=="Trematocarpus antarcticus")
L_UIJ<-subset(UI,Espece=="Leptophytum coulmanicum")
UnidUIJ<-subset(UI, Espece=="unidalgae")

S_Hi_UIJ<-colSums(Hi_UIJ [,-(1:2)])
P_Hi_UIJ<-c(S_Hi_UIJ*100/854301776)
S_M_UIJ<-colSums(D_UIJ[,-(1:2)])
P_M_UIJ<-c(S_M_UIJ*100/854301776)
S_D_UIJ<-colSums(D_UIJ[,-(1:2)])
P_D_UIJ<-c(S_D_UIJ*100/854301776)
S_I_UIJ<-colSums(I_UIJ[,-(1:2)])
P_I_UIJ<-c(S_I_UIJ*100/854301776)
S_P_UIJ<-colSums(P_UIJ[,-(1:2)])
P_P_UIJ<-c(S_P_UIJ*100/854301776)
S_T_UIJ<-colSums(T_UIJ[,-(1:2)])
P_T_UIJ<-c(S_T_UIJ*100/854301776)
S_L_UIJ<-colSums(L_UIJ[,-(1:2)])
P_L_UIJ<-c(S_L_UIJ*100/854301776)
SUIJ<-colSums(UnidUIJ[,-(1:2)])
PUIJ<-c(SUIJ*100/854301776)


	##Plot du taux de recouvrement des algues
setwd(dir="C:/Users/Utilisateur/Desktop/Memoire/Biigle/Recouvrement")
load("Pourcentage.RData")

Pourcentage[is.na(Pourcentage)] = 0
  
PlotRecouv<-ggplot(Pourcentage, aes(Station, Taux, fill=Espece), na.rm=TRUE)+
  geom_bar(stat="identity", alpha=0.9)+
  scale_y_continuous(limits=c(0, 100))+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  labs(x="Stations",
       y="Taux de recouvrement (en %)")

PlotRecouv






###TAUX DE RECOUVREMENT PAR LES PORIFERES

	##BI:
#Division des fichiers par espèce et par site:
DenBIJ<-subset(BI,Espece=="Dendrilla antarctica")
MyBIJ<-subset(BI, Espece == "Mycale acerata")
KiBIJ<-subset(BI, Espece=="Kirkpatrickia variolosa")
HoBIJ<-subset(BI, Espece=="Homaxinella balfourensis")

#49 photos donc 49 * 8 294 192 = 406 415 408  pixels en tout

sommeDenBIJ<-colSums(DenBIJ[,-(1:2)])
PourDenBIJ<-c(sommeDenBIJ*100/406415408)
sommeMyBIJ<-colSums(MyBIJ[,-(1:2)])
PourMyBIJ<-c(sommeMyBIJ*100/406415408)
sommeKiBIJ<-colSums(KiBIJ[,-(1:2)])
PourKiIJ<-c(sommeKiBIJ*100/406415408)
sommeHoBIJ<-colSums(HoBIJ[,-(1:2)])
PourHoIJ<-c(sommeHoBIJ*100/406415408)


	## FH:
# 44 photos analysées donc 44* 8 294 192 = 364944448 pixels

DenFHJ<-subset(FH,Espece=="Dendrilla antarctica")
MyFHJ<-subset(FH, Espece == "Mycale acerata")
KiFHJ<-subset(FH, Espece=="Kirkpatrickia variolosa")
HoFHJ<-subset(FH, Espece=="Homaxinella balfourensis")


sommeDenFHJ<-colSums(DenFHJ [,-(1:2)])
PourDenFHJ<-c(sommeDenFHJ*100/364944448)
sommeMyFHJ<-colSums(MyFHJ[,-(1:2)])
PourMyFHJ<-c(sommeMyFHJ*100/364944448)
sommeKiFHJ<-colSums(KiFHJ[,-(1:2)])
PourKiFHJ<-c(sommeKiFHJ*100/364944448)
sommeHoFHJ<-colSums(HoFHJ[,-(1:2)])
PourHoFHJ<-c(sommeHoFHJ*100/364944448)


	## GR:
# 79 photos analysées donc 79* 8 294 192 = 655241168 pixels

DenGRJ<-subset(GR,Espece=="Dendrilla antarctica")
MyGRJ<-subset(GR, Espece == "Mycale acerata")
KiGRJ<-subset(GR, Espece=="Kirkpatrickia variolosa")
HoGRJ<-subset(GR, Espece=="Homaxinella balfourensis")

S_DenGRJ<-colSums(DenGRJ [,-(1:2)])
P_DenGRJ<-c(S_DenGRJ*100/655241168)
S_MyGRJ<-colSums(MyGRJ[,-(1:2)])
P_MyGRJ<-c(S_MyGRJ*100/655241168)
S_KiGRJ<-colSums(KiGRJ[,-(1:2)])
P_KiGRJ<-c(S_KiGRJ*100/655241168)
S_HoGRJ<-colSums(HoGRJ[,-(1:2)])
P_HoGRJ<-c(S_HoGRJ*100/655241168)


	## HI:
# 88 photos analysées donc 88* 8 294 192 = 729888896 pixels

DenHIJ<-subset(HI,Espece=="Dendrilla antarctica")
MyHIJ<-subset(HI, Espece == "Mycale acerata")
KiHIJ<-subset(HI, Espece=="Kirkpatrickia variolosa")
HoHIJ<-subset(HI, Espece=="Homaxinella balfourensis")

S_DenHIJ<-colSums(DenHIJ [,-(1:2)])
P_DenHIJ<-c(S_DenHIJ*100/729888896)
S_MyHIJ<-colSums(MyHIJ[,-(1:2)])
P_MyHIJ<-c(S_MyHIJ*100/729888896)
S_KiHIJ<-colSums(KiHIJ[,-(1:2)])
P_KiHIJ<-c(S_KiHIJ*100/729888896)
S_HoHIJ<-colSums(HoHIJ[,-(1:2)])
P_HoHIJ<-c(S_HoHIJ*100/729888896)



	## MI:
# 192 photos analysées donc 192 * 8 294 192 = 1592484860 pixels

DenMIJ<-subset(MI,Espece=="Dendrilla antarctica")
MyMIJ<-subset(MI, Espece == "Mycale acerata")
KiMIJ<-subset(MI, Espece=="Kirkpatrickia variolosa")
HoMIJ<-subset(MI, Espece=="Homaxinella balfourensis")


S_DenMIJ<-colSums(DenMIJ [,-(1:2)])
P_DenMIJ<-c(S_DenMIJ*100/1592484860)
S_MyMIJ<-colSums(MyMIJ[,-(1:2)])
P_MyMIJ<-c(S_MyMIJ*100/1592484860)
S_KiMIJ<-colSums(KiMIJ[,-(1:2)])
P_KiMIJ<-c(S_KiMIJ*100/1592484860)
S_HoMIJ<-colSums(HoMIJ[,-(1:2)])
P_HoMIJ<-c(S_HoMIJ*100/1592484860)


	## NH:
# 166 photos analysées donc 166* 8 294 192 = 1376835870 pixels

DenNHJ<-subset(NH,Espece=="Dendrilla antarctica")
MyNHJ<-subset(NH, Espece == "Mycale acerata")
KiNHJ<-subset(NH, Espece=="Kirkpatrickia variolosa")
HoNHJ<-subset(NH, Espece=="Homaxinella balfourensis")


S_DenNHJ<-colSums(DenNHJ [,-(1:2)])
P_DenNHJ<-c(S_DenNHJ*100/1376835870)
S_MyNHJ<-colSums(MyNHJ[,-(1:2)])
P_MyNHJ<-c(S_MyNHJ*100/1376835870)
S_KiNHJ<-colSums(KiNHJ[,-(1:2)])
P_KiNHJ<-c(S_KiNHJ*100/1376835870)
S_HoNHJ<-colSums(HoNHJ[,-(1:2)])
P_HoNHJ<-c(S_HoNHJ*100/1376835870)


	## SK:
# 133 photos analysées donc 133* 8 294 192 = 1103127540 pixels

DenSKJ<-subset(SK,Espece=="Dendrilla antarctica")
MySKJ<-subset(SK, Espece == "Mycale acerata")
KiSKJ<-subset(SK, Espece=="Kirkpatrickia variolosa")
HoSKJ<-subset(SK, Espece=="Homaxinella balfourensis")


S_DenSKJ<-colSums(DenSKJ [,-(1:2)])
P_DenSKJ<-c(S_DenSKJ*100/1103127540)
S_MySKJ<-colSums(MySKJ[,-(1:2)])
P_MySKJ<-c(S_MySKJ*100/1103127540)
S_KiSKJ<-colSums(KiSKJ[,-(1:2)])
P_KiSKJ<-c(S_KiSKJ*100/1103127540)
S_HoSKJ<-colSums(HoSKJ[,-(1:2)])
P_HoSKJ<-c(S_HoSKJ*100/1103127540)


	## UI:
# 103 photos analysées donc 103* 8 294 192 = 854301776 pixels

DenUIJ<-subset(UI,Espece=="Dendrilla antarctica")
MyUIJ<-subset(UI, Espece == "Mycale acerata")
KiUIJ<-subset(UI, Espece=="Kirkpatrickia variolosa")
HoUIJ<-subset(UI, Espece=="Homaxinella balfourensis")

S_DenUIJ<-colSums(DenUIJ [,-(1:2)])
P_DenUIJ<-c(S_DenUIJ*100/854301776)	
S_MyUIJ<-colSums(MyUIJ[,-(1:2)])
P_MyUIJ<-c(S_MyUIJ*100/854301776)
S_KiUIJ<-colSums(KiUIJ[,-(1:2)])
P_KiUIJ<-c(S_KiUIJ*100/854301776)
S_HoUIJ<-colSums(HoUIJ[,-(1:2)])
P_HoUIJ<-c(S_HoUIJ*100/854301776)
P_HoUIJ



	##Plot du taux de recouvrement des porifères
load("Demonspongetot.RData")

PlotRecouv<-ggplot(Demonspongetot[-(33:40),], aes(Station, Taux, fill=Espece), na.rm=TRUE)+
  geom_bar(stat="identity", alpha=0.9)+
  scale_y_continuous(limits=c(0, 10))+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  labs(x="Stations",
       y="Taux de recouvrement (en %)")

PlotRecouv







### ANALYSE DES TRAITS BIOLOGIQUES###
	#Fuzzy correspondance analysis: 
fuzzy<-list(tab = tab,species.names = name,moda.names = modaname, col.blocks = col.blocks)

prep.fuzzy.var(fuzzy$tab, fuzzy$col.blocks)




###ANALYSE DE REDONDANCE###
load(file="Moyen.RData")
load(file="Faune55.RData")


	#Transformation logarithmique des paramètres sédimentaires:
env1 <- log(Moyen[,-1])

	#Transformation de hellinger des taxons:
sp1.hel <- decostand(Faune55[,-1], "hellinger")  

spe.rda <- rda(sp1.hel ~ ., env1)
summary(spe.rda, display=NULL)



VERIFIER:
perc <- round(100*(summary(spe.rda)$cont$importance[2, 1:2]), 2)

plot(spe.rda)

anova.cca(spe.rda, step=1000)
anova.cca(spe.rda, step=1000, by="axis")








###BIAIS D'OBSERVATEUR

	## Test de Shapiro:
#H0: l'échantillon suit une loi normale
#Si H0 vrai-> Ttest
#Si H0 faux ->Mann Witney

	#Données de Jeanne:
shapiro.test(Richessespe_8) #p=0.9751 donc >0.05 donc on accepte H0 donc ça suit une loi normale
	#Données de Léa:
shapiro.test(RichessespeLea) #p=0.3695 donc >0.05 donc on accepte H0 donc ça suit une loi normale


	##Egalité des variances:
#H0: les variances sont égales
var.test(Richessespe_8, RichessespeLea)


	## Test T de  Student:
#H0: pas de différence significative entre les richesses spécifiques des 2 observateurs
ttest = t.test(Richessespe_8, RichessespeLea, var.equal=F) 
ttest ## p = 0.454 > 0.05 donc pas de différence significative au seuil alpha 0.05 avec une proba de 95%








### COURBES CUMULATIVES D'ESPECES
load(file="Sp2_8.Rdata")

	#Division des taxons par site:
MI_8<-subset(Sp2_8, station == "MI")
BI_8<-subset(Sp2_8, station=="BI")
FH_8<-subset(Sp2_8, station=="FH")
GR_8<-subset(Sp2_8, station=="GR")
HI_8<-subset(Sp2_8, station=="HI")
NH_8<-subset(Sp2_8, station=="NH")
SK_8<-subset(Sp2_8, station=="SK")
UI_8<-subset(Sp2_8, station=="UI")


	#Courbes d'accumulation:
AccMI_8<-specaccum(MI_8[,-1])
AccBI_8<-specaccum(BI_8[,-1])
AccFH_8<-specaccum(FH_8[,-1])
AccGR_8<-specaccum(GR_8[,-1])
AccHI_8<-specaccum(HI_8[,-1])
AccNH_8<-specaccum(NH_8[,-1])
AccSK_8<-specaccum(SK_8[,-1])
AccUI_8<-specaccum(UI_8[,-1])


#Plot:
par(mfcol=c(4,2))
plotMI_8<-plot(AccMI_8, col="blue", ci.type="poly", ci.col="lightblue", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
     le site MI", ylim=c(0,23))


plotBI_8<-plot(AccBI_8, col="blue", ci.type="poly", ci.col="yellow", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site BI", ylim=c(0,17))

plotFH<-plot(AccFH_8, col="blue", ci.type="poly", ci.col="purple", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site FH", ylim=c(0,17))

plotGR<-plot(AccGR_8, col="blue", ci.type="poly", ci.col="red", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site GR", ylim=c(0,18))

plotHI<-plot(AccHI_8, col="blue", ci.type="poly", ci.col="green", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site HI",ylim=c(0,17))

plotNH<-plot(AccNH_8, col="blue", ci.type="poly", ci.col="orange", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site NH", ylim=c(0,17))

plotSK<-plot(AccSK_8, col="blue", ci.type="poly", ci.col="brown", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site SK", ylim=c(0,19))

plotUI<-plot(AccUI_8, col="blue", ci.type="poly", ci.col="darkgreen", ylab="Nombre cumulatif d'espèces", xlab="Nombre de d'images analysées", 
               main= "Courbe cumulative des espèces dans 
               le site UI", ylim=c(0,15))




 ###CALCUL DES ESTIMATEURS CHAO JACKKNIFE ET BOOTSTRAP ###

specpool(MI_8[,-1])
poolMI_8<-poolaccum(MI_8[,-1])

specpool(BI_8[,-1])
poolBI_8<-poolaccum(BI_8[,-1])

specpool(FH_8[,-1])
poolFH_8<-poolaccum(FH_8[,-1])

specpool(GR_8[,-1])
poolGR_8<-poolaccum(GR_8[,-1])

specpool(HI_8[,-1])
poolHI_8<-poolaccum(HI_8[,-1])

specpool(NH_8[,-1])
poolNH_8<-poolaccum(NH_8[,-1])

specpool(SK_8[,-1])
poolSK_8<-poolaccum(SK_8[,-1])

specpool(UI_8[,-1])
poolUI_8<-poolaccum(UI_8[,-1])
