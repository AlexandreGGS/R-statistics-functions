# Utilisation RMD en results = 'asis' pour génération HTML

if(!require("knitr")) install.packages("knitr", repos="http://cran.us.r-project.org") ;
library(knitr)
if(!require("viridis")) install.packages("viridis", repos="http://cran.us.r-project.org") ;
library(viridis)
if(!require("plotly")) install.packages("plotly", repos="http://cran.us.r-project.org") ;
library(plotly)
if(!require("devtools")) install.packages("devtools", repos="http://cran.us.r-project.org")
library(devtools) 
if (!require("DT")) devtools::install_github("rstudio/DT")
library(DT)
if(!require("survival")) install.packages("survival", repos="http://cran.us.r-project.org")
library(survival)
if(!require("readxl")) install.packages("readxl", repos="http://cran.us.r-project.org")
library(readxl)
if(!require("lubridate")) install.packages("lubridate", repos="http://cran.us.r-project.org")
library(lubridate)
if(!require("pyramid")) install.packages("pyramid", repos="http://cran.us.r-project.org")
library(pyramid)
if(!require("tidyverse")) install.packages("tidyverse", repos="http://cran.us.r-project.org")
library(tidyverse)
#if(!require("kableExtra")) install.packages("kableExtra", repos="http://cran.us.r-project.org")
# library(kableExtra)

# kableExtra à utiliser si RMD en HTML
# BALISE à ajouter dans un fichier RMD
#<style>
#div.color { background-color:#ebf2f9;
#font-family: Verdana;}
#</style>


# Dépendances -------------------------------------------------------------

confint_prop_binom <- function(vecteur, pourcent=FALSE) {
    # retourne moyenne, borne inf, borne sup (entre 0 et 1)
    if( sum(!(vecteur %in% c(0,1)))>0) {
        cat("Erreur : valeurs autres que 0 et 1\n") ;
        print(as.matrix(table(vecteur))) ;
        mean <- 0 ;
        low_bound <- 0 ;
        upp_bound <- 0 ;
    } else {
        x <- sum(vecteur) ;
        n <- length(vecteur) ;
        temp <- binom.test(x, n, alternative ="two.sided", conf.level = 0.95) ;  
        mean <- temp$estimate[[1]] ;
        low_bound <- temp$conf.int[[1]] ;
        upp_bound <- temp$conf.int[[2]] ;
    }
    if( pourcent ) {
        return(c(
            paste(round(mean*100,2),"%", sep=""), 
            paste(round(low_bound*100,2),"%", sep=""), 
            paste(round(upp_bound*100,2),"%", sep="")
        )) ;
    } else {
        return(c(mean, low_bound, upp_bound))
    }
}

# effet_reel est la variable à prédire, avec valeurs 0 ou 1
# reponse est le nombre (quantitatif continu) issu du test ou du score à évaluer
# montre_seuils permet d'afficher un tableau avec les Se, Sp et d² aux seuils spécifiés dans seuils
# si montre_seuils=TRUE mais seuils=c(), tous les seuils constatés dans les données sont montrés.
# plot_seuils permet d'afficher certains seuils sur le graphique de la courbe ROC
courbe_roc <- function(effet_reel, reponse, montre_seuils=TRUE, seuils=c(), plot_seuils=c()  ) {
    if(!require("ROCR")) install.packages("ROCR", repos="http://cran.us.r-project.org") ;
    library(ROCR) ;
    
    temp <- data.frame(outcome=effet_reel, newvar=reponse) ;
    temp <- na.omit(temp) ;
    pred <- prediction( temp$newvar, temp$outcome) ;
    perf1 <- performance(pred, "sens", "spec") ;
    plot(perf1, print.cutoffs.at=plot_seuils, xlim=c(0,1), ylim=c(0,1)) ;
    perf2 <- performance(pred, "auc")
    cat(perf2@y.name, ":", perf2@y.values[[1]], "\n"  ) ;
    boxplot(temp$newvar ~ temp$outcome, col="#AAAAFF")
    # représentation de tous les seuils
    if( montre_seuils) {
        if( length(seuils)==0) {
            seuils <- unique(temp$newvar) ;
            seuils <- seuils[order(seuils)] ;
        }
        cat("Seuil", "Se","Sp", "VPP", "VPN", "d^2", "\n", sep="    ") ;
        for( seuil in seuils ) {
            vp <- nrow(temp[ temp$newvar>=seuil & temp$outcome==1, ]) ;
            vn <- nrow(temp[ temp$newvar<seuil  & temp$outcome==0, ]) ;
            fp <- nrow(temp[ temp$newvar>=seuil & temp$outcome==0, ]) ;
            fn <- nrow(temp[ temp$newvar<seuil  & temp$outcome==1, ]) ;
            se <- vp/(vp+fn) ;
            sp <- vn/(vn+fp) ;
            vpp <- vp/(vp+fp) ;
            vpn <- vn/(vn+fn) ;
            d2 <- (1-se)^2 + (1-sp)^2 ;
            cat(seuil, se,sp,vpp,vpn,d2, "\n", sep="    ") ;
        }
    }
}


# gold_standard et test sont deux vecteurs binaires de même longueur.
evalue_test <- function(gold_standard, test  ) {
    if(!require("psy")) install.packages("psy", repos="http://cran.us.r-project.org")
    library(psy) ;
    
    temp <- data.frame(outcome=gold_standard, newvar=test) ;
    nrow_av <- nrow(temp) ;
    temp <- na.omit(temp) ;
    nrow_ap <- nrow(temp) ;
    manquants <- nrow_av - nrow_ap ;
    if( manquants>0 ) {
        cat("Valeurs manquantes : n=", manquants,"soit",100*manquants/nrow_av,"%.\n") ;
    }
    plot(table(gold_standard, test), col="#AAAAFF") ;
    print(table(gold_standard, test)) ;
    
    vp <- nrow(temp[ temp$newvar==1 & temp$outcome==1, ]) ;
    vn <- nrow(temp[ temp$newvar==0  & temp$outcome==0, ]) ;
    fp <- nrow(temp[ temp$newvar==1 & temp$outcome==0, ]) ;
    fn <- nrow(temp[ temp$newvar==0  & temp$outcome==1, ]) ;
    se <- vp/(vp+fn) ;
    sp <- vn/(vn+fp) ;
    vpp <- vp/(vp+fp) ;
    vpn <- vn/(vn+fn) ;
    d2 <- (1-se)^2 + (1-sp)^2 ;
    kappa <- ckappa(temp) ;
    kappa <-  kappa$kappa ;
    auc <- se*sp + (1-se)*sp/2 + (1-sp)*se/2 ;
    if( se==0 & sp==0) {
        fscore <- 0 ;
    } else {
        fscore <- 1/((1/se + 1/vpp)/2) ;
    }
    
    cat("Effectifs :\n") ;
    cat("      Vrais positifs :", vp, "\n") ;
    cat("      Vrais négatifs :", vn, "\n") ;
    cat("      Faux positifs :", fp, "\n") ;
    cat("      Faux négatifs :", fn, "\n") ;
    cat("Métriques :\n") ;
    cat("      Sensibilité (rappel) :", se, "\n") ;
    cat("      SPécificité :", sp, "\n") ;
    cat("      Valeur prédictive positive (précision) :", vpp, "\n") ;
    cat("      Valeur prédictive négative :", vpn, "\n") ;
    cat("      F-mesure :", fscore, "\n") ;
    cat("      Kappa :", kappa, "\n") ;
    cat("      AUC (avec 1 point) :", auc, "\n") ;
    cat("      d² :", d2, "\n") ;
    
    plot(
        x=c(0, sp, 1),
        y=c(1, se, 0),
        xlab="Specificity",
        ylab="Sensitivity",
        type="l"
    ) ;
}


# Fonction Jo pour tester graphiquement la normalité
normPlot <- function(vecteur, etiquette, ...) {
    stopifnot(is.numeric(vecteur))
    op <- par(mfrow=c(1, 2), cex.main=0.9, ...)
    rangenorm <- c(mean(vecteur) - 3.72 * sd(vecteur), mean(vecteur) + 3.72 * sd(vecteur))
    dens <- dnorm(pretty(rangenorm, 200), mean = mean(vecteur), sd = sd(vecteur))
    plot(density(vecteur), xlab="", main=paste("Courbe de densité de la variable", etiquette, "\navec superposition d'une courbe de moyenne/ds identiques"), frame=F, lwd=2)
    lines(pretty(rangenorm, 200), dens, lty=2, lwd=2, col="red")
    qqnorm(vecteur, main=paste("Normal QQ Plot de", etiquette), frame=F)
    qqline(vecteur, col="red", lwd=2)
    par(op)
}




# Fonctions qui enleve les valeurs extremes qui peuvent poser problèmes pour les fonctions graphiques
remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * (qnt[2] - qnt[1])
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
}



# DESCRIPTIF --------------------------------------------------------------

######## Fonctions Desc avec output HTML (nécessite kableExtra)

# BINAIRE



desc_binaire_html <- function(vector, name="Variable", table=TRUE, old_graph = FALSE, plotly = TRUE, ...) {
	  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")
  name_html = paste0("<b>",name,"</b>")
  cat("<br><br>------------------------------------------------------------------------------------<br>",name_html,"<br>---------------<br>")
  
  if( length(unique(na.omit( vector ))) <2 ) {
    cat("Cette colonne comporte au plus une valeur et ne sera pas analysée<br>")
    return( TRUE )
  }
  
  vector <- as.numeric(vector) ;
  if( sum(is.na(vector))==0) {
    cat("Aucune valeur manquante.<br>") ;
  } else {
    cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.<br>") ;
    vector <- vector[!is.na(vector)] ;
  }
  cat("Effectif analysé :", length(na.omit(vector)),"<br>") ;
  cat("------------------------------------------------------------------------------------<br>") ;
  if( table ) {
    temp <- as.data.frame(table(vector)) ;
    modalites <- temp[,1]
    df_tmp = as.data.frame(modalites)
    for( une_modalite in modalites) {
      temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE)
      eff <- table(vector)[une_modalite]
      mean <- temp[1]
      IC <- paste0("[",temp[2],";",temp[3],"]")
      df_tmp$eff[df_tmp$modalites == une_modalite] = eff
      df_tmp$mean[df_tmp$modalites == une_modalite] = mean
      df_tmp$IC[df_tmp$modalites == une_modalite] = IC
    }
    cat("</div>")
    names(df_tmp)=c("Modalité", "Effectif","Proportion","IC95%")
    print(kableExtra::kbl(df_tmp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  }
  if(old_graph){
  	pie(table(vector)/length(vector), main=name, col=c("white", "cornflowerblue")) ;
  }
  if(plotly){
    df = vector %>% as_tibble() %>% group_by(value) %>% tally() %>% arrange(desc(value)) %>% mutate(colors = c("#2A2D6F", "#ACAFE8")) 
    fig = df %>% plot_ly(labels = ~value, sort = F, values = ~n,marker = list(colors=~colors))
    fig %>% add_pie(hole = 0.6) %>%
      layout(title = name,  showlegend = T,
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  }
}


# QUALI

desc_quali_html <- function(vector, name="Variable", table=TRUE, sort="decroissant", limit_chart=Inf, tronque_lib_chart=20, plotly = TRUE, old_graph = FALSE, plot_eff = TRUE, ...) {
	  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")
  name_html = paste0("<b>",name,"</b>")
  cat("<br><br>------------------------------------------------------------------------------------<br>",name_html,"<br>---------------<br>")
  
  if( length(unique(na.omit( vector ))) <2 ) {
    cat("Cette colonne comporte au plus une valeur et ne sera pas analysée<br>") ;
    return( TRUE ) ;
  }
  if( sum(is.na(vector))==0) {
    cat("Aucune valeur manquante.<br>") ;
  } else {
    cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.<br>") ;
    vector <- vector[!is.na(vector)] ;
  }
  cat("Effectif analysé :", length(na.omit(vector)),"<br>") ;
  cat("------------------------------------------------------------------------------------<br></div>") ;
  
  temp <- as.data.frame(table(vector)) ;
  modalites <- temp[,1] ;
  if( sort=="alpha") {
    modalites <- modalites[order(modalites)] ;
  } else if( sort=="croissant") {
    modalites <- modalites[order(temp[,2])] ;
  } else if( sort=="decroissant") {
    	  vector <- fct_infreq(as.factor(vector))
  }
  # tableau de contingence
  if( table ) {
    
    temp <- as.data.frame(table(vector)) ;
    modalites <- temp[,1]
    df_tmp = as.data.frame(modalites)
    for( une_modalite in modalites) {
      temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE)
      eff <- table(vector)[une_modalite]
      mean <- temp[1]
      IC <- paste0("[",temp[2],";",temp[3],"]")
      df_tmp$eff[df_tmp$modalites == une_modalite] = eff
      df_tmp$mean[df_tmp$modalites == une_modalite] = mean
      df_tmp$IC[df_tmp$modalites == une_modalite] = IC
    }
    names(df_tmp)=c("Modalité", "Effectif","Proportion","IC95%")
    print(kableExtra::kbl(df_tmp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  }
  
  # graphique à l'ancienne
  if(old_graph == TRUE){
  par(mar=c(4, 10, 4, 2) + 0.1) ;
  vector <- factor(vector, levels = modalites)
  temp <- rev(prop.table(table(vector)))
  id <- max(nrow(temp)-limit_chart+1,1):nrow(temp) ;
  # la première modalité traçée est celle la plus proche du point (0,0), on inverse donc nos vecteurs
  temp <- rev(temp[rev(id)])
  modalites <- names(temp)
  barplot(temp, horiz=TRUE, main=name, las=2, col="cornflowerblue", names.arg = substring(modalites, 1, tronque_lib_chart)) ;
  par(mar=c(5, 4, 4, 2) + 0.1) ;
}
# Graph plotly
if(plotly == TRUE){
  # récup pourcentage en numeric
  df_tmp = tibble(modalites=vector) %>% group_by(modalites) %>%
   summarise(nb = n(),y_plot = round(nb/nrow(.),2)*100,.groups = 'drop')
if(plot_eff == TRUE) {
df_tmp$y_plot = df_tmp$nb
y_lab_plot = "Effectifs"
} else {y_lab_plot = "Proportions"}
  gg = ggplot(df_tmp,aes(x=fct_rev(modalites),y=y_plot)) +
    geom_bar(stat = "identity", fill = "cornflowerblue") +
    scale_fill_viridis(discrete = T) +
    ggtitle(name) + xlab("Modalités") + ylab(y_lab_plot) +
    coord_flip() + theme_bw()
  ggplotly(gg)
}

}
#QUANTI_DISC
desc_quanti_disc_html = function(vector, name="Variable", mean_ci=TRUE, table=TRUE, Sum = T, sort="alpha", xlim=NULL, plotly = TRUE, old_plot = FALSE, ...) {
  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")
  name_html = paste0("<b>",name,"</b>")
  cat("<br><br>------------------------------------------------------------------------------------<br>",name_html,"<br>---------------<br>")
  
  if( length(unique(na.omit( vector ))) <2 ) {
    cat("Cette colonne comporte au plus une valeur et ne sera pas analysée<br>") ;
    return( TRUE ) ;
  }
  
  vector <- as.numeric(vector) ;
  if( sum(is.na(vector))==0) {
    cat("Aucune valeur manquante.<br>") ;
  } else {
    cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.<br>") ;
    vector <- vector[!is.na(vector)] ;
  }
  cat("Effectif analysé:", length(na.omit(vector)),"<br>") ;
  cat("------------------------------------------------------------------------------------<br></div>") ;
  if (Sum) {vector_noNA = na.omit(vector)
  df = as.data.frame(rbind(as.matrix(summary(vector_noNA)), Sd = sd(vector_noNA,na.rm=T)))
  row.names(df) = c("Minimum",
                    "1er quartile",
                    "Médiane",
                    "Moyenne",
                    "3eme quartile",
                    "Maximum",
                    "Ecart type")
  names(df) = name
  print(kableExtra::kbl(df) %>% kable_styling(bootstrap_options = "striped", full_width = F))} ;
  
  if( table ) {
    temp <- as.data.frame(table(vector)) ;
    modalites <- temp[,1] ;
    if( sort=="alpha") {
      modalites <- modalites[order(modalites)] ;
    } else if( sort=="croissant") {
      modalites <- modalites[order(temp[,2])] ;
    } else if( sort=="decroissant") {
      modalites <- modalites[order(0-temp[,2])] ;
    }
    
    df_tmp = as.data.frame(modalites)
    for( une_modalite in modalites) {
      temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE)
      eff <- table(vector)[une_modalite]
      mean <- temp[1]
      IC <- paste0("[",temp[2],";",temp[3],"]")
      df_tmp$eff[df_tmp$modalites == une_modalite] = eff
      df_tmp$mean[df_tmp$modalites == une_modalite] = mean
      df_tmp$IC[df_tmp$modalites == une_modalite] = IC
    }
    names(df_tmp)=c("Modalité", "Effectif","Proportion","IC95%")
    cat("<br><b><center>Détail des modalités</center></b><br>")
    print(kableExtra::kbl(df_tmp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
    
  }
  if( mean_ci ) {
    sd <- sd(vector) ;
    mean <- mean(vector) ;
    n <- length(vector) ;
    cat( "\n<br><div class = \"color\">Moyenne et intervalle de confiance à 95% :",  
         round(mean,2),"[",round(mean-1.96*sd/sqrt(n),2),";",round(mean+1.96*sd/sqrt(n),2),"]</div><br>")
  }
  if(old_plot) {
    plot(table(vector)/length(vector), xlab=name, ylab="proportion", col="cornflowerblue", xlim=xlim)
  }
  if(plotly) {
    t = tibble(x = vector) %>% group_by(x) %>% tally() 
    ggplotly(
      ggplot(t, aes(x=x, y=n)) +
        geom_segment( aes(x=x, xend=x, y=0, yend=n), color="black") +
        geom_point( color="cornflowerblue", size=3) +
        theme_light() +
        theme(panel.grid.major.x = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank()) + xlab(name) + ylab("Effectif") +
        scale_x_continuous(breaks = as.numeric(names(table(t$x))))
    )
   }
}


# QUANTI_CONT
desc_quanti_cont_html <- function(vector, name="Variable", mean_ci=TRUE, alpha=0.05, def_breaks = "Sturges", plot_boxplot = FALSE, plotly = TRUE, old_graph=FALSE, density=TRUE, ...) {
    cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")
	cat("<br>------------------------------------------------------------------------------------<br><b>",name,"</b><br>---------------<br>") ;
  if( length(unique(na.omit( vector ))) <2 ) {
    cat("Cette colonne comporte au plus une valeur et ne sera pas analysée\n") ;
    return( TRUE ) ;
  }
  vector <- as.numeric(vector) ;
  if( sum(is.na(vector))==0) {
    cat("Aucune valeur manquante.<br>") ;
  } else {
    cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.<br>") ;
    vector <- vector[!is.na(vector)] ;
  }
  cat("Effectif analysé :", length(na.omit(vector)),"<br>") ;
  cat("------------------------------------------------------------------------------------<br></div>") ;
  vector_noNA = na.omit(vector)
  df = as.data.frame(rbind(as.matrix(summary(vector_noNA)), Sd = sd(vector_noNA,na.rm=T)))
  row.names(df) = c("Minimum",
                    "1er quartile",
                    "Médiane",
                    "Moyenne",
                    "3eme quartile",
                    "Maximum",
                    "Ecart type")
  names(df) = name
  print(kableExtra::kbl(df) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  if( mean_ci ) {
    sd <- sd(vector) ;
    mean <- mean(vector) ;
    n <- length(vector) ;
    cat( "<br><div class = \"color\">>Moyenne et intervalle de confiance à ",100*(1-alpha),"% :",  
         round(mean,2),"[",round(mean+qnorm(alpha/2)*sd/sqrt(n),2),";",round(mean+qnorm(1-alpha/2)*sd/sqrt(n),2),"]</div><br>")
  }
  if (old_graph) {
    hist(vector, col="cornflowerblue", freq = FALSE, main = name, xlab = name, breaks = def_breaks, ylim = c(0,max(hist(vector, plot = F)$density, density(vector)$y)), ...)
    if (density) {
      lines(density(vector), col="red") ;
    }
  }
  if (plot_boxplot) {
    t = tibble(x = vector)
    gg_box = ggplot(t, aes(x ="",y=x)) +
    geom_violin(fill = "cornflowerblue") +
    coord_flip() + #geom_boxplot(width=0.1) + 
    ggtitle(name) + ylab("Valeur") + theme_bw() + xlab("")
    ggplotly(gg_box)
  }
  if (plotly) {
    t = tibble(x = vector)
    gg = ggplot(t, aes(x)) +
    geom_histogram(aes(y = stat(density)), fill = "cornflowerblue", binwidth = function(x) 2 * IQR(x) / (length(x)^(1/3))) +
    ggtitle(name) + xlab("Valeur") + theme_bw() +
    geom_line(stat="density", color = "red")
   ggplotly(gg)
  }
}



desc_quanti_cont_delai <- function(vector, name="Variable", mean_ci=TRUE, alpha=0.05, def_breaks = "Sturges", plot_boxplot = F, graph=TRUE, density=TRUE, ...) {
    cat(name,"\n") ;
    if( length(unique(na.omit( vector ))) <2 ) {
        cat("Cette colonne comporte au plus une valeur et ne sera pas analysée\n") ;
        return( TRUE ) ;
    }
    vector <- as.numeric(vector) ;
    if( sum(is.na(vector))==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.\n") ;
        vector <- vector[!is.na(vector)] ;
    }
    cat("Effectif analysé :", length(na.omit(vector)),"\n") ;
    cat("------------------------------------------------------------------------------------\n") ;
     
    if( mean_ci ) {
        sd <- sd(vector) ;
        mean <- mean(vector) ;
        n <- length(vector) ;
        cat( "\nMoyenne et intervalle de confiance à ",100*(1-alpha),"% :",  
             round(mean,2),"[",round(mean+qnorm(alpha/2)*sd/sqrt(n),2),";",round(mean+qnorm(1-alpha/2)*sd/sqrt(n),2),"].\n",
             "\nCalcul des IC",100*(1-alpha),"% à partir du théorème central limite") ;
    }
    if (graph) {
        hist(vector, col="cornflowerblue", freq = FALSE, main = name, xlab = name, breaks = def_breaks, ylim = c(0,max(hist(vector, plot = F)$density, density(vector)$y)), ...)
        if (density) {
            lines(density(vector), col="red") ;
        }
    }
    if (plot_boxplot) {
        boxplot(vector, main = name, ylab = "", col = "cornflowerblue")
    }
}




# extrait_dataframe = dataframe sous la forme de variables binaires
# ex1 : mondataframe[,c("var1", "var2")]
# ex2 : cbind(vecteur1, vecteur2)
# limit_chart permet de limiter le nombre de sorties dans le graphique
desc_quali_multi <- function(extrait_dataframe, name="Variable", sort="decroissant", limit_chart=Inf, tronque_lib_chart=20, ...) {
    if(!require("lattice")) install.packages("lattice", repos="http://cran.us.r-project.org") ;
    library(lattice) ;
    list_na <- apply(extrait_dataframe, 1, function(x) {sum(is.na(x))>0})
    
    cat(name,"\n") ;
    
    if( sum(list_na)==0) {
        cat("Aucune valeur manquante.\n")
        extrait_df <- extrait_dataframe
    } else {
        cat("Au moins une observation possède une valeur manquante pour l'une des variables étudiées, celle(s)-ci  est(/sont) exclue(s) pour toutes les variables.\n")
        cat("Observations exclues : n=", sum(list_na),"soit",round(100*mean(list_na),2),"%.\n") ;
        extrait_df <- extrait_dataframe
        extrait_dataframe <- na.omit(extrait_dataframe)
            }
    cat("Effectif analysé :", nrow(extrait_df)-sum(list_na),"\n") ;
    cat("------------------------------------------------------------------------------------\n") ;
    
    liste_modalites <- names(extrait_dataframe) ;
    if( sort=="alpha") {
        liste_modalites <- names(extrait_dataframe) ;
        liste_modalites <- liste_modalites[order(liste_modalites)] ;
    } else {
        effectifs <- c() ;
        for( une_modalite in liste_modalites ) {
            effectifs <- c(effectifs, sum(extrait_dataframe[une_modalite])) ;
        }
        if( sort=="croissant" ) {
            liste_modalites <- liste_modalites[order(effectifs)] ;
        } else if( sort=="decroissant" ) {
            liste_modalites <- liste_modalites[order(effectifs, decreasing=TRUE)] ;
        }
    }
    cat("\nDans le tableau et le graphique ci-dessous, le total excède généralement 100% car un individu peut avoir plusieurs modalités simultanément.\n") ;
    cat("\nLe calcul des IC95% est realisé à l'aide d'une loi binomiale\n")
    cat("\nModalite\tEffectif\tProportion\tIC95%\n") ;
    nb <- c() ;
    nb_ib <- c() ;
    nb_sb <- c() ;
    for( une_modalite in liste_modalites) {
        binaire <- extrait_dataframe[,une_modalite] ;
        temp <- confint_prop_binom(vecteur=binaire, pourcent=TRUE) ;
        temp_nb <- confint_prop_binom(vecteur=binaire, pourcent=FALSE) ;
        cat(une_modalite,"\t",max(0, table(binaire)[ifelse(is.logical(binaire), "TRUE", "1")], na.rm=T),"\t",temp[1],"\t[",temp[2],";",temp[3],"]\n") ;
        nb <- c(nb, temp_nb[1]) ;
        nb_ib <- c(nb_ib, temp_nb[2]) ;
        nb_sb <- c(nb_sb, temp_nb[3]) ;
    }
    cat("------------------------------------------------------------------------------------\n") ;
    cat("\nEffectifs en fonction du nombre de modalités :\n") ;
    eff_moda <- as.data.frame(table(apply(extrait_dataframe,1,sum)))
    names(eff_moda) <- c("Nombre de modalités", "Effectif")
    print(eff_moda, row.names = FALSE)  
    # maintenant un graphe avec Lattice
    #label, nb, nb_ib, nb_sb, xlim, main, xlab, ref=-999, signif=NULL) 
    id <- 1:min(limit_chart, length(liste_modalites))
    #max(length(liste_modalites)-limit_chart+1,1):length(liste_modalites) ;
    nb <- nb[id]
    nb_ib <- nb_ib[id]
    nb_sb <- nb_sb[id]
    liste_modalites <- substr(liste_modalites,1,tronque_lib_chart)
    label2 <- factor(liste_modalites[id], levels = rev(liste_modalites[id]))
    xlim <- c(  0 , max(nb_sb)+0.1 )
    plot(dotplot(label2 ~ nb, xlim=xlim, panel=function(x,y) {    
        #scales = list(x = list(log = 2)),
        #panel.abline(v=ref,col='#CCCCCC', lty=2)
        panel.barchart(x=nb, y=y, col="cornflowerblue")
        #panel.xyplot(x,y,pch=16,cex=1,col='navy')
        panel.segments(nb_ib,as.numeric(y),nb_sb,as.numeric(y),lty=1,col='navy')
        panel.segments(nb_ib,as.numeric(y)+0.15,nb_ib,as.numeric(y)-0.15,lty=1,col='navy')
        panel.segments(nb_sb,as.numeric(y)+0.15,nb_sb,as.numeric(y)-0.15,lty=1,col='navy')
    }, xlab="proportion", main=name)) ;
}


# Multivariee -------------------------------------------------------------


# Valeurs pour method : pearson (paramétrique), spearman (non paramétrique)
bivarie_quanti_quanti_html = function(x, y, xname="Variable quantitative 1", yname="Variable quantitative 2", method="pearson", droite_reg=TRUE) {
  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style>")
  cat("<br><div class = \"color\">Analyse bivariée : Méthode du coefficient de corrélation de ",method," sur deux variables quantitatives<br>")
  coeff_r <- cor(x, y, method=method, use="complete.obs" ) ;
  coeff_r2 <- coeff_r^2 ;
  temp <- cor.test(x, y, method=method, alternative="two.sided") ;
  pval <- temp$p.value ;
  cat("Coefficient de correlation de", method, ": r=",coeff_r,", avec r²=", coeff_r2,
      "<br><strong> La p-value = ",pval,"<br></strong></div>") ;
  plot(x=x, y=y, type="p", pch=20, xlab=xname, ylab=yname) ;
  if( droite_reg ) {
    reg <- glm(y~x) ;
    b <- reg$coefficients[1] ;
    a <- reg$coefficients[2] ;
    cat("<br><div class = \"color\">Equation droite : y =",a,"* x +",b,"</div><br>") ;
    points_x <- c(min(na.omit(x)), max(na.omit(x))) ;
    points_y <- a*points_x + b ;
    lines(x=points_x, y=points_y, type="l", col="blue") ;
  }
}


# Valeurs pour method : chisq (semi-paramétrique), fisher (non paramétrique mais attention au nombre de modalités)
bivarie_quali_quali_html <- function(x, y, xname="Variable qualitative 1", yname="Variable qualitative 2", method="chisq", prop.table = TRUE, table = FALSE, mosaic = T,...) {
  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")
  nb_valide <- nrow(na.omit(data.frame(x,y)))
  nb_manquant <- length(x) - nb_valide
  cat(paste0(yname, " en fonction de ", xname, ".<br>")) ;
  if (nb_manquant==0) {
    cat("Aucune valeur manquante.<br>") ;
  } else {
    cat("Valeurs manquantes : n=", nb_manquant,"soit",100*nb_manquant/length(x),"%.<br>") ;
  }
  cat(paste0("Effectif analysé : ", nb_valide, ".<br>")) ;
  cat("------------------------------------------------------------------------------------<br>") ;
  
  if( method=="chisq") {
    cat("Analyse bivariée : Test du Chi 2 sur deux variables qualitatives<br>")
    obj <- chisq.test(table(x,y)) ;
    print(obj) ;
    pval <- obj$"p.value" ;
  } else if( method=="fisher") {
    cat("Analyse bivariée : Test de Fisher sur deux variables qualitatives<br>")
    obj <- fisher.test(table(x,y), workspace = 200000000) ;
    print(obj) ;
    pval <- obj$"p.value" ;
  } else {
    cat("Nom de méthode incorrect :", method,"<br>") ;
    return("Erreur !") ;
  }
  cat("<br>------------------------------------------------------------------------------------<br>")
  cat("<strong>La p-value (petit p) de ce test = ",pval,"</strong>")
  cat("<br>------------------------------------------------------------------------------------<br></div><br>")
  if (prop.table) {
    prop_table <-round(100*prop.table(table(y,x),2),2)  
    tableau_perc <- as.data.frame.matrix(prop_table)
    cat(paste0("<br><center><strong>Proportions de ", yname, " (lignes) par modalité de ",xname," (colonnes)</center></strong><br>"))
    print(kableExtra::kbl(tableau_perc) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  }
  if (table) {
    tableau_eff <- as.data.frame.matrix(table(y,x))
    cat("<br><center><strong>Tableau des effectifs</center></strong><br>")
    print(kableExtra::kbl(tableau_eff) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  }
  
  if (mosaic) {plot(table(x, y), xlab=xname, ylab=yname, main="Mosaicplot", col="cornflowerblue")} ;
}

# VERSION 2 AVEC UN MOSAICPLOT PLUS JOLI				   
bivarie_quali_quali_html_v2 <- function(
  x, y,
  xname = "Variable qualitative 1",
  yname = "Variable qualitative 2",
  method = "chisq",
  prop.table = TRUE,
  table = FALSE,
  mosaic = TRUE,                 # <- si TRUE, on affiche le ggplot (ex-mosaic)
  inner_max = 0.85,              # taille max (relative) du carré intérieur
  text_size = 5,                 # taille du % au centre
  axis_text_size = 18,           # taille du texte des axes
  ...
) {
  # Dépendances pour le graphique
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Veuillez installer ggplot2")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Veuillez installer dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Veuillez installer tidyr")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Veuillez installer viridis")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Veuillez installer scales")
  # Tables jolies (optionnelles, déjà présentes dans ta version)
  if ((prop.table || table) && !requireNamespace("kableExtra", quietly = TRUE)) {
    warning("kableExtra n'est pas installé : les tableaux seront affichés bruts.")
  }

  # --- En-tête HTML identique à ta version
  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")

  nb_valide <- nrow(na.omit(data.frame(x, y)))
  nb_manquant <- length(x) - nb_valide
  cat(paste0(yname, " en fonction de ", xname, ".<br>"))
  if (nb_manquant == 0) {
    cat("Aucune valeur manquante.<br>")
  } else {
    cat("Valeurs manquantes : n=", nb_manquant, " soit ", 100 * nb_manquant / length(x), "%.<br>")
  }
  cat(paste0("Effectif analysé : ", nb_valide, ".<br>"))
  cat("------------------------------------------------------------------------------------<br>")

  # --- Tests (inchangés)
  if (method == "chisq") {
    cat("Analyse bivariée : Test du Chi 2 sur deux variables qualitatives<br>")
    obj <- chisq.test(table(x, y))
    print(obj)
    pval <- obj[["p.value"]]
  } else if (method == "fisher") {
    cat("Analyse bivariée : Test de Fisher sur deux variables qualitatives<br>")
    obj <- fisher.test(table(x, y), workspace = 200000000)
    print(obj)
    pval <- obj[["p.value"]]
  } else {
    cat("Nom de méthode incorrect :", method, "<br>")
    return("Erreur !")
  }
  cat("<br>------------------------------------------------------------------------------------<br>")
  cat("<strong>La p-value (petit p) de ce test = ", pval, "</strong>")
  cat("<br>------------------------------------------------------------------------------------<br></div><br>")

  # --- Tableaux (inchangés)
  if (prop.table) {
    prop_table <- round(100 * prop.table(table(y, x), 2), 2)
    tableau_perc <- as.data.frame.matrix(prop_table)
    cat(paste0("<br><center><strong>Proportions de ", yname, " (lignes) par modalité de ", xname, " (colonnes)</center></strong><br>"))
    if (requireNamespace("kableExtra", quietly = TRUE)) {
      print(kableExtra::kbl(tableau_perc) %>% kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE))
    } else {
      print(tableau_perc)
    }
  }
  if (table) {
    tableau_eff <- as.data.frame.matrix(table(y, x))
    cat("<br><center><strong>Tableau des effectifs</center></strong><br>")
    if (requireNamespace("kableExtra", quietly = TRUE)) {
      print(kableExtra::kbl(tableau_eff) %>% kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE))
    } else {
      print(tableau_eff)
    }
  }

   # --- VISU (remplace le mosaic plot)
  if (isTRUE(mosaic)) {
    suppressPackageStartupMessages({
      library(dplyr); library(tidyr); library(ggplot2); library(viridis); library(scales)
    })

    d <- data.frame(x = x, y = y) %>% tidyr::drop_na()

    # Proportions EXACTEMENT comme tableau_perc
    prop_tbl <- prop.table(table(d$y, d$x), 2)
    lab_df <- as.data.frame(as.table(prop_tbl)) %>%
      dplyr::rename(y = Var1, x = Var2, p_exact = Freq) %>%
      dplyr::mutate(
        label = paste0(round(100 * p_exact, 0), "%")  # arrondi à 0 décimale
      )

    tab <- d %>%
      dplyr::count(x, y, name = "n") %>%
      tidyr::complete(x, y, fill = list(n = 0)) %>%
      dplyr::left_join(lab_df, by = c("x", "y")) %>%
      dplyr::mutate(
        side = sqrt(p_exact) * inner_max
      )

    p_plot <-
      ggplot2::ggplot(tab, ggplot2::aes(x = x, y = y)) +
      # Carré blanc fixe
      geom_tile(fill = "white", color = "grey70", width = 0.98, height = 0.98) +
      # Carré intérieur proportionnel et coloré
      geom_tile(aes(width = side, height = side, fill = p_exact)) +
      # Texte = arrondi à 0 décimale
      geom_text(
        aes(label = label, color = p_exact < 0.5),
        size = text_size, fontface = "bold", family = "Century Gothic"
      ) +
      # Palette viridis avec légende arrondie à 0 décimale
      scale_fill_viridis_c(
        option = "D", direction = -1,
        labels = function(x) paste0(round(100 * x, 0), "%"),
        name = "Part"
      ) +
      scale_color_manual(values = c("white", "black"), guide = "none") +
      labs(x = xname, y = yname) +
      coord_fixed() +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_text(family = "Century Gothic", face = "bold", size = axis_text_size),
        axis.title.y = element_text(family = "Century Gothic", face = "bold", size = axis_text_size),
        axis.text.x  = element_text(family = "Century Gothic", face = "bold", size = axis_text_size, angle = 45),
        axis.text.y  = element_text(family = "Century Gothic", face = "bold", size = axis_text_size),
        legend.title = element_text(family = "Century Gothic", face = "bold"),
        legend.text  = element_text(family = "Century Gothic")
      )

    print(p_plot)
  }

}
				   
				   
# Valeurs pour method : 
# - student (2groupes, paramétrique) ; 
# - anova (plus de 2 groupes, paramétrique) ;
# - wilcoxon (Mann-Whitney, 2 groupes, non paramétrique) ;
# - kruskal (Kruskal-Wallis, plus de 2 groupes, non paramétrique)
library(plotly)
bivarie_quali_quanti_html <- function(x, y, xname="Variable qualitative", yname="Variable quantitative", method="student", cond.app = F, graph = F, ggraph = T, ...) {
  cat("<style>
div.color { background-color:#ebf2f9;
font-family: Verdana;}
</style><br><div class = \"color\">")
  if(!is.factor(x)) {x=as.factor(x)}
  nb_valide <- nrow(na.omit(data.frame(x,y)))
  nb_manquant <- length(x)-nb_valide
  cat(paste0(yname, " en fonction de ", xname, ".<br>")) ;
  if (nb_manquant==0) {
    cat("Aucune valeur manquante.<br>") ;
  } else {
    cat("Valeurs manquantes : n=", nb_manquant,"soit",100*nb_manquant/length(x),"%.<br>") ;
  }
  cat(paste0("Effectif analysé : ", nb_valide, ".<br>")) ;
  cat("------------------------------------------------------------------------------------<br>") ;
  
  if (cond.app == T) {
    cat("Vérification des conditions d'application :<br>")
    cat("1. Normalité de la distribution dans chaque groupe (méthode graphique). Voir graphiques.<br>")
    var_min <- numeric()
    var_max <- numeric()
    for (une_modalite in unique(x)) {
      normPlot(y[x == une_modalite], une_modalite)
      ifelse(length(na.omit(y[x == une_modalite])) >= 30, cat("Note : l'Effectif de l'échantillon est supérieur ou égal à 30 (hors NA).<br>"), cat(""))
    }
    var_min <- min(c(var_min, sd(y[x == une_modalite])^2, na.rm = T))
    var_max <- max(c(var_max, sd(y[x == une_modalite])^2, na.rm = T))
    cat("2. Egalité des variances.<br>")
    cat(paste("Le rapport de la variance max / variance min = ", round(var_max/var_min, 2), "(doit être inférieur à 1,5)"))
    cat("<br>------------------------------------------------------------------------------------<br>")
  }
  
  if( method=="student") {
    cat("Analyse bivariée : Test t de Student d'une variable qualitative avec une variable quantitative<br>")
    obj <- t.test(y~x, var.equal=FALSE);
    pval <- as.numeric(obj$p.value);
    print(obj) ;
  } else if( method=="anova") {
    cat("Analyse bivariée : Test d'analyse de variances (ANOVA) d'une variable qualitative avec une variable quantitative<br>")
    obj <- aov(y~x);
    pval <- as.numeric(summary(obj)[[1]][["Pr(>F)"]][1])
    print(summary(obj)) ;
  } else if( method=="wilcoxon") {
    cat("Analyse bivariée : Test de Wilcoxon d'une variable qualitative avec une variable quantitative<br>")
    obj <-wilcox.test(y~x);
    pval<-obj$p.value;
    print(obj) ;
  } else if( method=="kruskal") {
    cat("Analyse bivariée : Test de Kruskal-Wallis d'une variable qualitative avec une variable quantitative<br>")
    obj <-kruskal.test(x=y, g=x, na.action=na.rm);
    pval<-obj$p.value;
    print(obj) ;
  } else {
    cat("Nom de méthode incorrect :", method,"<br>") ;
    return("Erreur !") ;
  }
  cat("<strong><br>------------------------------------------------------------------------------------<br>")
  cat("La p-value (petit p) de ce test = ",pval)
  cat("<br>------------------------------------------------------------------------------------<br></strong></div>")
  maxcar = max(nchar(as.character(unique(x))),na.rm = T)
  if (method %in% c("student", "anova")) {
    cat("<br><center><strong>Moyennes avec IC95% pour chaque modalité de la variable qualitative</center></strong><br>")
    df_temp = tibble(quali = NA,
                     n = NA,
                     mean_quanti = NA,
                     ic = NA)
    for (une_modalite in sort(unique(na.omit(x)))) {
      quali <- une_modalite
      n <- sum(!is.na(y[x == quali]))
      mean_quanti <- mean(y[x == quali], na.rm = TRUE)
      mean_quanti <- round(mean_quanti,2)
      low_bound <- round(mean_quanti-1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
      upp_bound <- round(mean_quanti+1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
      df_temp = bind_rows(df_temp,tibble(quali,
                                         n,
                                         mean_quanti,
                                         ic = paste0("[",low_bound,";",upp_bound,"]"))
      )
    }
    df_temp = df_temp %>% filter(!is.na(quali)) %>% select(`Modalité` = quali,
                                                           `Effectif` = n,
                                                           `Moyenne` = mean_quanti,
                                                           `IC95%` = ic)
    print(kableExtra::kbl(df_temp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  } else if (method %in% c("wilcoxon", "kruskal")) {
    cat("<br><center><strong>Médianes avec intervalle inter-quartile [Q1;Q3] pour chaque modalité de la variable qualitative</center></strong><br>")
    df_temp = tibble(quali = NA,
                     n = NA,
                     med_quanti = NA,
                     ic = NA)
    for (une_modalite in sort(unique(na.omit(x)))) {
      quali <- une_modalite
      n <- sum(!is.na(y[x == quali]))
      med_quanti <- median(y[x == quali], na.rm = TRUE)
      med_quanti <- round(med_quanti,2)
      low_bound <- round(quantile(y[x == quali], probs = 0.25, na.rm = TRUE), 2)
      upp_bound <- round(quantile(y[x == quali], probs = 0.75, na.rm = TRUE), 2)
      df_temp = bind_rows(df_temp,tibble(quali,
                                         n,
                                         med_quanti,
                                         ic = paste0("[",low_bound,";",upp_bound,"]"))
      )
    }
    df_temp = df_temp %>% filter(!is.na(quali)) %>% select(`Modalité` = quali,
                                                           `Effectif` = n,
                                                           `Médiane` = med_quanti,
                                                           `Intervalle interquartile` = ic)
    print(kableExtra::kbl(df_temp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
    
    cat("<br><center><strong>Moyennes avec IC95% pour chaque modalité de la variable qualitative</center></strong><br>")
    df_temp = tibble(quali = NA,
                     n = NA,
                     mean_quanti = NA,
                     ic = NA)
    for (une_modalite in sort(unique(na.omit(x)))) {
      quali <- une_modalite
      n <- sum(!is.na(y[x == quali]))
      mean_quanti <- mean(y[x == quali], na.rm = TRUE)
      mean_quanti <- round(mean_quanti,2)
      low_bound <- round(mean_quanti-1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
      upp_bound <- round(mean_quanti+1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
      df_temp = bind_rows(df_temp,tibble(quali,
                                         n,
                                         mean_quanti,
                                         ic = paste0("[",low_bound,";",upp_bound,"]"))
      )
    }
    df_temp = df_temp %>% filter(!is.na(quali)) %>% select(`Modalité` = quali,
                                                           `Effectif` = n,
                                                           `Moyenne` = mean_quanti,
                                                           `IC95%` = ic)
    print(kableExtra::kbl(df_temp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
  }
  
  if ( graph ) {
    boxplot(y~x, xlab=xname, ylab=yname, col="cornflowerblue", ...) ;
  }
  
  if( ggraph ) {
    df= tibble(x = x, y = y)
    df$x = fct_reorder(df$x,df$y,na.rm = T)
    ggplotly(
      ggplot(data = df, aes(x = x, y = y, fill = x)) +
             geom_boxplot() +
             xlab(xname) + ylab(yname) + theme_minimal() +
             theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 7)) +
             scale_fill_viridis_d() + coord_flip()
    ) %>% layout(showlegend = FALSE) # no legend
  }
}

## Test bivarié t de student pour groupes appariés
bivarie_quali_quanti_app <- function(x, y, xname="Variable qualitative", yname="Variable quantitative", method="student", graph = TRUE, ...) {
  if( method=="student") {
  cat("Analyse bivariée : Test t de Student d'une variable qualitative avec une variable quantitative avec appariement\n")
    obj <- t.test(y~x, var.equal=FALSE, paired=TRUE);
    pval <-as.numeric(obj$p.value);
    print(obj) ;
    cat("\n------------------------------------------------------------------------------------\n")
    cat("La p-value (petit p) de ce test = ",pval)
    cat("\n------------------------------------------------------------------------------------\n")
     } else if(method=="wilcoxon"){
    cat("Analyse bivariée : Test de Wilcoxon d'une variable qualitative avec une variable quantitative (Groupes appariés)\n")
    obj <-wilcox.test(y~x, paired = TRUE, conf.int = TRUE);
    pval<-obj$p.value;
    print(obj) ;
    cat("\n------------------------------------------------------------------------------------\n")
    cat("La p-value (petit p) de ce test = ",pval)
  }
     if ( graph ) {
  boxplot(y~x, xlab=xname, ylab=yname, col="#AAAAFF", ...) ;
 
}
}

## Test de MacNemar : Groupes appariés quali quali
bivarie_quali_quali_app <- function(x, y, xname="Variable qualitative", yname="Variable qualitative", graph = TRUE, ...) {
  cat("Analyse bivariée : Test du Chi2 de Mac Nemar de 2 variables qualitatives avec appariement\n")
  mac <- table(x,y)
  obj <- mcnemar.test(mac)
  pval <-as.numeric(obj$p.value)
 print(obj) ;
    cat("\n------------------------------------------------------------------------------------\n")
    cat("La p-value (petit p) de ce test = ",pval)
    cat("\n------------------------------------------------------------------------------------\n")
    cat("Tableau de contingence\n")
    print(table(x,y,dnn = c(xname,yname)))
     if ( graph ) {
  plot(table(x, y), xlab=xname, ylab=yname, main="Mosaicplot", col="#AAAAFF", ...)  ;
  }
}


# Ajout des bivar sur les multivalués = représentation graphique de plusieurs variables binaires croisées avec une autre variable
bivarie_quali_quali_multi <- function(extrait_dataframe, vecteur, yname="Variable multivaluee", vname = "Variable quali", sort="alpha", limit_chart=Inf, tronque_lib_chart_y=8, tronque_lib_chart_v=4, ...) {
    
    if(!require("lattice")) install.packages("lattice", repos="http://cran.us.r-project.org") ;
    library(lattice) ;
    list_na <- apply(cbind(extrait_dataframe, vecteur), 1, function(x) {sum(is.na(x))>0})
    
    cat(yname,"\n") ;
    
    liste_facteur <- unique(vecteur)
    n_mod_fact <- length(liste_facteur)
    n_var_bin <- ncol(extrait_dataframe)
    n_combinaison <- n_var_bin * length(liste_facteur)
    
    
    if( sum(list_na)==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Au moins une observation possède une valeur manquante pour l'une des variables étudiées, celle(s)-ci  est(/sont) exclue(s) pour toutes les variables.\n")
        cat("Observations exclues : n=", sum(list_na),"soit",round(100*mean(list_na),2),"%.\n") ;
        extrait_dataframe <- extrait_dataframe[!list_na,] ;
        vecteur <- vecteur[!list_na]
    }
    cat("Effectif analysé :", nrow(extrait_dataframe)-sum(list_na),"\n") ;
    cat("------------------------------------------------------------------------------------\n") ;
    
    liste_modalites <- names(extrait_dataframe) ;
    if( sort=="alpha") {
        liste_modalites <- names(extrait_dataframe) ;
        liste_modalites <- liste_modalites[order(liste_modalites)] ;
    } else {
        effectifs <- data.frame() ;
        for( une_modalite in liste_modalites ) {
            for (un_fact in unique(vecteur)) {
                selection_ligne <- vecteur %in% un_fact
                binaire <- extrait_dataframe[,une_modalite]
                effectifs <- c(effectifs, table(extrait_dataframe[selection_ligne,une_modalite], vecteur[selection_ligne])[ifelse(is.logical(binaire), "TRUE", "1"),])
            }
        }
        if( sort=="croissant" ) {
            liste_modalites <- liste_modalites[order(effectifs)] ;
        } else if( sort=="decroissant" ) {
            liste_modalites <- liste_modalites[order(effectifs, decreasing=TRUE)] ;
        }
    }
    cat("\nDans le tableau et le graphique ci-dessous, les pourcentages sont calculés pour chaque modalité de ", vname, ". Toutefois, le total peut excèder 100% car un individu peut avoir plusieurs modalités simultanément.\n") ;
    cat("\nLe calcul des IC95% est realisé à l'aide d'une loi binomiale\n")
    nb <- data.frame() ;
    nb_ib <- data.frame() ;
    nb_sb <- data.frame() ;
    for( une_modalite in liste_modalites) {
        cat("\nVariable binaire : ", une_modalite, ".\n")
        cat("\nModalite\tEffectif\tProportion\tIC95%\n") ;
        nb_i <- c() ;
        nb_ib_i <- c() ;
        nb_sb_i <- c() ;
        for(un_fact in unique(vecteur)) {
            selection_ligne <- vecteur %in% un_fact
            binaire <- extrait_dataframe[selection_ligne,une_modalite] ;
            temp <- confint_prop_binom(vecteur=binaire, pourcent=TRUE) ;
            temp_nb <- confint_prop_binom(vecteur=binaire, pourcent=FALSE) ;
            cat(un_fact,"\t",max(0, table(binaire)[ifelse(is.logical(binaire), "TRUE", "1")], na.rm=T),"\t",temp[1],"\t[",temp[2],";",temp[3],"]\n") ;
            nb_i <- c(nb_i, temp_nb[1]) ;
            nb_ib_i <- c(nb_ib_i, temp_nb[2]) ;
            nb_sb_i <- c(nb_sb_i, temp_nb[3]) ;
        }
        nb <- rbind(nb, data.frame(y=une_modalite, x=unique(vecteur), valeur=nb_i)) ;
        nb_ib <- rbind(nb_ib, data.frame(y=une_modalite, x=unique(vecteur), valeur=nb_ib_i)) ;
        nb_sb <- rbind(nb_sb, data.frame(y=une_modalite, x=unique(vecteur), valeur=nb_sb_i)) ;
    }
    cat("------------------------------------------------------------------------------------\n") ;
    
    # Tableau des effectifs en fonction du nombre de modalité
    cat("\nEffectifs en fonction du nombre de modalités :\n") ;
    eff_moda <- as.matrix(table(apply(extrait_dataframe,1,sum), vecteur, dnn = c("",vname)))
    print(eff_moda, row.names = FALSE)  
    
    # maintenant un graphe avec Lattice
    id <- 1:min(limit_chart, length(liste_modalites))
    nb <- nb[nb$y %in% liste_modalites[id],]
    nb_ib <- nb_ib[nb_ib$y %in% liste_modalites[id],]
    nb_sb <- nb_sb[nb_sb$y %in% liste_modalites[id],]
    liste_modalites <- substr(liste_modalites,1,tronque_lib_chart_y)
    liste_facteur <- substr(liste_facteur, 1, tronque_lib_chart_v)
    nb$y <- factor(substr(nb$y,1,tronque_lib_chart_y), levels = liste_modalites[id])
    xlim_par <- c(  0 , max(nb_sb$valeur)+0.1 )
    plotTop <- max(nb_sb$valeur) + 0.1
    tabbedNb <- tapply(nb$valeur, list(nb$x, nb$y), function(x) c(x = x))
    tabbedIb <- tapply(nb_ib$valeur, list(nb_ib$x, nb_ib$y), function(x) c(x = x))
    tabbedSb <- tapply(nb_sb$valeur, list(nb_sb$x, nb_sb$y), function(x) c(x = x))
    barCenters <- barplot(height = tabbedNb,
                          beside = TRUE, las = 1,
                          ylim = c(0, plotTop),
                          cex.names = (1-(n_var_bin/20))+0.05,
                          main = paste0(yname, " en fonction de ", vname),
                          ylab = "Proportion",
                          xlab = yname,
                          border = "black", axes = TRUE,
                          legend.text = TRUE,
                          args.legend = list(title = vname, 
                                             x = "top",
                                             cex = (1-(n_mod_fact/20))+0.05,
                                             horiz = T))
    
    arrows(barCenters, tabbedIb, barCenters,
           tabbedSb, lwd = 1.5, angle = 90,
           code = 3, length = 1/n_combinaison)
}


# Bivarié survie : Par rapport à desc_survie : ajout argument "vector_strate" et "name_strate"
# - Argument “colors_strate” : character, default to “pointille”,
#     couleurs pour les différentes courbes : “pointille”, “grays”, “rainbow”
# - Argument “scale_legend“ : default to 1, numeric
# - Argument “method”, character, default to “logrank”
#     “logrank” ou “cox”
# - Argument order_legend, logical, default to TRUE
#     Ordonne la légende par ordre décroissant de probabilités de survie. 
#     Peut être source d'erreurs, a désactiver si nécessaire.

bivarie_survie <- function(vector_evt, vector_delai, vector_strate,  name="Variable", 
                           name_strate = "strate", method = "logrank",
                           colors_strate = "pointille",
                           scale_legende = 1, order_legend = TRUE,
                           xmax=NULL, ymin=0, output_prop_survie=c()) {
    
    
    cat(name,"\n") ;
    vector_evt <- as.numeric(vector_evt) ;
    vector_delai <- as.numeric(vector_delai) ;
    vector_strate <- as.factor(vector_strate)
    
    manquant_strate <- is.na(vector_strate)
    vector_evt <- vector_evt[!manquant_strate]
    vector_delai <- vector_delai[!manquant_strate]
    vector_strate <- vector_strate[!manquant_strate]
    
    for (modalite_strate in levels(vector_strate)){
        manquant <- (is.na(vector_evt[vector_strate %in% modalite_strate]) + 
                         is.na(vector_delai[vector_strate %in% modalite_strate]))>0
        
        if( sum(manquant)==0 ) {
            cat("Aucune valeur manquante dans la strate", modalite_strate,"de la variable", name_strate,".\n") ;
        } else {
            cat("Valeurs manquantes dans la strate", modalite_strate,
                " de la variable ", name_strate,
                ": n=", sum(manquant),"soit",100*mean(manquant),"%.\n") ;
            vector_evt[vector_strate %in% modalite_strate] <- vector_evt[vector_strate %in% modalite_strate & !manquant] ;
            vector_delai[vector_strate %in% modalite_strate] <- vector_delai[vector_strate %in% modalite_strate & !manquant] ;
        }
        
    }
    
    
    
    
    effet_surv <- Surv(event=vector_evt, time=vector_delai, type="right") ;
    obj_survfit2 <- survfit( effet_surv~vector_strate, conf.int=TRUE) ;
    
    
    # Pour chaque strate on s'interesse aux résultats dans la strate en question
    resultats_totaux <- summary(obj_survfit2)
    resultats_totaux <- 
        data.frame(time = resultats_totaux$time,
                   n.risk = resultats_totaux$n.risk,
                   n.event = resultats_totaux$n.event,
                   surv = resultats_totaux$surv,
                   upper = resultats_totaux$upper,
                   lower = resultats_totaux$lower,
                   strata = sub(pattern = "vector_strate=", replacement = "",x = resultats_totaux$strata))
    
    
    if (order_legend){
        max_time <- ceiling(max(resultats_totaux$time) / 2)
        legend_survival <-  resultats_totaux[resultats_totaux$time < max_time,c("time","surv", "strata")]
        legend_survival_list <- split(x = legend_survival, f = legend_survival$strata)
        legend_survival_list2 <- lapply(legend_survival_list, function(x){
            x <- x[order(x$time, decreasing = TRUE),]
            x <- x[1,]
        })
        legend_survival2 <- do.call(rbind, legend_survival_list2)
        legend_survival2 <- as.vector(legend_survival2$strata[order(legend_survival2$surv, decreasing = TRUE)])
        
        order_legend <- legend_survival2
        
    } else {
        order_legend <- levels(vector_strate)
    }
    
    
    if(colors_strate %in% "pointille"){
        colors_strate <- "black"
        line_type <- 1: length(levels(vector_strate))
        
        colors_strate_legend <- "black"
        line_type_legend <- line_type
        names(line_type_legend) <- levels(vector_strate)
        line_type_legend <- line_type_legend[order_legend]
        
    } else if (colors_strate %in% "rainbow"){
        colors_strate <- rainbow(n = length(levels(vector_strate)))
        line_type <- 1
        
        colors_strate_legend <- colors_strate
        names(colors_strate_legend) <- levels(vector_strate)
        colors_strate_legend <- colors_strate_legend[order_legend]
        line_type_legend <- 1
        
    } else if (colors_strate %in% "grays"){
        colors_strate <- grey.colors(n = length(levels(vector_strate)), 
                                     start = 0.2, end = 0.7)
        line_type <- 1
        
        colors_strate_legend <- colors_strate
        names(colors_strate_legend) <- levels(vector_strate)
        colors_strate_legend <- colors_strate_legend[order_legend]
        line_type_legend <- 1
    } else {
        colors_strate <- "black"
        line_type <- 1: length(levels(vector_strate))
        
        colors_strate_legend <- "black"
        line_type_legend <- line_type
        names(line_type_legend) <- levels(vector_strate)
        line_type_legend <- line_type_legend[order_legend]
        
    }
    
    plot( obj_survfit2, xmax=xmax, ymin=ymin, main=name, 
          col = colors_strate, lty = line_type) ;
    
    legend("topright", legend = order_legend, 
           title = "Strates", 
           col = colors_strate_legend, lty = line_type_legend,
           border = F, cex = scale_legende)
    
    cat("Informations sur la survie (dont la médiane) :\n")
    print( obj_survfit2) ;
    # on peut ajouter des % de survie à telles dates...
    
    
    # strates <- sub(pattern = "vector_strate=", replacement = "", names(obj_survfit2$strata))
    
    for (strates_calcul in order_legend){
        resultats_strates <- resultats_totaux[resultats_totaux$strata %in% strates_calcul,]
        
        if( length(output_prop_survie)>0) {
            cat("\n Survies estimées par la méthode de Kaplan-Meier dans la strate '", strates_calcul,"' : \n", sep = "") ;
            # en général on a peu d'évenements donc il est intéressant de supprimer les censures
            y <- c() ;    y_upper <- c() ;    y_lower <- c() ;
            x <- c() ;
            for( i in length(resultats_strates$time) : 1 ) {
                if(resultats_strates$n.event[i]>=1) {
                    x <- c(x, resultats_strates$time[i]) ;
                    y <- c(y, resultats_strates$surv[i]) ;
                    y_upper <- c(y_upper, resultats_strates$upper[i]) ;
                    y_lower <- c(y_lower, resultats_strates$lower[i]) ;
                }
            }
            # maintenant on parcourt
            for( a_time in output_prop_survie) {
                if( a_time>max(resultats_strates$time)) {
                    cat("- à", a_time,": survie inconnue (suivi trop court)\n") ;
                } else {
                    start <- TRUE ;
                    for( i in 1:length(x) ) {
                        if(x[i]<a_time & !start) {
                            cat("- à t=", a_time,"(depuis t=", x[i], "), survie estimée à", y[i], "[",y_lower[i],";",y_upper[i],"]\n") ;
                            break ;
                        }
                        start <- FALSE ;
                    }
                }
            }
        }
    }
    
    if (method %in% "logrank") {
        
        cat("\n Test du Log Rank \n")
        log_rank <- survdiff(effet_surv~vector_strate)
        Effectifs <- sum(log_rank$n)
        tableau_recap <- data.frame(
            Modalites = levels(vector_strate),
            Effectifs = as.vector(log_rank$n),
            Observed = log_rank$obs,
            Expected = log_rank$exp
        )
        pvalue <- 1 - pchisq(log_rank$chisq,1)
        
        
        
        cat("\n Effectifs : ", Effectifs, "\n\n")
        print(tableau_recap)
        cat("\n pvalue (petit p) = ", pvalue,"\n")
        
        
        # print(survdiff(effet_surv~vector_strate))
        
        
    } else if (method %in% "cox") {
        cat("\n Modèle de Cox \n")
        cox <- coxph(effet_surv~vector_strate)
        summary_cox <- summary(cox)
        tableau_recap <- data.frame(
            Modalites = levels(vector_strate),
            "Hazard Ratio" = c(1, round(unname(summary_cox$conf.int[,"exp(coef)"]),2)),
            "lower conf.int" = c(" - ", round(unname(summary_cox$conf.int[,"lower .95"]),2)),
            "upper conf.int" = c(" - ", round(unname(summary_cox$conf.int[,"upper .95"]), 2)),
            pvalue = c("Reference", unname(summary_cox$coefficients[,"Pr(>|z|)"]))
        )
        
        cat("\n Effectifs = ", summary_cox$n, "; et nombre d'évènements = ", summary_cox$nevent, "\n")
        print(tableau_recap)
        
    }
}


# dataframe est l'extrait de dataframe à analyser : toutes les colonnes de cet argument seront incluses.
# method peut valoir "pearson", "kendall" ou "spearman"
# non_sig : que faut-il faire des coeff avec p>5% ? "montre", "cache", "barre". "cache" est par défaut, recommandé.

matrice_correlation <- function(dataframe, method="pearson", non_sig="cache") {
    if(!require("corrplot")) install.packages("corrplot", repos="http://cran.us.r-project.org") ;
    library(corrplot) ;
    
    nrow_av <- nrow(dataframe) ;
    dataframe <- na.omit(dataframe) ;
    nrow_ap <- nrow(dataframe) ;
    manquants <- nrow_av - nrow_ap ;
    if( manquants>0 ) {
        cat("Valeurs manquantes. Lignes ignorées pour la matrice de corrélation : n=", manquants,"soit",100*manquants/nrow_av,"%.\n") ;
    }
    M <- cor(dataframe, method=method) ;
    # maintenant on regarde et stocke la significativité de chaque coeff de corrélation
    conf.level <- 0.95 ;
    n <- ncol(M) ;
    p.mat <- matrix(NA, n, n) ; # lowCI.mat <- uppCI.mat <- matrix(NA, n, n) ;
    diag(p.mat) <- 0 ;   #  diag(lowCI.mat) <- diag(uppCI.mat) <- 1 ;
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test( M[, i], M[, j], conf.level=conf.level,  method=method) ;
            p.mat[i,j] <- p.mat[j,i] <- tmp$p.value ; # lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1] ; uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2] ;
        }
    }
    # Puis on trace le graphique
    cat("Le graphique ci-dessous représente le coefficient de corrélation de",method,"de chaque couple de variables. Le disque est d'autant plus grand et foncé que la corrélation est forte (en valeur absolue). La couleur bleue indique une corrélation positive. La couleur rouge indique une corrélation négative (variation en sens inverse).\n") ;
    if( non_sig=="montre") {
        cat("Ce graphique affiche tous les coefficients, même ceux non significativement différents de zéro.\n") ;
        corrplot(M, method="circle", type="lower", order="AOE", tl.cex=0.8) ;
        
    } else if( non_sig=="cache") {
        cat("Les coefficients non significativement différents de zéro sont masqués.\n") ;
        corrplot(M, p.mat=p.mat, insig="blank", method="circle", type="lower", order="AOE", tl.cex=0.8) ;
        
    } else if( non_sig=="barre") {
        cat("Les coefficients non significativement différents de zéro sont barrés.\n") ;
        corrplot(M, p.mat=p.mat, sig.level=0.05, method="circle", type="lower", order="AOE", tl.cex=0.8) ;
        
    }
    cat("Matrice de corrélation (coefficients de corrélation de", method,") :\n") ;
    print(M) ;
}

# Pyramide des âges

plot_pyramide_ages <- function(age, sexe, largeur = 5, vec_axis_x = NULL){
    if(!require("pyramid")) install.packages("pyramid", repos="http://cran.us.r-project.org") ;
    library(pyramid)
    hommes <- sexe == 1
    ages <- cut(age, breaks=c(-Inf, seq(from = largeur, to = 100, by = largeur), Inf), include.lowest = TRUE, right = FALSE, 
                labels=c("0 - 4",
                         "5 - 9",
                         "10 - 14",
                         "15 - 19",
                         "20 - 24",
                         "25 - 29",
                         "30 - 34",
                         "35 - 39",
                         "40 - 44",
                         "45 - 49",
                         "50 - 54",
                         "55 - 59",
                         "60 - 64",
                         "65 - 69",
                         "70 - 74",
                         "75 - 79",
                         "80 - 84",
                         "85 - 89",
                         "90 - 94",
                         "95 - 100",
                         "100 +")) ;
    repartition_hommes <- as.data.frame(table(ages[hommes]))
    repartition_femmes <- as.data.frame(table(ages[!hommes]))
    repartition <- merge(x = repartition_hommes, y = repartition_femmes, by = "Var1")  
    repartition <- rename(repartition, age = Var1, hommes = Freq.x, femmes = Freq.y)   
    nombre_max <- max(repartition$hommes, repartition$femmes, na.rm = TRUE)
    data <- data.frame(repartition_hommes[,2], repartition_femmes[,2], repartition_femmes[,1])
    pyramid(data, Lcol="cornflowerblue", Rcol="plum1", Laxis= vec_axis_x, Cadj =-0.008, Cgap = 0.18, Csize = 0.9, AxisFM = "d")
    
    
    #fenetres <- rbind(c(0, 0.45, 0, 1), c(0.45, 1, 0, 1))
    #split.screen(figs = fenetres)
    #screen(1)
    #barplot(repartition_hommes[,2], rep.int(5,21), axisnames=FALSE, horiz=TRUE,
    #xlim = c(nombre_max, 0), xlab = "Effectif")
    #title("Hommes")
    #screen(2)
    #barplot(repartition_femmes[,2],rep.int(5,21), horiz=TRUE, names.arg=repartition_femmes[,1],
    #xlim = c(0, nombre_max), xlab = "Effectif", cex.names = 0.7)
    #title("Femmes")
    #axis
    #close.screen(all=TRUE)
}


reg_cox_valide_risques_prop <- function(obj_cox) {
    if(!require("rms")) install.packages("rms", repos="http://cran.us.r-project.org") ;
    library(rms) ;
    cat("Pour valider la proportionnalité des risques, les p valeurs doivent être supérieures à 5% :\n")
    cox.zph(obj_cox)
}



# fonction pour tracer la moyenne mobile de Y en fonction de X.
# idéal pour explorer le caractère linéaire d'une relation, et proposer des seuils.
# y peut être binaire ou quantitative. X doit être quantitative.
bivarie_quanti_quanti_smooth <- function(x, y, halfwindow, xname="X", yname="Y") {
    selector <- complete.cases(x, y) ;
    x <- x[selector] ;
    y <- y[selector] ;
    if( length(x)>800) {
        pch <- 46 ;
    } else {
        pch <- 20 ;
    }
    
    plot(x=x, y=y, type="p", xlab=paste(xname," (n=",length(x), ")", sep=""), ylab=yname, col="#555555", pch=pch) ;
    smooth_x <- c();
    smooth_y <- c();
    smooth_y_inf <- c();
    smooth_y_sup <- c();
    values <- unique(x) ;
    values <- values[order(values)] ;
    for( i in values) {
        smooth_x <- c(smooth_x, i) ;
        # sample_x <- x[abs(x - i) <= halfwindow] ;
        sample_y <- y[abs(x - i) <= halfwindow] ;
        y_mean <- mean(sample_y) ;
        y_sd <- sd(sample_y) ;
        nb <- length(sample_y) ;
        smooth_y <- c(smooth_y, y_mean) ;
        smooth_y_inf <- c(smooth_y_inf, y_mean - 1.96*y_sd/sqrt(nb)) ;
        smooth_y_sup <- c(smooth_y_sup, y_mean + 1.96*y_sd/sqrt(nb)) ;
    }
    lines(x=smooth_x, y=smooth_y, type="l", col="blue", lwd=2) ;
    lines(x=smooth_x, y=smooth_y_inf, type="l", col="blue", lty="dotted", lwd=2) ;
    lines(x=smooth_x, y=smooth_y_sup, type="l", col="blue", lty="dotted", lwd=2) ;
}



arbre_decision_description <- function(obj_rpart, main="Arbre de décision") {
    if(!require("rpart")) install.packages("rpart", repos="http://cran.us.r-project.org") ;
    library(rpart) ;
    print(obj_rpart, digits=3, spaces=6) ;
    plot(obj_rpart) ; 
    text(obj_rpart, use.n=TRUE, col="blue", cex=0.7, digits=2, all=TRUE) ;
    
    # note that obj_rpart$where gives the number of the leaf, and obj_rpart$frame$yval the probability
    cases <- data.frame(ref=obj_rpart$y, leaf=obj_rpart$where)
    confiances <- data.frame(leaf=1:nrow(obj_rpart$frame), conf=obj_rpart$frame$yval) ;
    cases <- merge(x=cases, y=confiances, by.x="leaf", by.y="leaf")
    courbe_roc(effet_reel=cases$ref, reponse=cases$conf, montre_seuils=TRUE) ; 
}






# Fonction pour nettoyer un fichier importé avec read_execl

nettoyage_excel <- function(donnees){
    lignes_a_enlever <- apply(X = donnees, MARGIN = 1, FUN = function(x){all(is.na(x))})
    colonnes_a_enlever <- apply(X = donnees, MARGIN = 2, FUN = function(x){all(is.na(x))})
    
    donnees <-  donnees[!lignes_a_enlever, !colonnes_a_enlever]
    
    class(donnees) <- "data.frame"
    
    return(donnees)
}


# Fonction pour convertir les données chiffrées comprenant des espaces (séparateur millier) et convertir la "," en "." et le €

str_conv <- function(x) {
	x <- gsub("\u20AC","",as.character(x))
	x <- gsub(" ","",x)
	x <- gsub(",",".",x)
	#x <- trimws(x)
	x <- as.numeric(x)
	return(x)
}

# Fonction de connexion à PMSIpilot

session_pmsi_pilot <- function() {
	library(rvest)
	url       <-"https://pmsipilot/pmsipilot/"
ppsession <-html_session(url)               ## creation d'une session
ppform    <-html_form(ppsession)[[1]]       ## récupération des champs à remplir (pour rentrer le login)
filled_form <- set_values(ppform,
                          `utilisateur[login]` = "", 
                          `utilisateur[motdepasse]` = "") # Paramétrer les champs login et mot de passe PMSIPILOT
session <- submit_form(ppsession,filled_form)
return(session)
}

# Configuration du proxy et supression des certificats SSL

set_proxy <- function() {
	library(httr)
	set_config(use_proxy(url="http://proxy:8080/", port=8080, username="",password="")) # Proxy
	set_config( config( ssl_verifypeer = 0L ) ) # Pour éviter l'erreur : 'Peer certificate cannot be authenticated with given CA certificates'
} 


# PLOT pour regression

plot_odds<-function(x, title = NULL){
  temp <- summary(x)$coefficients
tmp<-data.frame(cbind(exp(coef(x)), exp(confint(x))),signif = ifelse(temp[,4]<0.05,1,0))
odds<-tmp[-1,]
names(odds)<-c("OR", "lower", "upper","signif")
odds$vars<-row.names(odds)
ticks<-c(seq(.1, 1, by =.1), seq(0, 10, by =1), seq(10, 100, by =10))
ggplot(odds, aes(y= OR, x = reorder(vars, OR),color=factor(signif) )) +
geom_point(size = 4) +
geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,size =1) +
scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
labs(title = title, x = "Variables", y = "OR") +
scale_color_manual(values=c("black","red")) +
guides(fill=FALSE, color=FALSE) +
theme_bw()
}
	
	
################### Copier un dataframe directement dans le clipboard pour colelr dans excel
copy_excel <- function(df, sep="\t", dec=",", max.size=(200*1000)){
    # Copy a data.frame to clipboard
    write.table(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
  }
	
################### Faire un spread comme tidyverse mais avec plusieurs variables
myspread <- function(df, key, value) {
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
        unite(temp, !!keyq, variable) %>%
        spread(temp, value)
}
		       
		 






