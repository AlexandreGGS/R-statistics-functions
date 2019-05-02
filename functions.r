if(!require("knitr")) install.packages("knitr", repos="http://cran.us.r-project.org") ;
library(knitr)
if(!require("broman")) install.packages("broman", repos="http://cran.us.r-project.org") ;
library(broman)   # utiliser myround(nombre)
opts_chunk$set(fig.width=6, fig.height=4, fig.path='figs/', warning=FALSE, message=FALSE, 
               echo=TRUE, warning=FALSE)
set.seed(53079239)
if(!require("qtl")) install.packages("qtl", repos="http://cran.us.r-project.org")
library(qtl) ;
if(!require("dplyr")) install.packages("dplyr", repos="http://cran.us.r-project.org") ;
library(dplyr) ;
if(!require("devtools")) install.packages("devtools", repos="http://cran.us.r-project.org")
library(devtools) ;
if (!require("DT")) devtools::install_github("rstudio/DT")
library(DT) ;
if(!require("survival")) install.packages("survival", repos="http://cran.us.r-project.org")
library(survival) ;
if(!require("ggplot2")) install.packages("ggplot2", repos="http://cran.us.r-project.org")
library(ggplot2)
#if(!require("rgeos")) install.packages("rgeos", repos="http://cran.us.r-project.org")
#library(rgeos)
if(!require("readxl")) install.packages("readxl", repos="http://cran.us.r-project.org")
library(readxl)
if(!require("lubridate")) install.packages("lubridate", repos="http://cran.us.r-project.org")
library(lubridate)
if(!require("pyramid")) install.packages("pyramid", repos="http://cran.us.r-project.org")
library(pyramid)
if(!require("pander")) install.packages("pander", repos="http://cran.us.r-project.org")
library(pander)
if(!require("kableExtra")) install.packages("kableExtra", repos="http://cran.us.r-project.org")
library(kableExtra)

#############################################################################################

# DEPENDANCES

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


# insert un saut de ligne dans une
saut_ligne <- function(chaine, limit = 30) {
    
}


#############################################################################################
#############################################################################################
#############################################################################################

# DATA MANAGEMENT

# terminologie : on gère :
# - pour l'ATIH, "grandes" terminologies : cim10, ccam, ghm, rghm (racines de ghm), cmd, ucd, lpp
# - pour l'ATIH, "petites" terminologies : atih_mode_entree, atih_mode_sortie, atih_prov, atih_dest,  atih_auto_um
# - pour l'ATIH : finess (utilisé de manière simplifiée)
# - autres : departement
# paste : si TRUE, dans le résultat on accole le code, un espace, et le libellé. Si FALSE, renvoit le libellé tout seul.
ajoute_libelle <- function(vector, terminologie="rien", paste=TRUE) {
    
    # on charge les terminologies complémentaires
    if(terminologie %in% c("cim10", "ccam", "ghm", "rghm", "cmd", "ucd", "ucd_dci", "lpp")) {
        charge_ressource("atih_termino_all") ;
        ajoute_libelle <- TRUE ;
    } else if(terminologie %in% c("atih_mode_entree", "atih_mode_sortie", "atih_prov", "atih_dest",  "atih_auto_um")) {
        charge_ressource("atih_terminologie_pmsi_mco") ;
        ajoute_libelle <- TRUE ;
    } else if(terminologie == "departement") {
        charge_ressource("insee_noms_departements") ;
        ajoute_libelle <- TRUE ;
    } else if(terminologie %in% c("finess", "status","secteur")) {
        charge_ressource("atih_identifiants_etablissements")
        ajoute_libelle <- TRUE ;
    } else {
        ajoute_libelle <- FALSE ;
    }
    
    # on ajoute libelle
    if(terminologie == "cim10") {
        libelle <- atih_terminologie_cim10[vector,mult="first",libelle] ;
    } else if(terminologie == "ccam") {
        libelle <- atih_terminologie_ccam[vector,mult="first",libelle] ;
    } else if(terminologie == "ghm") {
        libelle <- atih_terminologie_ghm[vector,mult="first",libelle] ;
    } else if(terminologie == "rghm") {
        libelle <- atih_terminologie_rghm[vector,mult="first",libelle] ;
    } else if(terminologie == "cmd") {
        libelle <- atih_terminologie_cmd[vector,mult="first",libelle] ;
    } else if(terminologie == "ucd") {
        libelle <- atih_terminologie_ucd[vector,mult="first",libelle] ;
    } else if(terminologie == "ucd_dci") {
        libelle <- atih_terminologie_ucd_dci[vector,mult="first",libelle] ;
    } else if(terminologie == "lpp") {
        libelle <- atih_terminologie_lpp[vector,mult="first",libelle] ;
    } else if(terminologie == "atih_mode_entree") {
        libelle <- atih_terminologie_pmsi_mco[list("entree_mode",vector),mult="first",libelle] ;
    } else if(terminologie == "atih_mode_sortie") {
        libelle <- atih_terminologie_pmsi_mco[list("sortie_mode",vector),mult="first",libelle] ;
    } else if(terminologie == "atih_prov") {
        libelle <- atih_terminologie_pmsi_mco[list("entree_prov",vector),mult="first",libelle] ;
    } else if(terminologie == "atih_dest") {
        libelle <- atih_terminologie_pmsi_mco[list("sortie_destination",vector),mult="first",libelle] ;
    } else if(terminologie == "atih_auto_um") {
        libelle <- atih_terminologie_pmsi_mco[list("autoum",vector),mult="last",libelle] ;
    } else if(terminologie == "departement") {
        libelle <- insee_noms_departements[vector,mult="first",libelle] ;
    } else if(terminologie == "finess") {
        libelle <- atih_identifiants_etablissements[vector,mult="last",raison_sociale] ;
    } else if(terminologie == "status") {
        libelle <- atih_identifiants_etablissements[vector,mult="last",status] ;
    } else if(terminologie == "secteur") {
        libelle <- atih_identifiants_etablissements[vector,mult="last",secteur] ;
    }
    if( paste ) {
        return(paste(vector, libelle)) ;
    } else {
        return(libelle) ;
    }
}
# Cette fonction prend un vecteur de variable qualitative multivaluée, et retourne un dataframe 
# d'autant de lignes, avec une variable binaire par modalité initiale.
# seules les modalités listées dans "values" seront traitées, les autres seront ignorées (la fonction ne crée pas une modalité "autres" par exemple)
# "newnames" doit préciser les noms des nouvelles variables binaires.
# Sinon, "newnames" peut préciser un préfixe (valeur unique), auquel on ajoutera "_" puis le nom de chaque modalité.
# Une valeur initiale "NA" donne autant de "0" (et non "NA") dans les variables binaires.
split_quali_multi <- function(vector, values, newnames="variable_quali", separateur="_") {
  if(newnames == "") {
    newnames <- as.character(values)
  }
    if(length(newnames)==1 | length(newnames) != length(values) ) {
        newnames <- paste(newnames, values, sep=separateur) ;    
    }
    values <- as.character(values) ;
    vector <- as.character(vector) ;
    if( sum(is.na(vector))>0) {
        cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.\n") ;
        cat("On considèrera que la réponse est valide et vaut toujours non. Ces enregistrement ne sont pas exclus.\n") ;
        vector[is.na(vector)] <- "" ;
    }
    new_df <- NULL ;
    for( i in 1:length(values)) {
        une_modalite <- values[i] ;
        nom_variable <- newnames[i] ;
        binaire <- as.numeric(regexpr(une_modalite,vector, ignore.case=TRUE, perl=FALSE)>-1) ;
        if( is.null(new_df)) {
            new_df <- as.data.frame(binaire);
            names(new_df) <- nom_variable ;
        } else {
            old_names <- names(new_df) ;
            if( nom_variable %in% names(new_df)  ) {
                binaire <- binaire + new_df[nom_variable] ;
                binaire <- as.numeric(binaire>0) ;
                new_df[nom_variable] <- binaire ;
            } else {
                new_df <- cbind(new_df, binaire) ;
                names(new_df) <- c(old_names, nom_variable) ;
            }
        }
    }
    return(new_df) ;
}

# Fonction permettant à partir d'un dataframe avec des modalités en colonnes, de sortir un vecteur conservant la modalité qui nous intéresse pour chaque ligne
# -- dataframe = input (tableau binaire avec les modalités en colonnes)
# -- sorted_modalities = vecteur comprenant les modalités triées (nom complet correspondant aux colonnes du df)
# -- new_modalities = vecteur comprenant les modalités triées (format de sortie)
aggregate_w_priority <- function(dataframe, sorted_modalities, new_modalities,  empty_as_na=FALSE) {
    new_vector <- rep("", nrow(dataframe)) ;
    for( i in length(sorted_modalities):1 ) {
        a_modality <- sorted_modalities[i]
        a_new_name <- new_modalities[i]
        new_vector[which(dataframe[,a_modality]==1)] <- a_new_name ;
    }
    if(empty_as_na) {
        new_vector[which(new_vector=="")] <- NA ;
    }
    return( new_vector) ;
}



# Fonction permettant d'effectuer les regroupements de codes à partir d'un fichier csv structuré en 2 colonnes (codes, regroupement)  
# L'argument variable_a_regrouper est un vecteur qui contient les codes à regrouper
# 
# L'argument chemin_fichier correspond au chemin du fichier .csv (separateur ;) 
# Le tableau doit au minimum contenir deux variables : "code" et "regroupement". La première contient la liste des codes et la seconde le nom du regroupement correspondant
# 
# L'argument v_a_r_multivalue indique si la variable à regrouper est multivaluée ou non
# Si oui, la sortie sera un dataframe de variables binaires ou un vecteur (parametre output = "df"/"vector"))
# Si non, la sortie sera un vecteur d'une variable qualitative à plusieurs modalitées
# MAJ : ajout d'une 3e colonne priority dans le tableau. Ainsi qu'un argument priority = T/F permettant de gérer les priorités de regroupement lorsque l'on travaille avec une 
# exact : Si TRUE, renvoie TRUE uniquement si variable_a_regrouper match exactement avec code. Si FALSE, renvoie TRUE également si variable_a_regrouper contient code.
regroupement_modalite <- function(variable_a_regrouper, chemin_fichier, v_a_r_multivalue = TRUE, output = "df", priority=F, exact = F, split_char = ";") {
    
    regroupements <- read.csv2(chemin_fichier, stringsAsFactors = FALSE)
    
    if (!v_a_r_multivalue | output=="vector") {
        data_return <- rep(NA, length(variable_a_regrouper))
    } else {
        list_var <- unique(regroupements$regroupement)
        data_return <- as.data.frame(matrix(data = NA, 
                                            nrow = length(variable_a_regrouper), 
                                            ncol = length(list_var), 
                                            dimnames = list(NULL, list_var)))
    }
    
    stopifnot(c("code","regroupement") %in% names(regroupements))
    
    # Recodage des modalités NA en autres
    regroupements$regroupement[is.na(regroupements$regroupement)] <- "Autres"
    
    # Priorisation des modalités
    if (!"priority" %in% names(regroupements)) {
        regroupements$priority <- 1
    }  
    
    temp <- aggregate(priority~regroupement, data = regroupements, FUN = min)
    list_mod <- temp[order(temp$priority, decreasing = T), "regroupement"]    
    if (v_a_r_multivalue) {
        variable_a_regrouper <- strsplit(variable_a_regrouper, split = split_char, fixed = T)
    }
    
    for (i in list_mod) {
        list_a_matcher <- regroupements$code[regroupements$regroupement == i]
        
        if (!exact) {
            list_a_matcher <- paste0("(", paste(list_a_matcher, collapse = ")|("), ")", collapse = "")
        }
        
        if (v_a_r_multivalue & output == "df") {
            if (exact) {
                temp_return <- lapply(variable_a_regrouper, function(x) x %in% list_a_matcher)
                temp_return <- lapply(temp_return, sum) > 0
                data_return[,i] <- temp_return
            } else {
                data_return[,i] <- grepl(list_a_matcher, variable_a_regrouper)
            }
            
        } else {
            if (exact) {
                if (v_a_r_multivalue & output == "vector") {
                    temp_return <- lapply(variable_a_regrouper, function(x) x %in% list_a_matcher)
                    temp_return <- lapply(temp_return, sum) > 0
                } else {
                    temp_return <- variable_a_regrouper %in% list_a_matcher
                }
            } else {
                temp_return <- grepl(list_a_matcher, variable_a_regrouper)
            }
            data_return[temp_return] <- i
        }
    }
    data_return[data_return == FALSE] <- 0
    data_return[data_return == TRUE] <- 1
    return(data_return)
}


#' @title Trouve department
#' @description Retourne le numéro de département, depuis un numéro de département, un code géographique ATIH, un code postal, ou un numéro FINESS
#' @param codegeo character A character vector with codegeo codes. 
#' @param type character Préciser ici le type de départ : cp, codegeo_atih, finess sont gérés à ce jour. A étendre si besoin. Tout autre valeur => non-modification du code.
#' @return A \code{character} vector with INSEE's department code.
#' @references 
#' \enumerate{
#'      \item \href{http://www.atih.sante.fr/mise-jour-2014-de-la-liste-de-correspondance-codes-postaux-codes-geographiques}{ATIH codegeo}
#'      \item \href{http://www.insee.fr/fr/methodes/nomenclatures/cog/telechargement/2014/txt/depts2014.txt}{INSEE's department codes}
#' }

trouve_departement <- function(geocode, type="departement") {
    if( ! "character" %in% class(geocode) ) geocode <- as.character(geocode) ;
    
    if( type=="finess") {
        
        nb_long_irreg <- sum(nchar(geocode)!=9 & nchar(geocode)!=8) ;
        if(nb_long_irreg>0) {
            cat("Attention,",nb_long_irreg,"codes ne faisaient pas 9 caractères (malgré l'ajout de zéros initiaux). Le résultat ci-dessous n'est pas garanti.\n") ;
        }
        selecteur <- which(nchar(geocode)==8) ;
        geocode[selecteur] <- paste("0",geocode[selecteur], sep="" ) ;
        dept <- substr(geocode, 1, 2) ;
        
        # Traitement spécial pour les DOM (dépt=97)
        dom <- dept == "97" 
        dept[dom] <- paste0(dept[dom], substr(geocode[dom], 4, 4))
        # traitement spécial pour Mayotte
        dept[dept == "98"] <- "976" ;
        
    } else if( type=="codegeo_atih" | type=="cp") {
        
        nb_long_irreg <- sum(nchar(geocode)!=5 & nchar(geocode)!=4) ;
        if(nb_long_irreg>0) {
            cat("Attention,",nb_long_irreg,"codes ne faisaient pas 5 caractères (malgré l'ajout de zéros initiaux). Le résultat ci-dessous n'est pas garanti.\n") ;
        }
        selecteur <- which(nchar(geocode)==4) ;
        geocode[selecteur] <- paste("0",geocode[selecteur], sep="" ) ;
        dept <- substr(geocode, 1, 2) ;
        
        # substitution des DOM
        dom_correspondance <- c("971"="9A", "972"="9B", "973"="9C", "974"="9D", "975"="9E", "976"="9F" ) ;
        for(i in seq_len(length(dom_correspondance))) {
            dept[dept == dom_correspondance[i]] <- names(dom_correspondance)[i] ;
        }
        # substitution des TOM
        tom <- paste0("9", LETTERS[7:26]) ;
        dept[dept %in% tom] <- "98" ;
        # Code 99 for foreign countries
        
    } else {
        
        # depuis un numéro de département... à vérifier a minima
        nb_long_irreg <- sum(nchar(geocode)!=1 & nchar(geocode)!=2) ;
        if(nb_long_irreg>0) {
            cat("Attention,",nb_long_irreg,"codes ne faisaient pas 2 caractères (malgré l'ajout de zéros initiaux). Le résultat ci-dessous n'est pas garanti.\n") ;
        }
        selecteur <- which(nchar(geocode)==1) ;
        geocode[selecteur] <- paste("0",geocode[selecteur], sep="" ) ;
        dept <- substr(geocode, 1, 3) ; # on laisse 3 à cause des DOM TOM
        
    }
    return(dept)
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




#############################################################################################
#############################################################################################
#############################################################################################

# UNIVARIEE 

desc_quanti_cont <- function(vector, name="Variable", mean_ci=TRUE, alpha=0.05, def_breaks = "Sturges", plot_boxplot = F, graph=TRUE, density=TRUE, ...) {
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
    print(rbind(as.matrix(summary(vector)), Sd = sd(vector))) ;
    
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
        boxplot(vector, main = name, ylab = "")
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

# Format : "aaaa-mm-jj" ou "jj/mm/aaaa" ou un vrai format R
desc_date <- function(vector, name="Variable", format="%Y-%m-%d", breaks="weeks", ...) {
    # Si le format est précisé à la Française, on peut rattrapper le coup.
    if( format=="aaaa-mm-jj") {
        format="%Y-%m-%d" ;
    } else if( format=="jj/mm/aaaa") {
        format="%d/%m/%Y" ;
    }
    vector <- as.Date(vector, format=format) ;
    
    cat(name,"\n") ;
    cat("Effectif analysé :", length(na.omit(vector)),"\n") ;
    cat("------------------------------------------------------------------------------------\n") ;
    if( sum(is.na(vector))==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.\n") ;
        vector <- vector[!is.na(vector)] ;
    }
    print(rbind(as.matrix(summary(vector)), Sd = sd(vector))) ;
    
    hist(vector, breaks=breaks, plot=TRUE, freq=FALSE, main=name, xlab=name, ylab ="Density", 
         start.on.monday = TRUE, col = "#AAAAFF")
    temp <- na.omit(as.numeric(vector)) ; 
    lines(density(temp), col = "red") ;
}

# Permet de tracer un histogramme de la distribution d'une variable de type heure
desc_heure <- function(vector, name="Variable", format="%H:%M", breaks="hours", ...) {
    # Si le format est précisé à la Française, on peut rattrapper le coup.
    if( format=="hh:mm") {
        format="%H:%M" ;
    }
    vector <- strptime(vector, format=format) ;
    
    cat(name,"\n") ;
    if( sum(is.na(vector))==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.\n") ;
        vector <- vector[!is.na(vector)] ;
    }
    cat("Effectif analysé :", length(na.omit(vector)),"\n") ;
    cat("------------------------------------------------------------------------------------\n") ;
    print(rbind(as.matrix(summary(vector)), Sd = sd(vector))) ;
    
    hist(vector, breaks=breaks, plot=TRUE, freq=FALSE, main=name, xlab=name, ylab ="Density", 
         start.on.monday = TRUE, col = "#AAAAFF")
    temp <- na.omit(as.numeric(vector)) ; 
    lines(density(temp), col = "red") ;
}


# valeurs de "sort" : "alpha", "croissant", "decroissant" ||| xlim = limite de x pour le graph
# sum = booléen. Afficher le summary() ?
desc_quanti_disc <- function(vector, name="Variable", mean_ci=TRUE, table=TRUE, Sum = T, sort="alpha", xlim=NULL, ...) {
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
    if (Sum) {print(rbind(as.matrix(summary(vector)), Sd = sd(vector)))} ;
    
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
        cat("\nLes calculs des IC95% sont réalisés à partir de la loi binomiale\n")
        cat("\nModalite\tEffectif\tProportion\tIC95%\n") ;
        for( une_modalite in modalites) {
            temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE) ;
            cat(une_modalite,"\t",table(vector)[une_modalite],"\t",temp[1],"\t[",temp[2],";",temp[3],"]\n") ;
        }
    }
    if( mean_ci ) {
        sd <- sd(vector) ;
        mean <- mean(vector) ;
        n <- length(vector) ;
        cat( "Moyenne et intervalle de confiance à 95% :",  
             round(mean,2),"[",round(mean-1.96*sd/sqrt(n),2),";",round(mean+1.96*sd/sqrt(n),2),"].\n",
             "Calcul des IC95% à partir du théorème central limite") ;
    }
    plot(table(vector)/length(vector), xlab=name, ylab="proportion", col="cornflowerblue", xlim=xlim) ;
    
}


desc_binaire <- function(vector, name="Variable", table=TRUE, ...) {
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
    if( table ) {
        temp <- as.data.frame(table(vector)) ;
        modalites <- temp[,1] ;
        cat("\nModalite\tEffectif\tProportion\tIC95%\n") ;
        for( une_modalite in modalites) {
            temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE) ;
            cat(une_modalite,"\t",table(vector)[une_modalite],"\t",temp[1],"\t[",temp[2],";",temp[3],"]\n") ;
            knitr::kable
        }
        cat("\nCalcul des IC95% à l'aide d'une loi binomiale")
    }
    pie(table(vector)/length(vector), main=name, col=c("white", "cornflowerblue")) ;
}



# limit_chart permet de limiter le nombre de sorties dans le graphique
# tronque_lib_chart permet de limiter la taille des étiquettes (code et libellé compris)
# sort >>>> alpha / decroissant / croissant
# tronque_lib_tab 0 par défaut (affiche tout) sinon n'affiche que les x premiers
# Pour afficher une modalité avec un effectif nul, il faut que l'argument vector soit un facteur contenant l'ensemble des levels à afficher
desc_quali <- function(vector, name="Variable", table=TRUE, sort="alpha", limit_chart=Inf, tronque_lib_chart=20, tronque_lib_tab=0, ...) {
    cat(name,"\n") ;
    if( length(unique(na.omit( vector ))) <2 ) {
        cat("Cette colonne comporte au plus une valeur et ne sera pas analysée\n") ;
        return( TRUE ) ;
    }
    if( sum(is.na(vector))==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", sum(is.na(vector)),"soit",100*mean(is.na(vector)),"%.\n") ;
        vector <- vector[!is.na(vector)] ;
    }
    cat("Effectif analysé :", length(na.omit(vector)),"\n") ;
    cat("------------------------------------------------------------------------------------\n") ;
    
    temp <- as.data.frame(table(vector)) ;
    modalites <- temp[,1] ;
    if( sort=="alpha") {
        modalites <- modalites[order(modalites)] ;
    } else if( sort=="croissant") {
        modalites <- modalites[order(temp[,2])] ;
    } else if( sort=="decroissant") {
        modalites <- modalites[order(0-temp[,2])] ;
    }
    # tableau de contingence
    if( table ) {
        cat("\nModalite\tEffectif\tProportion\tIC95%\n") ;
        if ( tronque_lib_tab == 0 ) {
            for( une_modalite in modalites) {
                temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE) ;      
                cat(une_modalite,"\t",table(vector)[une_modalite],"\t",temp[1],"\t[",temp[2],";",temp[3],"]\n") ;
            }
        } else {
            for( une_modalite in modalites[1:tronque_lib_tab]) {
                temp <- confint_prop_binom(vecteur=(vector==une_modalite), pourcent=TRUE) ;      
                cat(une_modalite,"\t",table(vector)[une_modalite],"\t",temp[1],"\t[",temp[2],";",temp[3],"]\n") ;
            }
        }
        cat("Calcul des IC95% à l'aide d'une loi binomiale")
    }
    # graphique maintenant
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


#############################################################################################
#############################################################################################
#############################################################################################

# MULTIVARIEE


# dans output_prop_survie, préciser éventuellement des dates pour afficher les % de survie
desc_survie <- function(vector_evt, vector_delai, name="Variable", xmax=NULL, ymin=0, output_prop_survie=c()) {
    cat(name,"\n") ;
    vector_evt <- as.numeric(vector_evt) ;
    vector_delai <- as.numeric(vector_delai) ;
    manquant <- (is.na(vector_evt) + is.na(vector_delai))>0 ;
    if( sum(manquant)==0 ) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", sum(manquant),"soit",100*mean(manquant),"%.\n") ;
        vector_evt <- vector_evt[!manquant] ;
        vector_delai <- vector_delai[!manquant] ;
    }
    effet_surv <- Surv(event=vector_evt, time=vector_delai, type="right") ;
    obj_survfit2 <- survfit( effet_surv~1, conf.int=TRUE) ;
    plot( obj_survfit2, xmax=xmax, ymin=ymin, main=name ) ;
    cat("Informations sur la survie (dont la médiane) :\n")
    print( obj_survfit2) ;
    # on peut ajouter des % de survie à telles dates...
    if( length(output_prop_survie)>0) {
        cat("Survies estimées par la méthode de Kaplan-Meier :\n") ;
        # en général on a peu d'évenements donc il est intéressant de supprimer les censures
        y <- c() ;    y_upper <- c() ;    y_lower <- c() ;
        x <- c() ;
        for( i in length(obj_survfit2$time) : 1 ) {
            if(obj_survfit2$n.event[i]>=1) {
                x <- c(x, obj_survfit2$time[i]) ;
                y <- c(y, obj_survfit2$surv[i]) ;
                y_upper <- c(y_upper, obj_survfit2$upper[i]) ;
                y_lower <- c(y_lower, obj_survfit2$lower[i]) ;
            }
        }
        # maintenant on parcourt
        for( a_time in output_prop_survie) {
            if( a_time>max(obj_survfit2$time)) {
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


# analyse des résidus d'une régression linéaire multiple
# regression est l'objet de régression (ou de step), et "variables" liste les noms des variables X à traiter
# Y est toujours traité en premier, variables_x peut donc rester vide
# le dataframe sera récupéré dans l'objet régressino
reg_lin_analyse_residus <- function(regression, variables_x = names(regression$model)[-1]) {
    # histogramme des résidus
    hist(regression$residuals, col="#AAAAFF", freq=FALSE, xlab="Residus", main="Histogramme des résidus") ;
    # résidus en fonction de Y
    y <- row.names(attr(regression$terms, which = "factors"))[1]
    plot(x=regression$model[,y], y=regression$residuals, xlab="Y", ylab="Residus") ;
    
    for(un_x in variables_x) {
        # résidus en fonction de chaque X
        if (is.integer(regression$model[,un_x]) | is.numeric(regression$model[,un_x]) | is.factor(regression$model[,un_x])) {
            plot(x=regression$model[,un_x], y=regression$residuals, xlab=un_x, ylab="Residus") ;
        } else if (is.character(regression$model[,un_x])) {
            plot(x=factor(regression$model[,un_x]), y=regression$residuals, xlab=un_x, ylab="Residus") ;
        } 
    }
}


# distance de Cook d'une régression linéaire
# nom_colonne_id est le nom de la colonne à utiliser pour retrouver l'identifiant. 
# à défaut, on prendra le numéro de la ligne.
reg_lin_dist_cook <- function(regression, nom_colonne_id=NULL) {
    temp <- cooks.distance(regression) ;
    if( is.null(nom_colonne_id)) {
        numeros <- 1:nrow(regression$model) ;
    } else {
        numeros <- regression$model[,nom_colonne_id] ;
    }
    plot(temp, type="h", xlab="numero individu", ylab="distance de Cook", xlim=c(1,length(temp))) ;
    temp2 <- data.frame(
        index=1:nrow(regression$model),
        id=numeros,
        cook=temp
    ) ;
    temp2 <- temp2[order(-temp2$cook),] ;
    cat("Les cinq pires individus quant à la distance de Cook :\n") ;
    print(temp2[1:5,]) ;
    text(labels=temp2[1:5, "id"], y=temp2[1:5, "cook"], x=temp2[1:5, "index"], cex=0.75, col="blue") ;
}


# coefficient de détermination d'une régression linéaire multiple
reg_lin_R2 <- function(regression) {
    # méthode universelle
    coefficient <- (summary(regression)$null.deviance - summary(regression)$deviance)/summary(regression)$null.deviance  ;
    # méthode simple pour modèle linéaire avec intercept, pour information
    # coefficient <- cor(regression$fitted.values, regression$y)^2
    cat("Coefficient de détermination :", coefficient, "\n") ;
    # return( coefficient) ;
}

## Graphique des OR pour régression logistique
reg_log_graphe_or <- function(obj_reglog, main="Logistic regression", xlab="Odds ratios", xlim=NULL) {
    temp <- summary(obj_reglog)$coefficients ; # 1 estimate, 2 DS, 4 p val
    label <- row.names(temp) ;
    prop <- exp(temp[,1]) ;
    prop_ib <- exp(temp[,1]-1.96*temp[,2]) ;
    prop_sb <- exp(temp[,1]+1.96*temp[,2]) ;
    pval <- temp[,4] ;
    signif <- pval<0.05 ; # noter que ce test est meilleur que l'IC construit.
    if(is.null(xlim)) {
        # on prend les limites sans l'intercept, qui fausse tout
        prop_ib2 <- prop_ib[label!="(Intercept)"] ;
        prop_sb2 <- prop_sb[label!="(Intercept)"] ;
        xlim <- c( min(c(0, prop_ib2))-0.1  , max(prop_sb2)+0.1) ;
    }
    graphe_ratio(label=label, nb=prop, nb_ib=prop_ib, nb_sb=prop_sb, xlim=xlim, main=main, xlab=xlab, ref=1, signif=signif, log=TRUE) ;
}

## Fonction pour vérifier la loglinéarité des variables (Coralie Delettrez)
reg_log_verif_loglineaire <- function(x, y, xname="Variable quantitative", yname="Variable binomiale (var à expliquer du modele de reglog)", titre="titre") {
    y <- as.numeric(as.character(y))
    propll <- aggregate(y,list(x),mean,na.rm=TRUE)
    names(propll) <- c("x","prop")
    ll<-cbind(merge(cbind(x,y),propll, by="x",all=TRUE),temp_deciles=as.integer(rownames(data)))
    ll$temp_deciles1<-as.numeric(cut(ll$temp_deciles,breaks = c(seq(from = 1, to = length(ll$temp_deciles), by = length(ll$temp_deciles)/10),length(ll$temp_deciles)+1),include.lowest = TRUE))
    verifll<-aggregate(ll,list(ll$temp_deciles1),mean,na.rm=TRUE)
    logit<-log((1-verifll$prop)/verifll$prop)
    droite <- lm(log((1-prop)/prop) ~ x, data = verifll)
    plot(verifll$x,log((1-verifll$prop)/verifll$prop),xlab=xname,ylab=yname,main=titre)
    abline(droite,col="red")}

## Fonction de validation d'un modèle logistique
reg_log_validation <- function(mod, desc_eve = T, roc = T, calib = T, quant.hoslem = 10) {
    if(!require("ResourceSelection")) {install.packages("ResourceSelection")}
    # 1. Nombre d'évènements suffisants
    if (desc_eve) {
        nb_eve <- min(table(mod$y))
        nb_var <- length(names(mod$model)) - 1
        cat("1. Nombre d'évènements suffisant\n")
        cat(paste("Le nombre d'évènements (ou de non-évènement) total est de :", nb_eve, "\n"))
        cat(paste("Il devrait être supérieur à", 5*nb_var, "-", 10*nb_var, "évènements.\n\n\n"))
        
        min_gpe <- integer()
        for (i in names(mod$model)[-1]) {
            min_gpe <- min(min_gpe, table(mod$y, mod$model[,i]))
        }
        cat(paste("Le plus petit groupe du croisement de Y avec une variable est de :", min_gpe, "(doit être supérieur à 5).\n\n\n"))
        
    }
    
    # 2. Qualité de prédiction (Courbe ROC, pseudo R2 de Mc Fadden)
    if (roc) {
        cat("2. Qualité de prédiction\n")
        plot(jitter(mod$fitted.values, factor = 2000), mod$y, main = "Y observé en fonction de Y prédit", xlab = "Prédit", ylab = "Observé")
        cat("La courbe ROC permet d'étudier le pouvoir discriminant du modèle, ainsi que de déterminer le meilleur seuil.\n")
        courbe_roc(mod$y, mod$fitted.values)
        
    }
    
    # 3. Calibration du modèle : test de Hosmer-Lemeshow
    if (calib) {
        test <- ResourceSelection::hoslem.test(mod$y, mod$fitted.values, g = quant.hoslem)
        cat("\n\n3. Calibration du modèle (test du Khi2 d'Hosmer-Lemeshow)\n")
        cat("La probabilité prédite est-elle représentative de la probabilité observée ?\n")
        if (test$statistic < 0.05) {
            cat(paste("Rejet de H0 : le modèle ne reflète pas la vrai probabilité ( p = ",round(test$statistic, 2), "), il n'est pas calibré.\n\n"))
        } else {
            cat(paste("Non rejet de H0 : le modèle reflète la vrai probabilité ( p = ",round(test$statistic, 2), "), il est bien calibré.\n\n"))
        }
    }
}

##### Regression logistique avec selection de variable (imputation multiple sur les NA) (Coralie Delettrez)
# y est la variable à expliquer
# vecteur correspond aux variables à expliquer "c("var_1","var_2",...,"var_n")
# La selection de variable se fait en backward

RegLogLin_IM <- function(imp,y,vecteur)  {
 
  liste_fit <- list()                              #Initialisation de la liste qui contiendra les différentes Régressions logisitiques          
  pc <- c()
  pc$pvalue <- 1                                   #Initialisation de la boucle
  i <- 1                                              
 
  while (pc$pvalue > 0.2) {                        #La sélection de variable s'arrête quand le test du rapport de vraisemblance est inf à 0.2
    vecteur <- as.matrix(vecteur)
    fit <- with(imp,glm((fmla <- as.formula(paste(y, "~", paste(vecteur, collapse= "+")))), family="binomial"))   #Réalisation de la Régression logistique sur les variables sélectionnées
    liste_fit[[i]] <- fit                          #Stockage des différentes Régressions logistiques
    aa <- summary(pool(fit))[-1,]                  #Stockage des p_valeurs de toutes les variables restantes
   
    if(length(vecteur)>1){                         #Critère d'arrêt s'il n'y a plus de variable
      a_sortir <- which.max(aa[,5])
      vecteur <- vecteur[-a_sortir] }             #Rejet de la variable ayant la plus grande p_value
   

    #Pour la réalisation du test de rapport de vraisemblance : on stocke dans comp0 toutes les variables retenues
    #et dans comp1 *dans le même ordre* (pour pouvoir utiliser pool.compare) toutes les variables retenues moins celle ayant la plus grande p-value
    comp0 <- with(imp,glm((fmla1 <- as.formula(paste(y, "~", paste(vecteur, collapse= "+")))), family="binomial"))
   
    if (substring(names(a_sortir),nchar(names(a_sortir)),nchar(names(a_sortir)))==2) {
      comp1 <- with(imp,glm(fmla2 <- as.formula(paste(y, "~", paste(c(vecteur,substring(names(a_sortir),1,nchar(names(a_sortir))-1)), collapse= "+"))), family="binomial"))
      } else {
      comp1 <- with(imp,glm(fmla2 <- as.formula(paste(y, "~", paste(c(vecteur,names(a_sortir)), collapse= "+"))), family="binomial"))}
   
   
    pc <- pool.compare(comp1, comp0, method='likelihood', data=imp)
    i <- i+1
  }
  return(summary(pool(comp1)))
}



## Graphique des HR d'un modèle de cox
reg_cox_graphe_hr <- function(obj_cox, main="Survival", xlab="Hazard ratios", xlim=NULL) {
    temp_confint <- summary(obj_cox)$conf.int ;
    temp_pval <- summary(obj_cox)$coefficients ;
    label <- row.names(temp_confint) ;
    prop <- temp_confint[,1] ;
    prop_ib <- temp_confint[,3] ;
    prop_sb <- temp_confint[,4] ;
    pval <- temp_pval[,5] ;
    signif <- pval<0.05 ; # noter que ce test est meilleur que l'IC construit.
    if(is.null(xlim)) {
        # pas d'intercept ici
        xlim <- c( min(c(0, prop_ib))-0.1  , max(prop_sb)+0.1) ;
    }
    graphe_ratio(label=label, nb=prop, nb_ib=prop_ib, nb_sb=prop_sb, xlim=xlim, main=main, xlab=xlab, ref=1, signif=signif, log=TRUE) ;
}


# Cette fonction est utilisée par les fonctions reg_log_graphe_or, reg_cox_graphe_hr et reg_lin_graphe_coeff.
# Ne cherchez pas à l'invoquer directement.
graphe_ratio <- function(label, nb, nb_ib, nb_sb, xlim, main, xlab, ref=-999, signif=NULL, log=FALSE) {
    if(!require("lattice")) install.packages("lattice", repos="http://cran.us.r-project.org") ;
    library(lattice) ;
    label2 <- as.factor(label) ;
    # label2 <- reorder(label2, nb) ;
    if( is.null(signif)) {
        signif <- rep(FALSE, length(nb)) ;
    }
    #xlim <- c(  min(0, nb_ib) , max(1, nb_sb) )
    plot(dotplot(label ~ nb, xlim=xlim, panel=function(x,y) {    
        if( FALSE ) {
            xscale.components.default(lim=xlim, scales=list(x=list(log=2))) ;
            scales=list(x=list(log=2)) ;
        }
        panel.abline(v=ref,col='#CCCCCC', lty=2 )
        panel.xyplot(x,y,pch=16,cex=1,col=ifelse(signif,'#990000','navy'))
        panel.segments(nb_ib,as.numeric(y),nb_sb,as.numeric(y),lty=1,col=ifelse(signif,'#990000','navy'))
    }, xlab=xlab, main=main)) ;
}


reg_lin_graphe_coeff <- function(regression, main="Regression", xlab="Coefficients", xlim=NULL) {
    temp <- summary(regression)$coefficients ;
    label <- row.names(temp) ;
    prop <- temp[,1] ;
    prop_ib <- temp[,1]-1.96*temp[,2] ;
    prop_sb <- temp[,1]+1.96*temp[,2] ;
    signif <- temp[,4]<0.05 ;
    if(is.null(xlim)) {
        # on prend les limites sans l'intercept, qui fausse tout
        prop_ib2 <- prop_ib[label!="(Intercept)"] ;
        prop_sb2 <- prop_sb[label!="(Intercept)"] ;
        xlim <- c( min(c(0, prop_ib2))-0.1  , max(prop_sb2)+0.1) ;
    }
    graphe_ratio(label=label, nb=prop, nb_ib=prop_ib, nb_sb=prop_sb, main=main, xlab=xlab, ref=0, signif=signif, xlim=xlim, log=FALSE) ;
}


# Valeurs pour method : pearson (paramétrique), spearman (non paramétrique)
bivarie_quanti_quanti <- function(x, y, xname="Variable quantitative 1", yname="Variable quantitative 2", method="pearson", droite_reg=TRUE) {
    cat("Analyse bivariée : Méthode du coefficient de corrélation de ",method," sur deux variables quantitatives\n")
    coeff_r <- cor(x, y, method=method, use="complete.obs" ) ;
    coeff_r2 <- coeff_r^2 ;
    temp <- cor.test(x, y, method=method, alternative="two.sided") ;
    pval <- temp$p.value ;
    cat("Coefficient de correlation de", method, ": r=",coeff_r,", avec r²=", coeff_r2," (p=",pval,")\n") ;
    plot(x=x, y=y, type="p", pch=20, xlab=xname, ylab=yname) ;
    if( droite_reg ) {
        reg <- glm(y~x) ;
        b <- reg$coefficients[1] ;
        a <- reg$coefficients[2] ;
        cat("Equation droite : y =",a,"* x +",b,"\n") ;
        points_x <- c(min(na.omit(x)), max(na.omit(x))) ;
        points_y <- a*points_x + b ;
        lines(x=points_x, y=points_y, type="l", col="blue") ;
    }
}


# Valeurs pour method : chisq (semi-paramétrique), fisher (non paramétrique mais attention au nombre de modalités)
bivarie_quali_quali <- function(x, y, xname="Variable qualitative 1", yname="Variable qualitative 2", method="chisq", bin = FALSE, prop.table = TRUE, table = FALSE, mosaic = T,...) {
    nb_valide <- nrow(na.omit(data.frame(x,y)))
    nb_manquant <- length(x) - nb_valide
    cat(paste0(yname, " en fonction de ", xname, ".\n")) ;
    if (nb_manquant==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", nb_manquant,"soit",100*nb_manquant/length(x),"%.\n") ;
    }
    cat(paste0("Effectif analysé : ", nb_valide, ".\n")) ;
    cat("------------------------------------------------------------------------------------\n") ;
    
    if( method=="chisq") {
        cat("Analyse bivariée : Test du Chi 2 sur deux variables qualitatives\n")
        obj <- chisq.test(table(x,y)) ;
        print(obj) ;
        pval <- obj$"p.value" ;
    } else if( method=="fisher") {
        cat("Analyse bivariée : Test de Fisher sur deux variables qualitatives\n")
        obj <- fisher.test(table(x,y)) ;
        print(obj) ;
        pval <- obj$"p.value" ;
    } else {
        cat("Nom de méthode incorrect :", method,"\n") ;
        return("Erreur !") ;
    }
    cat("\n------------------------------------------------------------------------------------\n")
    cat("La p-value (petit p) de ce test = ",pval)
    cat("\n------------------------------------------------------------------------------------\n")
    
    if( bin ) {
        cat("\nProportions de 1 et effectifs de la variable binaire pour chaque modalité de l'autre variable")
        cat("\nModalite\tProportion de 1\t\tEffectif\n")
        for(une_modalite in sort(unique(na.omit(x)))) {
            n <- sum(!is.na(y[x == une_modalite]))
            cat(une_modalite,"\t\t",round(sum(y[x == une_modalite], na.rm = TRUE)/n*100,2),"%\t\t",sum(y[x == une_modalite], na.rm = TRUE) ,"\n")
        } 
    }
    if (prop.table) {
        prop_table <-round(100*prop.table(table(y,x),2),2)  
        tableau_perc <- kable(prop_table, caption=paste0("Proportions de y (", yname, ") par modalité de x"), row.names = T)
        print(tableau_perc)
    }
    if (table) {
        tableau_eff <- kable(table(y,x), caption="Tableau des effectifs", row.names = T)
        print(tableau_eff)
    }
    
    if (mosaic) {plot(table(x, y), xlab=xname, ylab=yname, main="Mosaicplot", col="cornflowerblue")} ;
}


# Valeurs pour method : 
# - student (2groupes, paramétrique) ; 
# - anova (plus de 2 groupes, paramétrique) ;
# - wilcoxon (Mann-Whitney, 2 groupes, non paramétrique) ;
# - kruskal (Kruskal-Wallis, plus de 2 groupes, non paramétrique)

bivarie_quali_quanti <- function(x, y, xname="Variable qualitative", yname="Variable quantitative", method="student", cond.app = F, graph = TRUE, ...) {
    if(!is.factor(x)) {x=as.factor(x)}
    nb_valide <- nrow(na.omit(data.frame(x,y)))
    nb_manquant <- length(x)-nb_valide
    cat(paste0(yname, " en fonction de ", xname, ".\n")) ;
    if (nb_manquant==0) {
        cat("Aucune valeur manquante.\n") ;
    } else {
        cat("Valeurs manquantes : n=", nb_manquant,"soit",100*nb_manquant/length(x),"%.\n") ;
    }
    cat(paste0("Effectif analysé : ", nb_valide, ".\n")) ;
    cat("------------------------------------------------------------------------------------\n") ;
    
    if (cond.app == T) {
        cat("Vérification des conditions d'application :\n")
        cat("1. Normalité de la distribution dans chaque groupe (méthode graphique). Voir graphiques.\n")
        var_min <- numeric()
        var_max <- numeric()
        for (une_modalite in unique(x)) {
            normPlot(y[x == une_modalite], une_modalite)
            ifelse(length(na.omit(y[x == une_modalite])) >= 30, cat("Note : l'Effectif de l'échantillon est supérieur ou égal à 30 (hors NA).\n"), cat(""))
        }
        var_min <- min(c(var_min, sd(y[x == une_modalite])^2, na.rm = T))
        var_max <- max(c(var_max, sd(y[x == une_modalite])^2, na.rm = T))
        cat("2. Egalité des variances.\n")
        cat(paste("Le rapport de la variance max / variance min = ", round(var_max/var_min, 2), "(doit être inférieur à 1,5)"))
        cat("\n------------------------------------------------------------------------------------\n")
    }
    
    if( method=="student") {
        cat("Analyse bivariée : Test t de Student d'une variable qualitative avec une variable quantitative\n")
        obj <- t.test(y~x, var.equal=FALSE);
        pval <- as.numeric(obj$p.value);
        print(obj) ;
    } else if( method=="anova") {
        cat("Analyse bivariée : Test d'analyse de variances (ANOVA) d'une variable qualitative avec une variable quantitative\n")
        obj <- aov(y~x);
        pval <- as.numeric(summary(obj)[[1]][["Pr(>F)"]][1])
        print(summary(obj)) ;
    } else if( method=="wilcoxon") {
        cat("Analyse bivariée : Test de Wilcoxon d'une variable qualitative avec une variable quantitative\n")
        obj <-wilcox.test(y~x);
        pval<-obj$p.value;
        print(obj) ;
    } else if( method=="kruskal") {
        cat("Analyse bivariée : Test de Kruskal-Wallis d'une variable qualitative avec une variable quantitative\n")
        obj <-kruskal.test(x=y, g=x, na.action=na.rm);
        pval<-obj$p.value;
        print(obj) ;
    } else {
        cat("Nom de méthode incorrect :", method,"\n") ;
        return("Erreur !") ;
    }
    cat("\n------------------------------------------------------------------------------------\n")
    cat("La p-value (petit p) de ce test = ",pval)
    cat("\n------------------------------------------------------------------------------------\n")
    
    if (method %in% c("student", "anova")) {
        cat("\nMoyennes avec IC95% pour chaque modalité de la variable qualitative\n")
        cat("\nModalite\tEffectif\tMoyenne\tIC95%\n", sep = "")
        for (une_modalite in sort(unique(na.omit(x)))) {
            quali <- une_modalite
            n <- sum(!is.na(y[x == quali]))
            mean_quanti <- mean(y[x == quali], na.rm = TRUE)
            mean_quanti <- round(mean_quanti,2)
            low_bound <- round(mean_quanti-1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
            upp_bound <- round(mean_quanti+1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
            cat(une_modalite,"\t", n, "\t", mean_quanti,"\t[",low_bound,";",upp_bound,"]\n", sep = "")
        }
    } else if (method %in% c("wilcoxon", "kruskal")) {
        cat("\nMédianes avec intervalle inter-quartile pour chaque modalité de la variable qualitative\n")
        cat("\nModalite\tEffectif\tMédiane\tIQ\n", sep = "")
        for (une_modalite in sort(unique(na.omit(x)))) {
            quali <- une_modalite
            n <- sum(!is.na(y[x == quali]))
            med_quanti <- median(y[x == quali], na.rm = TRUE)
            med_quanti <- round(med_quanti,2)
            low_bound <- round(quantile(y[x == quali], probs = 0.25, na.rm = TRUE), 2)
            upp_bound <- round(quantile(y[x == quali], probs = 0.75, na.rm = TRUE), 2)
            cat(une_modalite,"\t", n, "\t", med_quanti,"\t[",low_bound,";",upp_bound,"]\n", sep = "")
        }
        cat("\nMoyennes avec IC95% pour chaque modalité de la variable qualitative\n")
        cat("\nModalite\tEffectif\tMoyenne\tIC95%\n", sep = "")
        for (une_modalite in sort(unique(na.omit(x)))) {
            quali <- une_modalite
            n <- sum(!is.na(y[x == quali]))
            mean_quanti <- mean(y[x == quali], na.rm = TRUE)
            mean_quanti <- round(mean_quanti,2)
            low_bound <- round(mean_quanti-1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
            upp_bound <- round(mean_quanti+1.96*sd(y[x == quali], na.rm = TRUE)/sqrt(n),2)
            cat(une_modalite,"\t", n, "\t", mean_quanti,"\t[",low_bound,";",upp_bound,"]\n", sep = "")
        }
    }
    
    if ( graph ) {
        boxplot(y~x, xlab=xname, ylab=yname, col="cornflowerblue", ...) ;
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
acp <- function(dataframe) {
    if(!require("FactoMineR")) install.packages("FactoMineR", repos="http://cran.us.r-project.org") ;
    library(FactoMineR) ;
    
    nrow_av <- nrow(dataframe) ;
    dataframe <- na.omit(dataframe) ;
    nrow_ap <- nrow(dataframe) ;
    manquants <- nrow_av - nrow_ap ;
    if( manquants>0 ) {
        cat("Valeurs manquantes. Lignes ignorées pour l'ACP : n=", manquants,"soit",100*manquants/nrow_av,"%.\n") ;
    }
    # avec le paramètre graph=T, cette méthode nous tracerait aussi le plan factoriel des variables et celui des individus
    mon_acp <- PCA(dataframe, scale.unit=TRUE, ncp=3, graph=F) ;
    
    # pour CP1-CP2 et CP1-CP3, plan factoriel des variables et plan des individus
    for( y in 2:3) {
        ylab <- paste("CP",y, sep="") ;
        col <- "black" ;
        plot(x=mon_acp$var$cor[,1], y=mon_acp$var$cor[,y], xlab="CP1", ylab=ylab, pch=".", cex=6, xlim=c(-2,2), ylim=c(-1,1), bty="n", col="red") ;
        lines(x=c(-1,1), y=c(0,0)) ;
        lines(y=c(-1,1), x=c(0,0)) ;
        cercle_x <- seq(from=-1, to=1, by=0.02) ; cercle_y <- sqrt(1-cercle_x^2) ; lines(y=cercle_y, x=cercle_x) ; cercle_y <- 0-sqrt(1-cercle_x^2) ; lines(y=cercle_y, x=cercle_x) ;
        text(x=mon_acp$var$cor[,1], y=mon_acp$var$cor[,y], labels=row.names(mon_acp$var$cor), col="blue", pos=4, cex=0.75) ;
        plot(x=mon_acp$ind$coord[,1],y=mon_acp$ind$coord[,y], xlab="CP1", ylab=ylab, pch=".", cex=4, col=col)
    }
    # valeurs propres
    plot(x=1:nrow(mon_acp$eig), y=mon_acp$eig[,"percentage of variance"]/100, type="h", ylim=c(0,1), xlab="composante", ylab="part de variance expliquée") ;
    plot(x=1:nrow(mon_acp$eig), y=mon_acp$eig[,"cumulative percentage of variance"]/100, type="b", ylim=c(0,1), xlab="composante", ylab="part de variance expliquée cumulée") ;
    cat("Part de variance expliquée cumulative :\n") ;
    print(mon_acp$eig[,"cumulative percentage of variance"]) ;
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


#' @title Plot pyramide des ages
#' @description Dessine une pyramide des ages en fonction des ages et du sexe
#' @param age integer Ages en années
#' @param sexe boolean Sexe masculin = 1 sinon est féminin
#  package cowplot pour placer correctement les labels des axes (AGE) => avec fonction ggdraw
#' @import ggplot2 grid dplyr

plot_pyramide_ages2 <- function(age, sexe, largeur = 5) {
    if(!require("ggplot2")) install.packages("ggplot2", repos="http://cran.us.r-project.org") ;
    library(ggplot2) ;
    if(!require("grid")) install.packages("grid", repos="http://cran.us.r-project.org") ;
    library(grid) ;
    if(!require("dplyr")) install.packages("dplyr", repos="http://cran.us.r-project.org") ;
    library(dplyr) ;
    #if(!require("cowplot")) install.packages("cowplot", repos="http://cran.us.r-project.org") ;
    #library(cowplot) ;
    hommes <- sexe == 1
    # Arrondir par classes d'age
    ages <- cut(age, breaks=c(-Inf, seq(from = largeur, to = 100, by = largeur), Inf)) ;
    repartition_hommes <- as.data.frame(table(ages[hommes]))
    repartition_femmes <- as.data.frame(table(ages[!hommes]))
    repartition <- merge(x = repartition_hommes, y = repartition_femmes, by = "Var1")  
    repartition <- rename(repartition, age = Var1, hommes = Freq.x, femmes = Freq.y)   
    nombre_max <- max(repartition$hommes, repartition$femmes)
    gghommes <- 
        ggplot(repartition) +
        aes(x = age, y = hommes) +
        geom_bar(stat = "identity", fill = "cornflowerblue") +
        coord_flip(ylim = c(0, nombre_max)) +
        scale_y_reverse() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #         axis.text.y = element_blank(),
            #         axis.ticks.y = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "mm")
        ) +
        theme_bw() 
    ggfemmes <- 
        ggplot(repartition) +
        aes(x = age, y = femmes) +
        geom_bar(stat = "identity", fill = "lightpink") +
        coord_flip(ylim = c(0, nombre_max)) +
        theme(
            axis.title.x = element_text(size = 0),
            axis.title.y = element_text(size = 0),
            #         axis.text.y element_text(size = 0),
            #         axis.ticks.y element_text(size = 0),
            #         plot.margin = unit(c(1, -1, 1, 0), "mm")
            plot.margin = unit(c(0, 0, 0, 0), "mm")
        ) +
        theme_bw()
    grid.newpage()
    # cadre_gauche <- viewport(x = 0.4, width = 0.2, name = "cadre_gauche")
    # cadre_droit <- viewport(x = 0.6, width = 0.3, name = "cadre_droit")
    # pushViewport(cadre_gauche)
    pushViewport(
        viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = c(0.5, 0.5)))
    )
    print(gghommes, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(ggdraw(switch_axis_position(ggfemmes + theme_bw(), axis = 'y')), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    return(invisible(NULL))
}



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





#############################################################################################
#############################################################################################
#############################################################################################

# CARTOGRAPHIE

## CARTO (alex georges)

carte_de_france <- function(departements, name="Nombre par département") {
    
    if(!require("maptools")) install.packages("maptools", repos="http://cran.us.r-project.org") ;
    library(maptools) ;
    if(!require("RColorBrewer")) install.packages("RColorBrewer", repos="http://cran.us.r-project.org") ;
    library(RColorBrewer) ;
    if(!require("classInt")) install.packages("classInt", repos="http://cran.us.r-project.org") ;
    library(classInt) ;
    if(!require("rgeos")) install.packages("rgeos", repos="http://cran.us.r-project.org") ;
    library(rgeos) ;
    
    map_france <- readShapeSpatial("geofla/DEPARTEMENT/DEPARTEMENT.SHP")
    
    dpt <- as.data.frame(table(departements))
    names(dpt) <- c("dpt","nombre")
    #Jointure entre le fond de carte et les données
    #map_france@data <- merge(map_france@data, dpt, by.x="CODE_DEPT", by.y="dpt", all.x=TRUE) ### BUG
    map_france@data <- data.frame(map_france@data, dpt[match(map_france@data[,"CODE_DEPT"], dpt[,"dpt"]),])
    
    #choix du nombre de couleurs/classes
    nclr <- 6
    #palette clr
    plotclr <- brewer.pal(nclr, "Greys")
    #discrétisation 
    distr <- classIntervals(map_france$nombre,nclr,style="quantile")$brks
    
    #attribution des couleurs aux régions
    colMap <- plotclr[(findInterval(map_france$nombre,distr,all.inside=TRUE))]
    #Affichage de la carte
    plot(map_france, col=colMap)
    #affichage de la légende
    legend("bottomleft", legend=myLeg(distr,2), fill = plotclr, cex = 0.6, bty = "n")
    #Titre
    title(main=name)  
}

## Fonction pour légender (arrondir les classes) sur la fonction de cartographie
myLeg <- function (vec, arrond) {  
    x <- vec  
    lx <- length(x)  
    if (lx < 3)    
        stop("pas suffisamment de classes")  
    res <- character(lx - 1)  
    res  
    for (i in 1:(lx - 1))    
    {res[i] <- paste(round(x[i],arrond), round(x[i + 1],arrond),sep=" - ")
    }  
    res  
}



## Fonction standardisation (Thibaut B) (Work in progress)

standardisation <- function(pop_etude, periode, pop_ref="france_entiere", sexe = T, age = T, departement = T, type = "directe") {
    
    stopifnot(is.data.frame(pop_etude), is.logical(sexe), is.logical(age))
    
    if (pop_ref == "france_entiere") {
        load("ressources/insee_population.RData")
    } else if (pop_ref == "france") {
        load("ressources/insee_population.RData", envir=globalenv() )
        insee_population <- insee_population[!insee_population$departement %in% c("971", "972", "973", "974", "975", "976"),]
    } else if (pop_ref == "monde") {
        
    } else {
        stopifnot(is.data.frame(pop_ref))
    }
    
    list_departement <- unique(pop_etude$departement)
    # Controles et listing des variables + exclusion des départements qui ne sont pas présents dans pop_etude
    list_var_join <- character()
    if (age) {
        list_var_join <- c(list_var_join, "age_num")
        stopifnot(sum(grepl(pattern = "age", names(pop_etude))) == 1)
    }
    if (sexe) {
        list_var_join <- c(list_var_join, "sexe")
        stopifnot(sum(grepl(pattern = "sexe", names(pop_etude))) == 1)
        insee_population <- insee_population[insee_population$sexe != "ensemble",]
    }
    if (departement) {
        list_var_join <- c(list_var_join, "departement")
        stopifnot(sum(grepl(pattern = "departement", names(pop_etude))) == 1)
    }
    stopifnot(length(!names(pop_etude) %in% list_var_join)>1)
    names(pop_etude)[!names(pop_etude) %in% list_var_join] <- "var"
    
    # filtre sur la période
    insee_population <- insee_population[insee_population$annee %in% periode,]
    population_ref_totale <- sum(insee_population$population[insee_population$departement %in% unique(pop_etude$departement)], na.rm = T) / length(periode)
    
    if (length(periode)>1) {
        insee_population <- aggregate(formula = population~age_num+departement+sexe , FUN = sum ,data = insee_population)
    } else {
        insee_population <- aggregate(formula = population~age_num+departement+sexe , FUN = sum ,data = insee_population)
    }
    
    # Jointure des deux tables
    if (is.null(list_var_join)) {
        warning("Pas de standardisation demandée !")
    } else {
        temp_pop_etude <- merge(x = pop_etude, y = insee_population, by = list_var_join, all = F, all.y = T, sort = F)
        temp_pop_etude$var[is.na(temp_pop_etude$var)] <- 0
    }
    
    # Calcul des taux bruts
    temp_pop_etude$tx_brut <- temp_pop_etude$var / temp_pop_etude$population
    
    # Calcul des taux standardisés
    if (type == "directe") {
        temp_pop_ref <- aggregate(population~sexe+age_num, data= insee_population, FUN = sum)
        names(temp_pop_ref) <- c("sexe" ,"age_num", "pop_ref")
        temp_pop_ref$pop_ref <- temp_pop_ref$pop_ref/length(periode)
        temp_result <- merge(temp_pop_etude, temp_pop_ref, by = c("sexe", "age_num"), all.x = T, sort = F)
        temp_result$nb_theorique <- temp_result$tx_brut * temp_result$pop_ref
        
        temp_result_ts <- aggregate(nb_theorique~departement, data = temp_result, FUN = sum)
        temp_result_ts$tx_std <- temp_result_ts$nb_theorique/population_ref_totale
        
    } else if (type == "indirecte") {
        
    }
    
    return(list(taux_stand = temp_result_ts[temp_result_ts$departement %in% list_departement,c("departement", "tx_std")], taux_bruts = temp_result[temp_result$departement %in% list_departement,]))
}



#############################################################################################
#############################################################################################
#############################################################################################

# DATA VISUALISATION

# fonction de réalisation de graphs avant-après
before_after <- function(x, y, xtxt, ytxt) {
    par(pin=c(1.5,4)) # 1,5 x 4 pouces 
    plot(
        rep(1,NROW(x)),
        x,
        type="n",
        xlim = c(1,2),
        ylim = range(x,y), # Affichage y du plus petit au plus grand x ou y 
        xlab = xtxt, # Etiquette x 
        ylab = ytxt,
        axes = F 
    )
    points(
        rep(2,NROW(y)),
        y,
        type="n"
    )
    axis(2)
    axis(1,labels = c("Avant","Après"), at = 1:2)
    segments(rep(1,NROW(x)), x, rep(2,NROW(y)), y, lty=1, lwd=1, col="blue")
    segments(rep(1,1), c(median(x)), rep(2,1), c(median(y)), lty=2, lwd=2, col="black")
}





#############################################################################################
#############################################################################################
#############################################################################################

# AUTRES

# Cette fonction est utilisée pour sauvegarder sur disque et charger en mémoire certaines ressources, comme les départements ou libellés CIM10 
# MAJ : 24/10/2016 Ajout de la population INSEE de 2006 à 2014 (mettre à jour la période en fonction de la diffusion des données sur le site)
charge_ressource <- function(ressource) {
    if(!require("data.table")) install.packages("data.table", repos="http://cran.us.r-project.org") ;
    library(data.table) ;
    
    if( ressource=="insee_noms_departements") {
        
        # source_dept <- "http://www.insee.fr/fr/methodes/nomenclatures/cog/telechargement/2015/txt/depts2015.txt" ;
        fichier  <- file.path("ressources", "insee_noms_departements.txt") ;
        fichier2  <- file.path("ressources", "insee_noms_departements.RData") ;
        if (!file.exists(fichier2) ) {
            nomenclature_dept <- fread(input=fichier, colClasses = "character", sep="\t", header=TRUE) ;
            setnames(nomenclature_dept,"DEP","code") ;
            setnames(nomenclature_dept,"NCCENR","libelle") ;
            nomenclature_dept <- nomenclature_dept[,list(code,libelle)]
            setkey(x=nomenclature_dept, code) ;
            insee_noms_departements <<- nomenclature_dept
            save(insee_noms_departements, file=fichier2 ) ;
        } else if( !exists("insee_noms_departements")) {
            load(file=fichier2, envir=globalenv() ) ;
        }
        # ***********************************************************************
        
    } else if( ressource=="insee_population") {
        
        # source_dept <- "http://www.insee.fr/fr/methodes/nomenclatures/cog/telechargement/2015/txt/depts2015.txt" ;
        lien <- "http://www.insee.fr/fr/themes/tableau_local_tsv.asp?ref_id=POP1B&nivgeo=DEP&codgeo=01&millesime=2013&niveau=2"
        fichier  <- file.path("ressources", "insee_population.RData") ;
        if (!file.exists(fichier) ) {
            
            if (!dir.exists(file.path("ressources", "Insee"))) dir.create(file.path("ressources", "Insee"))
            # Telechargements des fichiers population
            for (i in 2006:2014) {
                lien_temp <- sub(pattern = "(&millesime=)[0-9]{4}", replacement = paste0("\\1", i), x = lien)
                for (j in c(01:19, 21:95, "2A", "2B", "971", "972", "973", "974")) {
                    if (nchar(j)==1) j=paste0("0",j)
                    lien_temp <- sub(pattern = "(&codgeo=)[0-9]{1,3}[A-B]?", replacement = paste0("\\1", j), x = lien_temp)
                    nom_fichier <- paste("POP1B", i, j, sep = "_")
                    print(nom_fichier)
                    print(lien_temp)
                    download.file(lien_temp, destfile =  file.path("ressources", "Insee", paste0(nom_fichier, ".csv")), quiet = T)
                }
            }
            # fusion des fichiers
            for (i in list.files(chemin, pattern = ".csv$")) {
                nom_objet <- sub(pattern = ".csv", replacement = "", i)
                temp_obj <- scan(file = file.path(chemin, i), what = "character", sep = "\r", skip = 10)
                temp_obj <- read.table(text = temp_obj, header = F, stringsAsFactors = F, sep = "\t", nrows = 101, blank.lines.skip = T, fill = T)
                temp_obj <- temp_obj[,-5]
                names(temp_obj) <- c("age", "hommes", "femmes", "ensemble")
                assign(x = nom_objet, value = temp_obj)
            }
            
            insee_population <- data.frame()
            for (i in ls(pattern = "POP1B")) {
                temp_annee <- sub(pattern = ".+?_([0-9]{4})_.+", replacement = "\\1", i)
                temp_departement <- sub(pattern = ".+?_([0-9]{1,3}[A-B]?)$", replacement = "\\1", i)
                
                temp_obj <- get(i)
                temp_obj <- cbind(annee = rep(temp_annee, nrow(temp_obj)),
                                  departement = as.character(rep(temp_departement, nrow(temp_obj))),
                                  temp_obj)
                
                insee_population <- rbind(insee_population, temp_obj)
            }
            
            insee_population$age_num <- 0
            insee_population$age_num[nchar(insee_population$age) %in% 4:5] <- substring(insee_population$age[nchar(insee_population$age) %in% 4:5], 1, 1)
            insee_population$age_num[nchar(insee_population$age) == 6] <- substring(insee_population$age[nchar(insee_population$age) == 6], 1, 2)
            insee_population$age_num[nchar(insee_population$age) == 12 ] <- 0
            insee_population$age_num[nchar(insee_population$age) == 15] <- substring(insee_population$age[nchar(insee_population$age) == 15], 1, 3)
            
            insee_population$age_num <- as.integer(insee_population$age_num)
            insee_population$annee <- as.integer(insee_population$annee)
            
            insee_population <- reshape(insee_population, varying = c("hommes", "femmes", "ensemble"), v.names = "population", direction = "long", times = c("hommes", "femmes", "ensemble"), timevar = "sexe")
            row.names(insee_population) <- 1:nrow(insee_population)
            insee_population <- insee_population[,-7]
            
            setkey(x=as.data.table(insee_population), departement, age_num, sexe) ;
            insee_population <<- insee_population
            save(insee_population, file=fichier ) ;
        } 
        if( !exists("insee_population")) {
            load(file=fichier, envir=globalenv() ) ;
        }
        # ***********************************************************************
        
    }else if( ressource=="ign_fonds_cartes"  ) {
        
        if(!require("maptools")) install.packages("maptools", repos="http://cran.us.r-project.org") ;
        library(maptools) ;
        gpclibPermit() ; # pour débloquer certaines fonctionnalités
        
        fichier2  <- file.path("ressources", "ign_fonds_cartes.RData") ;
        if ( !file.exists(fichier2) ) {
            map_france_metropolitaine <<- readShapeSpatial("./ressources/ign_fonds_cartes/GEOFLA_2-0_SHP_LAMB93_FR-ED141_DEPARTEMENT") # métropole
            map_guadeloupe            <<- readShapeSpatial("./ressources/ign_fonds_cartes/GEOFLA_2-0_SHP_UTM20W84GUAD_D971-ED141_DEPARTEMENT") #971
            map_martinique            <<- readShapeSpatial("./ressources/ign_fonds_cartes/GEOFLA_2-0_SHP_UTM20W84MART_D972-ED141_DEPARTEMENT") #972
            map_guyane                <<- readShapeSpatial("./ressources/ign_fonds_cartes/GEOFLA_2-0_SHP_UTM22RGFG95_D973-ED141_DEPARTEMENT") # 973
            map_reunion               <<- readShapeSpatial("./ressources/ign_fonds_cartes/GEOFLA_2-0_SHP_RGR92UTM40S_D974-ED141_DEPARTEMENT") #974
            map_mayotte               <<- readShapeSpatial("./ressources/ign_fonds_cartes/GEOFLA_2-0_SHP_RGM04UTM38S_D976-ED141_DEPARTEMENT") #98
            # il manque 975=Saint Pierre et Miquelon
            save(map_france_metropolitaine, map_guadeloupe, map_martinique, map_guyane, map_reunion, map_mayotte, file=fichier2) ;
        } else if( !exists("map_france_metropolitaine")) {
            load(file=fichier2, envir=globalenv() ) ;
        }
        # ***********************************************************************
        
    } else if( ressource=="ign_geofla_communes"  ) {
        
        if(!require("maptools")) install.packages("maptools", repos="http://cran.us.r-project.org") ;
        library(maptools) ;
        gpclibPermit() ; # pour débloquer certaines fonctionnalités
        
        fichier2  <- file.path("ressources", "ign_geofla_communes.RData") ;
        if ( !file.exists(fichier2) ) {
            map_com_france_metropolitaine <<- readShapeSpatial("./ressources/ign_fonds_carte/GEOFLA_2-1_COMMUNE_SHP_LAMB93_FXX_2015-12-01") # métropole
            map_com_guadeloupe            <<- readShapeSpatial("./ressources/ign_fonds_carte/GEOFLA_2-1_COMMUNE_SHP_UTM20W84GUAD_D971_2015-12-01") #971
            map_com_martinique            <<- readShapeSpatial("./ressources/ign_fonds_carte/GEOFLA_2-1_COMMUNE_SHP_UTM20W84MART_D972_2015-12-01") #972
            map_com_guyane                <<- readShapeSpatial("./ressources/ign_fonds_carte/GEOFLA_2-1_COMMUNE_SHP_UTM22RGFG95_D973_2015-12-01") # 973
            map_com_reunion               <<- readShapeSpatial("./ressources/ign_fonds_carte/GEOFLA_2-1_COMMUNE_SHP_RGR92UTM40S_D974_2015-12-01") #974
            map_com_mayotte               <<- readShapeSpatial("./ressources/ign_fonds_carte/GEOFLA_2-1_COMMUNE_SHP_RGM04UTM38S_D976_2015-12-01") #98
            # il manque 975=Saint Pierre et Miquelon
            save(map_com_france_metropolitaine, map_com_guadeloupe, map_com_martinique, map_com_guyane, map_com_reunion, map_com_mayotte, file=fichier2) ;
            
        } else if( !exists("map_dep_france_metropolitaine")) {
            load(file=fichier2, envir=globalenv() ) ;
        }
        # ***********************************************************************
        
    } else if( ressource=="atih_terminologie_pmsi_mco"  ) {
        
        # il faudra encore faire import_CIM10, et créer de même les termino CCAM UCD et LPP.
        
        # fichier à encoder en ANSI
        fichier1 <- "ressources/atih_terminologie_pmsi_mco.csv" ;
        fichier2 <- "ressources/atih_terminologie_pmsi_mco.RData" ;
        if ( !file.exists(fichier2) ) {
            nomenclature_pmsi_mco <- fread(input=fichier1, colClasses = "character", sep=";") ;
            setkey(x=nomenclature_pmsi_mco, variable, valeur) ;
            atih_terminologie_pmsi_mco <<- nomenclature_pmsi_mco ;
            save( atih_terminologie_pmsi_mco, file = fichier2 ) ;
        } else if( !exists("atih_terminologie_pmsi_mco")) {
            load(file=fichier2, envir=globalenv() ) ;
        }
        # ***********************************************************************
        
    } else if( ressource=="atih_identifiants_etablissements") {
        
        fichier2 <- "ressources/atih_identifiants_etablissements.RData" ;
        if ( !file.exists(fichier2) ) {
            
            if(!require("stringdist")) install.packages("stringdist", repos="http://cran.us.r-project.org") ;
            library(stringdist) ;
            chemin_dossier <- "ressources/atih_liste_etablissements/";
            pattern_fichier <- "liste_etab.+\\.csv$";
            maxDist <- 3 ; # pour les recherches approximatives de noms de colonnes dans les fichiers
            
            file_list <- list.files( path = chemin_dossier, pattern = pattern_fichier, full.names = TRUE) ;
            for (un_fichier in file_list) {
                #   ***************************** ci-dessous : traitement d'un fichier donné
                # Récuperer l'année et le secteur Ã  partir du nom de fichier
                nom_fichier <- basename(un_fichier) ;
                regex_annee <-  "(?<=_)[[:digit:]]{4}(?=.csv)" ;
                m_annee <- regexpr(regex_annee, nom_fichier, perl = TRUE) ;
                annee <- regmatches(nom_fichier, m_annee) ;
                regex_secteur <-"(?<=etab_)[[:upper:]]+(?=_)" ;
                m_secteur <- regexpr(regex_secteur, nom_fichier, perl = TRUE) ;
                secteur <-  regmatches(nom_fichier, m_secteur) ;
                # Importer le fichier dans un data.frame (fichiers codés en latin-1)
                raw_etablissement <- read.csv2(
                    file=un_fichier, fileEncoding="latin1", stringsAsFactors=FALSE, check.names=FALSE,
                    strip.white=TRUE, colClasses="character", header=TRUE, quote="", row.names=NULL
                ) ;
                # On utilise la fonction amatch() du package stringdist car les noms de colonnes peuvent légèrement varier !!
                # La colonne "status" s'appelle "type valor" dans les fichiers MCO
                if (secteur == "MCO") {
                    nom_col_status <- "type valor"
                } else {
                    nom_col_status <- "status"
                }
                noms_colonnes <- colnames(raw_etablissement)
                indice_colonne <- function ( nom_colonne, liste=noms_colonnes, distance=maxDist) {
                    return(amatch(x = nom_colonne, table = liste, maxDist = distance))
                }
                # Créer un vecteur qui servira de dictionnaire : position des colonnes
                col_select <- c(
                    "finess" = indice_colonne("finess"),
                    "raison_sociale" = indice_colonne("raison sociale"),
                    "status" = indice_colonne(nom_col_status),
                    "region" = indice_colonne("region")
                ) ;
                # Sélectionner les bonnes colonnes
                clean_etablissement <- raw_etablissement[, col_select] ;
                # Changer les noms de colonnes par des noms standardisés
                colnames(clean_etablissement) <- names(col_select)
                # Ajouter l'année et le secteur, calculés plus haut
                clean_etablissement$annee <- as.integer(as.character(annee))
                clean_etablissement$secteur <- secteur
                # Calculer tout de suite le département
                clean_etablissement$dept <- trouve_departement(clean_etablissement$finess, type="finess") ;
                # Remettre le tout dans l'ordre
                clean_etablissement <- clean_etablissement[, c( "finess", "raison_sociale", "secteur", "status", "dept", "region", "annee")] ;
                donnees_un_etablissement <- clean_etablissement ;
                #   *****************************
                # à réincorporer
                if (exists("atih_identifiants_etablissements")) {
                    atih_identifiants_etablissements <<- rbind(atih_identifiants_etablissements, donnees_un_etablissement) ;
                } else {
                    atih_identifiants_etablissements <<- donnees_un_etablissement
                }
                atih_identifiants_etablissements <<- as.data.table(atih_identifiants_etablissements) ;
                # setkey(atih_identifiants_etablissements, finess, annee, secteur) ; # Pour simplifier : prendre la dernière valeur
                setkey(atih_identifiants_etablissements, finess) ;
            }
            save(atih_identifiants_etablissements, file=fichier2) ;
            
        } else if( !exists("atih_identifiants_etablissements")) {
            
            load(file=fichier2, envir=globalenv() ) ;
            
        }
        # ***********************************************************************
        
    } else if( ressource=="atih_termino_all") {
        
        fichier2  <- file.path("ressources", "atih_termino_all.RData") ;
        
        if ( !file.exists(fichier2) ) {
            rtrim <- function (x) { sub("\\s+$", "", x); }
            # CIM10 ATIH
            nomenclature_cim10 <- fread(input="./ressources/atih_terminologies/LIBCIM10.TXT", colClasses = "character", sep="|", header=FALSE) ;
            nomenclature_cim10 <- nomenclature_cim10[,list(V1,V4)]
            nomenclature_cim10[,V1 := rtrim(V1)] ;
            setnames(nomenclature_cim10,"V1","code") ;
            setnames(nomenclature_cim10,"V4","libelle") ;
            setkey(x=nomenclature_cim10, code) ;
            # CCAM ATIH
            nomenclature_ccam <- fread(input="./ressources/atih_terminologies/LIBCCAM.TXT", colClasses = "character", sep="|", header=FALSE) ;
            nomenclature_ccam <- nomenclature_ccam[,list(V1,V3)]
            nomenclature_ccam[,V1 := substring(V1,1,7)] ;
            setnames(nomenclature_ccam,"V1","code") ;
            setnames(nomenclature_ccam,"V3","libelle") ;
            setkey(x=nomenclature_ccam, code) ;
            # GHM MCO ATIH
            nomenclature_ghm <- fread(input="./ressources/atih_terminologies/LIBGHMFG.TXT", colClasses = "character", sep="|", header=FALSE) ;
            setnames(nomenclature_ghm,"V1","code") ;
            setnames(nomenclature_ghm,"V2","libelle") ;
            setkey(x=nomenclature_ghm, code) ;
            # Racines de GHM MCO ATIH
            nomenclature_rghm <- fread(input="./ressources/atih_terminologies/LIBRGHM.TXT", colClasses = "character", sep="\t", header=TRUE) ;
            setkey(x=nomenclature_rghm, code) ;
            # CMD MCO ATIH
            nomenclature_cmd <- fread(input="./ressources/atih_terminologies/LIBCMDFG.TXT", colClasses = "character", sep="|", header=FALSE) ;
            setnames(nomenclature_cmd,"V1","code") ;
            setnames(nomenclature_cmd,"V2","libelle") ;
            setkey(x=nomenclature_cmd, code) ;
            # LPP MCO ATIH
            nomenclature_lpp <- fread(input="./ressources/atih_terminologies/LIBLPP.TXT", colClasses = "character", sep="\t", header=TRUE) ;
            nomenclature_lpp <- nomenclature_lpp[code!="",]
            setkey(x=nomenclature_lpp, code) ;
            # UCD MCO ATIH
            nomenclature_ucd <- fread(input="./ressources/atih_terminologies/LIBUCD.TXT", colClasses = "character", sep="\t", header=TRUE) ;
            nomenclature_ucd <- nomenclature_ucd[code!="",]
            setkey(x=nomenclature_ucd, code) ;
            # UCD MCO ATIH DCI
            nomenclature_ucd_dci <- fread(input="./ressources/atih_terminologies/LIBUCDDCI.txt", colClasses = "character", sep="\t", header=TRUE) ;
            nomenclature_ucd_dci <- nomenclature_ucd_dci[code!="",]
            setkey(x=nomenclature_ucd_dci, code) ;
            # non faits mais fichiers présents en SSR : CM, GME, RGME, CM
            
            atih_terminologie_cim10 <<- nomenclature_cim10 ;
            atih_terminologie_ccam <<- nomenclature_ccam ;
            atih_terminologie_ghm <<- nomenclature_ghm ;
            atih_terminologie_cmd <<- nomenclature_cmd ;
            atih_terminologie_rghm <<- nomenclature_rghm ;
            atih_terminologie_lpp <<- nomenclature_lpp ;
            atih_terminologie_ucd <<- nomenclature_ucd ;
            atih_terminologie_ucd_dci <<- nomenclature_ucd_dci ;
            
            save( atih_terminologie_cim10, atih_terminologie_ccam, atih_terminologie_ghm, atih_terminologie_cmd, atih_terminologie_rghm, atih_terminologie_lpp, atih_terminologie_ucd, atih_terminologie_ucd_dci, file = fichier2 ) ;
        } else if( !exists("atih_terminologie_cim10")) {
            load(file=fichier2, envir=globalenv() ) ;
        }
        # ***********************************************************************    
    }
}

## Fonction calcul de taux de mortalité avec IC95% / année
## Obtention d'un dataframe

taux_morta_ic_years <- function(data = data) {
    stopifnot(is.data.frame(data))
    df_taux <- data.frame(annee = NULL, taux_morta = NULL, low_bound = NULL, upp_bound = NULL)
    for (une_annee in unique(data$annee)) {
        annee <- une_annee
        dc_bin <- ifelse(data$sortie_mode == 9 & data$annee == une_annee, 1, 0)
        x <- sum(dc_bin, na.rm=TRUE)
        n <- length(dc_bin[data$annee == une_annee])
        temp <- binom.test(x, n, alternative ="two.sided", conf.level = 0.95) ;  
        df_temp <- data.frame(annee = une_annee, 
                              taux_morta = temp$estimate[[1]],
                              low_bound = temp$conf.int[[1]],
                              upp_bound = temp$conf.int[[2]]
        )
        df_taux <- rbind(df_taux, df_temp)
    }
    return(df_taux)
}


# ***********************************************************************
# ***********************************************************************
# ***********************************************************************   
# *********************************************************************** 

# Fonction permettant d'effectuer les regroupements de codes à partir d'un fichier csv structuré en 2 colonnes (codes, regroupement)  

regroupement_modalite <- function(variable_a_regrouper, chemin_fichier, v_a_r_multivalue = TRUE) {
  
  # L'argument variable_a_regrouper est un vecteur qui contient les codes à regrouper
  # 
  # L'argument chemin_fichier correspond au chemin du fichier .csv (separateur ;) 
  # Le tableau doit au minimum contenir deux variables : "code" et "regroupement". La première contient la liste des codes et la seconde le nom du regroupement correspondant
  # 
  # L'argument v_a_r_multivalue indique si la variable à regrouper est multivaluée ou non
  # Si oui, la sortie sera un dataframe de variables binaires
  # Si non, la sortie sera un vecteur d'une variable qualitative à plusieurs modalitées
  # 
  regroupements <- read.csv2(chemin_fichier, stringsAsFactors = FALSE)
  
  if (!v_a_r_multivalue) {
    data_return <- rep(NA, length(variable_a_regrouper))
  } else {
    list_var <- unique(regroupements$regroupement)
    data_return <- as.data.frame(matrix(data = NA, 
                                        nrow = length(variable_a_regrouper), 
                                        ncol = length(list_var), 
                                        dimnames = list(NULL, list_var)))
  }
  
  stopifnot(c("code","regroupement") %in% names(regroupements))
  
  # Recodage des modalités NA en autres
  regroupements$regroupement[is.na(regroupements$regroupement)] <- "Autres"
  
  for (i in unique(regroupements$regroupement)) {
    list_a_matcher <- regroupements$code[regroupements$regroupement == i]
    if (v_a_r_multivalue) {
     list_a_matcher <- paste0("(", paste(list_a_matcher, collapse = ")|("), ")", collapse ="")
     data_return[,i] <- grepl(list_a_matcher, variable_a_regrouper)
    } else {
      data_return[variable_a_regrouper %in% list_a_matcher] <- i
    }
  }
  data_return[data_return == FALSE] <- 0
  data_return[data_return == TRUE] <- 1
  return(data_return)
}

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################


# Fonction de fusion sur les dataframes comprennant des colonnes avec le même nom et un contenu différent
# La fonction récupère le contenu des colonnes et renvoie un dataframe sans doublon de names()
## En input : df_split est le dataframe obtenu après l'utilisation de fonction split_quali_multi
## Les valeurs en double (2 x le même nom avec un code différent ne compte que pour 1)

aggreg_df <- function(df_split) {
# On renomme les colonnes en NA en "Autres"
names(df_split)[is.na(names(df_split)) == TRUE] <- "Autres"
# récupération des colonnes en double ou +
doublons <- (names(df_split)[which(table(names(df_split))>1)])
for( nom_col in unique(doublons)) {
vec_bol <- names(df_split) == nom_col
new_col <- apply(df_split[vec_bol], 1,sum)
df_split[nom_col] <- new_col
vec_new_df <- ifelse(duplicated(names(df_split)),0,1)
df_split <- df_split[,vec_new_df == 1]
df_split[df_split > 1] <- 1
return(df_split)
  }
}

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################


# Fonction renvoyant dataframe avec les noms de colonnes remplacés par new_names
change_col_names <- function(dataframe, new_names) {
  if(ncol(dataframe) == length(new_names)) {
    df_temp <- as.data.frame(0)
    for (i in 1:length(new_names)) {
      df_temp[new_names[i]] <- dataframe[1, i]
    }
    for (j in 1:(nrow(dataframe)-1)) {
      for (i in 1:length(new_names)) {
        df_temp[j+1, new_names[i]] <- dataframe[j+1, i]
      }
    }
    return(df_temp[,2:ncol(df_temp)]) ;
  }
  else {
    print("Err : Le vecteur ne contient pas le même nombre de noms qu'il y a de colonnes dans dataframe. ")
    return(dataframe) ;
  }
}


######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

######### Ajout NDE ##############


# Fonction pour remplacer des NA non reperés dans un fichier
remplacement_NA <- function(donnees, vecteur_NA){
    for (i in vecteur_NA){
        donnees <-  data.frame(sapply(donnees, function(x) {gsub(pattern = i, replacement = NA, x)}), stringsAsFactors = FALSE)
    }
  
  suppressWarnings(
    variable_numerique <- sapply(donnees, function(x){all(is.na(x) == is.na(as.numeric(x)))})
  )
  donnees[,variable_numerique] <- lapply(donnees[,variable_numerique], as.numeric)
  
  return(donnees)
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
	
	
######## Fonctions Desc avec output HTML (nécessite kableExtra)

# BINAIRE
desc_binaire_html <- function(vector, name="Variable", table=TRUE, ...) {
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
    names(df_tmp)=c("Modalité", "Effectif","Proportion","IC95%")
    print(knitr::kable(df_tmp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
    
    
    cat("\n<br>Calcul des IC95% à l'aide d'une loi binomiale\n<br>")
  }
  
  pie(table(vector)/length(vector), main=name, col=c("white", "cornflowerblue")) ;
}

	
# QUALI

desc_quali_html <- function(vector, name="Variable", table=TRUE, sort="alpha", limit_chart=Inf, tronque_lib_chart=20, ...) {
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
    cat("------------------------------------------------------------------------------------<br>") ;
    
    temp <- as.data.frame(table(vector)) ;
    modalites <- temp[,1] ;
    if( sort=="alpha") {
        modalites <- modalites[order(modalites)] ;
    } else if( sort=="croissant") {
        modalites <- modalites[order(temp[,2])] ;
    } else if( sort=="decroissant") {
        modalites <- modalites[order(0-temp[,2])] ;
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
    print(knitr::kable(df_tmp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
    }
    
        cat("\n<br>Calcul des IC95% à l'aide d'une loi binomiale\n<br>")
    
        
    # graphique maintenant
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

#QUANTI_DISC
desc_quanti_disc_html <- function(vector, name="Variable", mean_ci=TRUE, table=TRUE, Sum = T, sort="alpha", xlim=NULL, ...) {
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
    cat("Effectif analysé :", length(na.omit(vector)),"<br>") ;
    cat("------------------------------------------------------------------------------------<br>") ;
    if (Sum) {print (kable(rbind(as.matrix(summary(vector)), Sd = sd(vector) ))) } ;
    
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
       # cat("\nLes calculs des IC95% sont réalisés à partir de la loi binomiale\n")
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
    print(knitr::kable(df_tmp) %>% kable_styling(bootstrap_options = "striped", full_width = F))
    
    }
    if( mean_ci ) {
        sd <- sd(vector) ;
        mean <- mean(vector) ;
        n <- length(vector) ;
        cat( "\n<br>Moyenne et intervalle de confiance à 95% :",  
             round(mean,2),"[",round(mean-1.96*sd/sqrt(n),2),";",round(mean+1.96*sd/sqrt(n),2),"].<br>",
             "Calcul des IC95% à partir du théorème central limite") ;
    }
    plot(table(vector)/length(vector), xlab=name, ylab="proportion", col="cornflowerblue", xlim=xlim) ;
    
}
	
