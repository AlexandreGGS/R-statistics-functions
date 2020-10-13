# modification fonction Table1::make.table() buggée 

library(Table1)

mt = function (dat, cat.varlist = NULL, cat.header = names(dat[, cat.varlist]), 
    cat.rownames = lapply(lapply(dat[, cat.varlist], factor), 
        levels), cat.ptype = "None", cont.varlist = NULL, 
    cont.header = names(dat[, cont.varlist]), cont.ptype = "None", 
    strat = NULL, cat.rmstat = "None", cont.rmstat = "None", 
    dec = 2, pname = TRUE, colnames = NULL, output = "plain", 
    vspace = TRUE, varorder = "data", stripe = TRUE, stripe.col = "#F7F7F7", 
    header.style = "bold", factor.style = "bold", 
    stat.style = "plain", nowrap = TRUE, caption = NULL, footer = NULL, 
    tspanner = NULL, n.tspanner = NULL, cgroup = NULL, n.cgroup = NULL, col.columns = "none") 
{
    if (missing(dat) | (!is.data.frame(dat))) {
        warning("A data frame object must be provided in dat = .")
    }
    if (!(length(cat.varlist) == length(cat.header))) {
        warning("cat.varlist and cat.header must be the same length.")
    }
    if (!(length(cat.varlist) == length(cat.rownames))) {
        warning("cat.varlist and cat.rownames must be the same length.")
    }
    if (!(length(cont.varlist) == length(cont.header))) {
        warning("cont.varlist and cont.header must be the same length.")
    }
    if (is.null(strat)) {
        cat.strat = rep(list(strat), length(cat.varlist))
        cont.strat = rep(list(strat), length(cont.varlist))
        strat.miss = 0
        strat.rem = 0
    }
    else if (!is.null(strat)) {
        cat.strat = rep(list(interaction(sapply(strat, FUN = get, 
            pos = dat, simplify = FALSE, USE.NAMES = TRUE))), 
            length(cat.varlist))
        cont.strat = rep(list(interaction(sapply(strat, FUN = get, 
            pos = dat, simplify = FALSE, USE.NAMES = TRUE))), 
            length(cont.varlist))
        strat.miss = lapply(sapply(strat, FUN = get, pos = dat, 
            simplify = FALSE, USE.NAMES = TRUE), function(x) sum(is.na(x)))
        strat.rem = sum(is.na(interaction(sapply(strat, FUN = get, 
            pos = dat, simplify = FALSE, USE.NAMES = TRUE))))
    }
    if (!is.null(strat) & strat.rem > 0) {
        print(paste("Total observations removed from table:", 
            strat.rem, sep = " "))
        print("Summary of total missing stratification variable(s):")
        print(strat.miss)
        footer.miss <- paste(strat.rem, "observations removed due to missing values in \n                         stratifying variable(s)", 
            sep = " ")
    }
    else {
        footer.miss <- paste("")
    }
    if (is.null(cont.varlist)) {
        cats <- mapply(cat.var, var = sapply(cat.varlist, FUN = get, 
            pos = dat, simplify = FALSE, USE.NAMES = TRUE), strat = cat.strat, 
            cat.rmstat = cat.rmstat, dec = dec, rownames = cat.rownames, 
            header = cat.header, ptype = cat.ptype, pname = pname, 
            vspace = TRUE, SIMPLIFY = FALSE)
        if (varorder == "data") {
            tab <- do.call(rbind, (c(cats))[order(match(names(c(cats)), 
                names(dat)))])
        }
        else if (varorder == "abc") {
            tab <- do.call(rbind, (c(cats))[order(match(names(c(cats)), 
                sort(names(dat))))])
        }
        if (all(cat.ptype == "None")) {
            tab <- tab[, -dim(tab)[2]]
        }
        tab[grepl("NaN", tab)] <- "-"
        tab[grepl("NA", tab)] <- "-"
        tab[grepl("-Inf", tab)] <- "-"
    }
    else if (is.null(cat.varlist)) {
        conts <- mapply(cont.var, var = sapply(cont.varlist, 
            FUN = get, pos = dat, simplify = F, USE.NAMES = T), 
            strat = cont.strat, cont.rmstat = cont.rmstat, dec = dec, 
            header = cont.header, ptype = cont.ptype, pname = pname, 
            vspace = TRUE, SIMPLIFY = FALSE)
        if (varorder == "data") {
            tab <- do.call(rbind, (c(conts))[order(match(names(c(conts)), 
                names(dat)))])
        }
        else if (varorder == "abc") {
            tab <- do.call(rbind, (c(conts))[order(match(names(c(conts)), 
                sort(names(dat))))])
        }
        if (all(cont.ptype == "None")) {
            tab <- tab[, -dim(tab)[2]]
        }
        tab[grepl("NaN", tab)] <- "-"
        tab[grepl("NA", tab)] <- "-"
        tab[grepl("-Inf", tab)] <- "-"
    }
    else {
        cats <- mapply(cat.var, var = sapply(cat.varlist, FUN = get, 
            pos = dat, simplify = FALSE, USE.NAMES = TRUE), strat = cat.strat, 
            cat.rmstat = cat.rmstat, dec = dec, rownames = cat.rownames, 
            header = cat.header, ptype = cat.ptype, pname = pname, 
            vspace = TRUE, SIMPLIFY = FALSE)
        conts <- mapply(cont.var, var = sapply(cont.varlist, 
            FUN = get, pos = dat, simplify = FALSE, USE.NAMES = TRUE), 
            strat = cont.strat, cont.rmstat = cont.rmstat, dec = dec, 
            header = cont.header, ptype = cont.ptype, pname = pname, 
            vspace = TRUE, SIMPLIFY = FALSE)
        if (varorder == "data") {
            tab <- do.call(rbind, (c(cats, conts))[order(match(names(c(cats, 
                conts)), names(dat)))])
        }
        else if (varorder == "abc") {
            tab <- do.call(rbind, (c(cats, conts))[order(match(names(c(cats, 
                conts)), sort(names(dat))))])
        }
        if (all(cat.ptype == "None" & all(cont.ptype == 
            "None"))) {
            tab <- tab[, -dim(tab)[2]]
        }
        tab[grepl("NaN", tab)] <- "-"
        tab[grepl("NA", tab)] <- "-"
        tab[grepl("-Inf", tab)] <- "-"
    }
    if (!is.null(colnames)) {
        colnames <- colnames
    }
    else if (is.null(colnames)) {
        colnames <- colnames(tab)
    }
    if (output == "plain") {
        out.plain(tab, colnames = colnames)
    }
    else if (output == "html") {
        out.html(tab, colnames = colnames, stripe = stripe, stripe.col = stripe.col, 
            header.style = header.style, factor.style = factor.style, 
            stat.style = stat.style, caption, footer, tspanner, 
            n.tspanner, cgroup, n.cgroup, col.columns = col.columns, 
            nowrap = nowrap)
    }
    else if (output == "latex") {
        out.latex(tab, colnames = colnames, header.style = header.style, 
            factor.style = factor.style, stat.style = stat.style)
    }
}
