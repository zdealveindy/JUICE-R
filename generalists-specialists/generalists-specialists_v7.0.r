# JUICE-R JUICE.table JUICE.species.data R.exe no.transf http://www.davidzeleny.net/juice-r/doku.php?id=scripts:generalists-specialists

lib.admin <- function (packages, install.string = NULL)
{
  require (tcltk)
  if (is.null (install.string)) 
  {
    install.string <- list (NA)
    for (i in seq (1, length (packages))) install.string [[i]] <- NA
  }
  for (i in seq (1, length (packages))) if (is.na (install.string [[i]][1])) install.string [[i]] <- as.character (as.expression (call ("install.packages", pkgs = packages [i], repos = "http://cran.rstudio.com")))
  
  admin.list <- list ()
  admin.label <- list ()
  admin.button <- list ()
  cancel.admin<- tclVar (0)
  check.install <- function (packages) lapply (packages, FUN = function (x) require (x, character.only = T))
  reconfigure <- function (packages, admin.list, admin.label, admin.button) 
  {
    installed.packages <- lapply (packages, FUN = function (x) require (x, character.only = T))
    for (i in seq (1, length (packages))) tkconfigure (admin.button[[i]], state = if (installed.packages[[i]]) 'disabled' else 'normal')
    for (i in seq (1, length (packages))) tkconfigure (admin.label[[i]], text = paste (packages[i], ' is ', ifelse (installed.packages[[i]], '', 'not '), 'installed', sep = ''))
  }
  
  installed.packages <- check.install (packages)
  # create the tcltk wizard  
  admin <- tktoplevel ()
  tkwm.title (admin, 'Installed libraries')
  tkpack (tklabel (admin, text = 'Administration of R packages'))
  for (i in seq (1, length (packages))) admin.list[[i]] <- tkframe (admin)
  for (i in seq (1, length (packages))) admin.label[[i]] <- tklabel (admin.list[[i]], width = 20, text = paste (packages[i], ' is ', ifelse (installed.packages[[i]], '', 'not '), 'installed', sep = ''))  
  for (i in seq (1, length (packages))) admin.button[[i]] <- tkbutton (admin.list[[i]], width = 15, text = paste ('Install', packages[i], sep = ' '), command = function (x) {eval (parse (text = install.string[[i]])); reconfigure (packages[i], admin.list[i], admin.label[i], admin.button[i])})
  for (i in seq (1, length (packages))) tkconfigure (admin.button[[i]], state = if (installed.packages[[i]]) 'disabled' else 'normal')
  for (i in seq (1, length (packages))) tkpack (admin.label[[i]], admin.button[[i]], side = 'right', anchor = 'w', side = 'left' )
  for (i in seq (1, length (packages))) tkpack (admin.list [[i]])
  tkpack (tkbutton (admin, width = 20, text = 'Continue', command = function() {tclvalue(cancel.admin)<-1; tkdestroy (admin)} ))
  tkbind (admin, '<Destroy>', function() {tclvalue(cancel.admin)<-1; tkdestroy (admin)})
  tkwait.variable (cancel.admin)
}


packages <- c('devtools', 'parallel', 'betapart', 'vegan', 'genspe')
if (!all (unlist (lapply (packages, FUN = function (x) require (x, character.only = T))))) lib.admin (packages,  install.string = list (NA, NA, NA, NA, as.character (as.expression (c (call ('require', 'devtools'), call ('install_github', 'zdealveindy/genspe'))))))

require (genspe)
require (tcltk)
beals.file <- NA

input.matrix <- as.matrix (JUICE.table)
species.data <- JUICE.species.data[,c('Species_Name', 'Layer')]
rownames (species.data) <- JUICE.species.data[, 'SpecName_Layer']
names (species.data) <- c('full.sci.name', 'layer')

calculate.theta.tcltk (input.matrix = input.matrix, species.data = species.data, juicer = T)