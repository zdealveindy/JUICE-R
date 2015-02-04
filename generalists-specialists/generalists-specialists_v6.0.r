# JUICE-R JUICE.table R.exe no.transf http://www.davidzeleny.net/juice-r/doku.php?id=scripts:generalists-specialists

# dependent libraries: parallel and betapart
packages <- c('parallel', 'betapart') 


lib.admin <- function (packages, install.string = NULL)
{
  require (tcltk)
  if (is.null (install.string)) install.string <- rep (NA, length (packages))
  for (i in seq (1, length (packages))) if (is.na (install.string [i])) install.string [i] <- as.character (as.expression (call ('install.packages', pkgs = packages [i], repos = "http://cran.rstudio.com")))
  
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
  for (i in seq (1, length (packages))) admin.button[[i]] <- tkbutton (admin.list[[i]], width = 15, text = paste ('Install', packages[i], sep = ' '), command = function (x) {eval (parse (text = install.string[i])); reconfigure (packages, admin.list, admin.label, admin.button)})
  for (i in seq (1, length (packages))) tkconfigure (admin.button[[i]], state = if (installed.packages[[i]]) 'disabled' else 'normal')
  for (i in seq (1, length (packages))) tkpack (admin.label[[i]], admin.button[[i]], side = 'right', anchor = 'w', side = 'left' )
  for (i in seq (1, length (packages))) tkpack (admin.list [[i]])
  tkpack (tkbutton (admin, width = 20, text = 'Continue', command = function() {tclvalue(cancel.admin)<-1; tkdestroy (admin)} ))
  tkbind (admin, '<Destroy>', function() {tclvalue(cancel.admin)<-1; tkdestroy (admin)})
  tkwait.variable (cancel.admin)
}

beals.2 <- function (x, include = TRUE) # version of beals from vegan, for only p/a data and with progress bar
{
win.pb2 <- winProgressBar (title = 'Beals', label = 'start', min = 1, max = nrow (x)+2, initial = 0, width = 300)
x <- as.matrix(x)
x [x > 0] <- 1
refX <- x
incSp <- include
refX <- as.matrix(refX)
setWinProgressBar (win.pb2, 1, label = 'Crossprod')
M <- crossprod(refX, refX)
C <- diag(M)
setWinProgressBar (win.pb2, 1, label = 'First sweep')
M <- sweep(M, 2, replace(C, C == 0, 1), "/")
if (!incSp) for (i in 1:ncol(refX)) M[i, i] <- 0
S <- rowSums(x)
b <- x
for (i in 1:nrow(x)) {
  setWinProgressBar (win.pb2, i+1, label = i)
  b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
  }                       
SM <- rep(S, ncol(x))
if (!incSp) SM <- SM - x
b <- b/replace(SM, SM == 0, 1)
close (win.pb2)
b
}

calculate.theta <- function (input.matrix, sp, select.spp, remove.out, psample, reps, version, parallel, win.pb)
{
  if (parallel) write (paste (sp, '\n'), file = 'GS-progress.r', append = T) else setWinProgressBar (win.pb, sp, label = paste ("Species no. ", sp))
  temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,]
  temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  
  if (remove.out)
  {
    veg.dist <- as.matrix (dist (temp.matrix))
    diag (veg.dist) <- NA
    distances <- rowMeans (veg.dist, na.rm = T)
    outliers <- distances > (mean (distances) + 2*sd (distances))
    temp.matrix <- temp.matrix[!outliers,]
    temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  }
  
  if (!nrow (temp.matrix) < psample)  # not sure why there is this condition
  {
    rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = dim (temp.matrix)[1], byrow = F)
    sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))
    
    mc.mat <- array(0,dim=c(psample,dim (temp.matrix)[2],reps))  
    for(i in 1:reps) mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
    total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
    mean.alpha <- colMeans (apply (mc.mat, c(1,3), sum))
    
    if (version == "multiplicative" | version == "beals") Wbeta.vec <- total.rich/mean.alpha
    if (version == "additive") Wbeta.vec <- total.rich-mean.alpha 
    if (version == "pairwise.jaccard") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'jaccard')$beta.jac)))
    if (version == "multiple.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SOR))
    
    GS <- mean(Wbeta.vec)  		#mean Whittaker beta value for all reps (= THETA: G-S metric)
    GS.sd <- sd(Wbeta.vec)			#s.d. of above
    meanco <- mean(total.rich)			#mean # cooccurrences in "psample" plots
    meanco.sd <- sd(total.rich)		#s.d. of above
    
    sci.name <- labels (select.spp[sp])	#scientific name
    local.avgS <- mean(mean.alpha)				#approximate mean local richness
    occur.freq <- as.vector (select.spp[sp])							#total number of plots
    
    meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
    meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
    result <- c(sci.name,local.avgS,occur.freq,meanco,meanco.sd,meanco.u,meanco.l,GS,GS.sd)
    return (result)
  }
}

GS.function <- function (input.matrix, reps, psample, version = "multiplicative", beals.file, parallel = parallel, no.cores = no.cores, remove.out = remove.out) # the function calculating theta value of species specialization
{
if ( is.na (reps) || reps < 2) {
  tkmessageBox (type = "ok", message = "Number of random subsamples must be integer >= 2")
  end.end <- F
  } else
  if (is.na (psample) || psample < 2) 
    {
    tkmessageBox (type = "ok", message = "Minimal frequency of species must be integer >= 2")
    end.end <- F
    } else
    {

input.matrix <- ifelse (input.matrix > 0, 1, 0)
Nplots <- dim (input.matrix)[1]
plots.per.spp <- colSums (input.matrix)
select.spp <- plots.per.spp[plots.per.spp >= psample]
Nspp <- length (select.spp)

if (version == "beals")
{
  if (is.na (beals.file))
  {
  beals.matrix <- beals.2 (input.matrix, include = T)
  win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix), initial = 0, width = 300) 
  for (co in seq (1, ncol (input.matrix)))
    {
    setWinProgressBar (win.pb, co + 1, label = paste ("Prepare beals smoothed table: species ", co))
    beals.temp <- beals.matrix[,co][as.logical (input.matrix[,co])]
    stats.temp <- fivenum (beals.temp)
    iqr <- diff (stats.temp [c(2,4)])
    beals.thresh <- min (beals.temp[!(beals.temp < stats.temp[2] - 1.5 * iqr)])
    input.matrix[,co] <- as.numeric (beals.matrix[,co] >= beals.thresh)
    }
  close (win.pb)
  write.table (input.matrix, file = 'beals-data.txt', sep = '\t', col.names = TRUE)
  } else  
  {
    win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix)+1, initial = 0, width = 300) 
    setWinProgressBar (win.pb, 1, label = "Reading Beals smoothed table")
    beals.matrix <- as.matrix (read.delim (file = beals.file, row.names = 1, check.names = F))
    if (! all (dim (beals.matrix) == dim (input.matrix))) {tkmessageBox (type = "ok", message = paste ("Selected Beals matrix has different size than species matrix! \nYou need to calculate new beals smoothed species pool data using current species data. Close the JUICE-R application and run it again from JUICE, and calculate the Multiplicative beta on species pool analysis without selecting the Beals smoothing table.")); stop ('Beals matrix has different size than species matrix!')}
    input.matrix <- beals.matrix
    close (win.pb)
  }
}

if (!parallel) 
{
win.pb <- winProgressBar (title = "Calculation progress", label = paste ("Species no. ", 1), min = 1, max = Nspp, initial = 0, width = 300) 
temp.res <- lapply (1:Nspp, FUN = function (sp) calculate.theta (input.matrix = input.matrix, sp = sp, select.spp = select.spp, remove.out = remove.out, psample = psample, reps = reps, version = version, parallel = parallel, win.pb = win.pb))
close (win.pb)
theta.out <- as.data.frame (matrix (unlist (temp.res), ncol = 9, byrow = T, dimnames = list (NULL, c('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'meanco.sd', 'meanco.u', 'meanco.l', 'GS', 'GS.sd'))))
write.table (theta.out, file = 'theta_out.txt', sep = '\t', qmethod = 'double', col.names = NA)
}

if (parallel)
{
require (parallel)
workers <- makeCluster (no.cores)
if (file.exists ('GS-progress.r')) file.remove ('GS-progress.r')
clusterExport (workers, c('calculate.theta', 'tclvalue', 'input.matrix', 'select.spp', 'remove.out', 'psample', 'reps', 'version', 'parallel'), envir = environment ())
temp.res <- parLapplyLB (workers, 1:Nspp, fun = function (sp) calculate.theta (input.matrix = input.matrix, sp = sp, select.spp = select.spp , remove.out = remove.out, psample = psample, reps = reps, version = version, parallel = parallel, win.pb = NULL))
theta.out <- as.data.frame (matrix (unlist (temp.res), ncol = 9, byrow = T, dimnames = list (NULL, c('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'meanco.sd', 'meanco.u', 'meanco.l', 'GS', 'GS.sd'))))
write.table (theta.out, file = 'theta_out.txt', sep = '\t', qmethod = 'double', col.names = NA)
stopCluster (workers)
}

spd <- matrix (unlist (strsplit (as.character(theta.out[,1]), split = "_")), ncol = 2, byrow = T)
spaces <- 50-unlist (lapply ((strsplit(spd[,1], split = "" )), length))
rep.spaces <- function (times) paste (rep (" ", times), collapse = "")
write.table (file = "theta_import.species.data.txt", paste (spd[,1], unlist (apply (as.matrix (spaces), 1, rep.spaces)), theta.out[,8]), quote = F, row.names = F, col.names = F)

tkmessageBox (type = "ok", message = paste ("Result files were saved into", getwd (), "\n\nYou can use the file theta_import.species.data.txt to import the species theta values directly to JUICE (File > Import > Species Data)."))
end.end <- T
  }
cancel <- tclVar (1)
}



library (tcltk)
if (!all (unlist (lapply (packages, FUN = function (x) require (x, character.only = T))))) lib.admin (packages)
cancel <- tclVar (0)
end.end <- F
beals.file <- NA

GSmethods <- c ("Additive beta diversity (Fridley et al. 2007)", "Multiplicative beta diversity (Zeleny 2009)", "Multiplicative beta on species pool (Botta-Dukat 2012)", "Pairwise Jaccard dissimilarity (Manthey & Fridley 2009)", "Multiple Sorensen dissimilarity (Manthey & Fridley (2009)")
base <- tktoplevel()
tkwm.title(base, "Generalists-specialists")

  spec.frm <- tkframe (base, borderwidth=2)
  frame.title <- tkframe (spec.frm, relief = 'groove', borderwidth = 2, background = 'grey')
  frame.a <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.b <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.c <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.d <- tkframe (spec.frm, borderwidth = 2)
  frame.e <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.e.1 <- tkframe (frame.e)
  frame.e.2 <- tkframe (frame.e)
  frame.f <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.f.1 <- tkframe (frame.f)
  frame.g <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  
  
  GSmethod <- tclVar ("additive")
  label.radio <- tklabel (frame.a, text = "Which beta diversity algorithm to use?")
  radio1 <- tkradiobutton (frame.a, text = GSmethods[1], value = "additive", variable = GSmethod)
  radio2 <- tkradiobutton (frame.a, text = GSmethods[2], value = "multiplicative", variable = GSmethod)
  radio3 <- tkradiobutton (frame.a, text = GSmethods[3], value = "beals", variable = GSmethod)
  radio4 <- tkradiobutton (frame.a, text = GSmethods[4], value = "pairwise.jaccard", variable = GSmethod)
  radio5 <- tkradiobutton (frame.a, text = GSmethods[5], value = "multiple.sorensen", variable = GSmethod)
  tk.psample <- tclVar (5)
  tk.reps <- tclVar (10)
  parallel <- tclVar (0)
  no.cores <- tclVar (2)
  remove.out <- tclVar (0)
  label.entry1 <- tklabel (frame.b, text = "Minimal frequency of species ")
  entry1 <- tkentry (frame.b, width = 5, textvariable = tk.psample)
  
  label.entry2 <- tklabel (frame.c, text = "Number of random subsamples ")
  entry2 <- tkentry (frame.c, width = 5, textvariable = tk.reps)
  
  button1 <- tkbutton (frame.d, text = "Calculate", command = function () end.end <- GS.function (as.matrix (JUICE.table), psample = as.numeric (tclvalue (tk.psample)), reps = as.numeric (tkget (entry2)), version = as.character (tclvalue (GSmethod)), beals.file = beals.file, parallel = as.logical (as.numeric (tclvalue (parallel))), no.cores = as.numeric (tclvalue (no.cores)), remove.out = as.logical (as.numeric (tclvalue (remove.out)))))


  choose.label <- tklabel (frame.e.2, text = 'Select the file with beals smoothed data')
  choose.button <- tkbutton (frame.e.1, text = 'Select', command = function () assign ('beals.file', choose.files (), env = .GlobalEnv))
  tkpack (choose.button)
  tkpack (choose.label)
  tkpack (tklabel (frame.e, text = 'Beals smoothing (included in method of Botta-Dukat 2012)'), anchor = 'w')
  tkpack (frame.e.1, frame.e.2, side = 'left',ipady = 5, ipadx = 5, padx = 5, pady = 5)

  
  parallel.label <- tklabel (frame.f, text = 'Parallel calculation (enable only if you have more than one core)')
  parallel.no.cores.label <- tklabel (frame.f.1, text = 'number of cores: ')
  parallel.no.cores.entry <- tkentry (frame.f.1, width = 2, textvariable = no.cores)
  parallel.checkbutton <- tkcheckbutton (frame.f.1, text = 'use parallel computing,', variable = parallel)
  
  tkpack (tklabel (frame.g, text = 'Outlier analysis (McCune & Mefford 1999, suggested by Botta-Dukat 2012)'), tkcheckbutton (frame.g, text = 'remove outlier samples (with very different species composition)', variable = remove.out), anchor = 'w')
  
  tkpack (label.radio, radio1, radio2, radio3, radio4, radio5, anchor = 'w')
  tkpack (label.entry1, entry1, anchor = 'w', side = 'left')
  tkpack (label.entry2, entry2, anchor = 'w', side = 'left')
  tkpack (button1)
  tkpack (parallel.checkbutton, parallel.no.cores.label, parallel.no.cores.entry, side = 'left')
  tkpack (parallel.label,  frame.f.1, anchor = 'w')
  
  tkpack (tklabel (frame.title, text = 'Calculation of generalists and specialists using co-occurrence species data \n Author: David Zeleny (zeleny.david@gmail.com) \n JUICE-R application (www.bit.ly/habitat-specialists)'), ipady = 10, ipadx = 10, padx = 10, pady = 10)
  tkpack (frame.title, side = 'top', expand = T, fill = 'both')
  tkpack (frame.a, side = "top", ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.e, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.f, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.g, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.b, frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10, expand = T, fill = 'both')
  tkpack (frame.d, side = "bottom", pady = 10, padx = 10, expand = T, fill = 'both')

  tkpack (spec.frm)
  tkbind (base, "<Destroy>", function() tclvalue(cancel)<-2)  
  
  tkraise (base)
  tkwait.variable (cancel)


