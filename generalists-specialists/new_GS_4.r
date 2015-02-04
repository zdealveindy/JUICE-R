#JUICE-R JUICE.table R.exe no.transf

beals.2 <- function (x, type = 0, include = TRUE) 
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

GS.function <- function (input.matrix, reps, psample, version = "multiplicative", beals.file) # the function calculating theta value of species specialization
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

input.matrix[input.matrix > 0] <- 1
write.csv2 (input.matrix, file = 'original-data.csv')
if (version == "beals")
{
  if (is.na (beals.file))
  {
  require (vegan)
  beals.matrix <- beals.2 (input.matrix, include = F)
  } else {
  win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix)+1, initial = 0, width = 300) 
  setWinProgressBar (win.pb, 1, label = "Reading Beals smoothed table calculated by JUICE")
  beals.matrix <- read.csv2 (file = beals.file, head = F)
  setWinProgressBar (win.pb, 1, label = "Transforming JUICE Beals smoothing file into matrix format")
  beals.matrix <- beals.matrix [,c(1,2,4)]
  beals.matrix[,3] <- as.numeric (as.character (beals.matrix[,3]))
  beals.matrix <- matrix (beals.matrix[,3], ncol = ncol (input.matrix), byrow = T)
  close (win.pb)
  }
  if (! all (dim (beals.matrix) == dim (input.matrix))) stop ('Beals matrix has different size than species matrix!')
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
  write.csv2 (input.matrix, file = 'beals-data.csv')
}


Nplots <- dim (input.matrix)[1]

plots.per.spp <- colSums (input.matrix)
select.spp <- plots.per.spp[plots.per.spp >= psample]

Nspp <- length (select.spp)

win.pb <- winProgressBar (title = "Calculation progress", label = paste ("Species no. ", 1), min = 1, max = Nspp, initial = 0, width = 300) 

sci.name <- rep(0,Nspp)
meanco <- rep(0,Nspp)
meanco.sd <- rep(0,Nspp)
local.avgS <- rep(0,Nspp)
occur.freq <- rep(0,Nspp)
GS <- rep(0,Nspp)
GS.sd <- rep(0,Nspp)

for (sp in 1:Nspp)
{
setWinProgressBar (win.pb, sp, label = paste ("Species no. ", sp))

temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == labels (select.spp[sp])]>0,]
temp.matrix <- temp.matrix[,colSums (temp.matrix > 0) > 0]

#spp.matrix <- as.data.frame (matrix (nrow = reps, ncol = 2, dimnames = list (NULL, c ("mean.alpha", "total.rich"))))

rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = dim (temp.matrix)[1], byrow = F)
sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))

mc.mat <- array(0,dim=c(psample,dim (temp.matrix)[2],reps))	
for(i in 1:reps) {
	mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
}
total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
mean.alpha <- colMeans (apply (mc.mat, c(1,3), sum))


if (version == "multiplicative" | version == "beals") Wbeta.vec <- total.rich/mean.alpha else Wbeta.vec <- total.rich-mean.alpha # here is the option to choose between multiplicative Zeleny"s algorithm and additive Fridley"s algorithm

GS[sp] <- mean(Wbeta.vec)			#mean Whittaker beta value for all reps (= THETA: G-S metric)
GS.sd[sp] <- sd(Wbeta.vec)			#s.d. of above
meanco[sp] <- mean(total.rich)			#mean # cooccurrences in "psample" plots
meanco.sd[sp] <- sd(total.rich)		#s.d. of above

sci.name[sp] <- labels (select.spp[sp])	#scientific name
local.avgS[sp] <- mean(mean.alpha)				#approximate mean local richness
occur.freq[sp] <- as.vector (select.spp[sp])							#total number of plots
}

close (win.pb)

#Output matrix
meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
theta.out <- data.frame(sci.name,local.avgS,occur.freq,meanco,meanco.sd,meanco.u,meanco.l,GS,GS.sd)
tkdestroy (base)
write.csv2 (file = "theta_out.csv", theta.out)

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
cancel <- tclVar (0)
end.end <- F
beals.file <- NA

GSmethods <- c ("Multiplicative beta on species pool (Botta-Dukat 2011)", "Multiplicative beta diversity (Zeleny 2009)", "Additive beta diversity (Fridley et al. 2007)")
base <- tktoplevel()
tkwm.title(base, "Calculation of generalists-specialists metric")

  spec.frm <- tkframe(base,borderwidth=2)
  frame.a <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.b <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.c <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.d <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.e <- tkframe(spec.frm, relief = 'groove', borderwidth = 2)
  
  GSmethod <- tclVar ("multiplicative")
  label.radio <- tklabel (frame.a, text = "Which beta diversity algorithm to use?")
  radio1 <- tkradiobutton (frame.a, text = GSmethods[1], value = "beals", variable = GSmethod)
  radio2 <- tkradiobutton (frame.a, text = GSmethods[2], value = "multiplicative", variable = GSmethod)
  radio3 <- tkradiobutton (frame.a, text = GSmethods[3], value = "additive", variable = GSmethod)
  
  tk.psample <- tclVar (5)
  tk.reps <- tclVar (10)
  label.entry1 <- tklabel (frame.b, text = "Minimal frequency of species")
  entry1 <- tkentry (frame.b, width = 5, textvariable = tk.psample)
  
  label.entry2 <- tklabel (frame.c, text = "Number of random subsamples")
  entry2 <- tkentry (frame.c, width = 5, textvariable = tk.reps)
  
  button1 <- tkbutton (frame.d, text = "Calculate", command = function () end.end <- GS.function (as.matrix (JUICE.table), psample = as.numeric (tclvalue (tk.psample)), reps = as.numeric (tkget (entry2)), version = as.character (tclvalue (GSmethod)), beals.file = beals.file))

  tkpack (tklabel (frame.e, text = 'Choose file with beals smoothed table\n (if available and needed)'), tkbutton (frame.e, text = 'Choose', command = function () assign ('beals.file', choose.files (), env = .GlobalEnv)))
  
  tkpack (label.radio, radio1, radio2, radio3)
  tkpack (label.entry1, entry1)
  tkpack (label.entry2, entry2)
  tkpack (button1)
  
  tkpack (frame.a, side = "top", ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w")
  tkpack (frame.e, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w")
  tkpack (frame.b, frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10)
  tkpack (frame.d, side = "bottom", pady = 10, padx = 10)

  tkpack (spec.frm)
  tkbind (base, "<Destroy>", function() tclvalue(cancel)<-2)  
  
  windows ()
  bringToTop (dev.cur (), stay = F)
  tkraise (base)
  #dev.off ()
  tkwait.variable (cancel)


