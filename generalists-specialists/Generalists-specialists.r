# Calculating the measure of species habitat specialization from cooccurrence vegetation data
# Author: David Zeleny, zeleny@sci.muni.cz
# Version 1.0
# Notes: offers calculation of Theta-value according to Fridley et al. 2007 and Zeleny 2009


# Definition of functions
convert.table.2c <- function (veg.data)   # converts sample x species vegetation table into column form
{
c2 <- as.data.frame.table (as.table (as.matrix (veg.data)))
as.matrix(c2[c2$Freq != 0, c(2,1)])
}

GS.function <- function (result.matrix, reps, psample, version = 'multiplicative') # the function calculating theta value of species specialization
{

if ( is.na (reps) || reps < 2) {
  tkmessageBox (type = 'ok', message = 'Number of random subsamples must be integer >= 2')
  end.end <- F
  } else
  if (is.na (psample) || psample < 2) 
    {
    tkmessageBox (type = 'ok', message = 'Minimal frequency of species must be integer >= 2')
    end.end <- F
    } else
    {

GOSmat <- cbind (result.matrix[,2], factor (result.matrix[,1]))
SppMat <- cbind (1:length(levels(factor(result.matrix[,1]))), levels (factor(result.matrix[,1])))

#Appendix S1: R code implementing generalist-specialist metric based on co-occurrences (Fridley et al.)

#Jason Fridley, January 2007; email: fridley@syr.edu
#Implemented in R 2.3.1

############################################################################
#INPUT DATA
#
#import co-occurrence dataset as a 2-column text file ("list" format): column 1 = plot, column 2 = species code (ie species must be NUMERIC)
#GOSmat <-
#
#import labels dataset as a 2-column text file: column 1 = numeric species code (links to above), column 2 = species label (name, etc)
#SppMat <-
#
#set parameters
reps <- reps	#number of random samples per species
psample <- psample	#min plot occurrences for a species to be analyzed
#
############################################################################

#OUTPUT is "Theta.out" table

#Do not edit below this line------------------------------------------------

plotID <- factor(GOSmat[,1])		#factorized plot vector
SppID <- GOSmat[,2]			#species per plot as numbers or codes
Nplots <- length(levels(plotID))	#number of plots
richness <- tapply(SppID,plotID,length)	#vector of number of species in each plot
max.rich <- max(richness)		#maximum local richness value
metacom <- table(plotID,SppID)		#plot x species matrix

#Select subset of species for analysis that occur in "plot.cut" # of plots or more
plots.per.spp <- tapply(plotID,SppID,length)				#vector of number plot occurrences for each species
Species <- as.numeric(labels(plots.per.spp[plots.per.spp>=psample])[[1]])	#vector of selected species
Nspp <- length(Species)							#number of selected species

#SPECIES LOOP

sci.name <- rep(0,Nspp)
meanco <- rep(0,Nspp)
meanco.sd <- rep(0,Nspp)
local.avgS <- rep(0,Nspp)
tot.cooccur <- rep(0,Nspp)
occur.freq <- rep(0,Nspp)
GS <- rep(0,Nspp)
GS.sd <- rep(0,Nspp)

win.pb <- winProgressBar (title = 'Calculation progress', label = paste ('Species no. ', 1), min = 1, max = Nspp, initial = 0, width = 300) 

for(sp in 1:Nspp) {

setWinProgressBar (win.pb, sp, label = paste ('Species no. ', sp))

#Plot selection
lab <- as.numeric(labels(metacom)[2][[1]])
xlab <- c(1:dim(metacom)[2])
metacol <- xlab[lab==Species[sp]]
sp.plots <- as.logical(metacom[,metacol])
sp.metacom <- metacom[sp.plots,]
Np <- dim(sp.metacom)[1]
wide <- length(xlab)

#Monte Carlo procedure
rpmat <- matrix(c(1:Np),reps,Np,byrow=T)				#"reps" rows of plot sequences
rpmat <- t(apply(rpmat,1,function(x)sample(x,psample))	)		#randomize plot sequence orders, taking "psample" plots
mc.mat <- array(0,dim=c(psample,wide,reps))				#monte carlo matrix: psamples (eg 20) x #allspecies (178) x #reps (eg 100)
for(i in 1:reps) {
	mc.mat[,,i] <- sp.metacom[rpmat[i,],]
}

#-----------
colsum <- apply(mc.mat,c(2,3),sum)		#sum columns of each rep, reps become columns
colsum[colsum>0] <- 1				#convert >0 values to ones
rich.vec <- colSums(colsum)-1			#vector of # cooccurrences for each rep
mc.mat[mc.mat>0] <- 1				#convert species numbers to ones
rowsum <- apply(mc.mat,c(1,3),sum)		#sum rows of each rep, reps become columns
Walpha.vec <- colMeans(rowsum)			#vector of "avg local richness" (Whittaker's alpha) for each rep
if (version == 'multiplicative') Wbeta.vec <- rich.vec/Walpha.vec else Wbeta.vec <- rich.vec-Walpha.vec # here is the option to choose between multiplicative Zeleny's algorithm and additive Fridley's algorithm
GS[sp] <- mean(Wbeta.vec)			#mean Whittaker beta value for all reps (= THETA: G-S metric)
GS.sd[sp] <- sd(Wbeta.vec)			#s.d. of above
meanco[sp] <- mean(rich.vec)			#mean # cooccurrences in "psample" plots
meanco.sd[sp] <- sd(rich.vec)		#s.d. of above

sci.name[sp] <- as.character(SppMat[,2][SppMat[,1]==Species[sp]])	#scientific name
local.avgS[sp] <- mean(rowSums(sp.metacom))				#approximate mean local richness
occur.vec <- colSums(sp.metacom)
tot.cooccur[sp] <- length(occur.vec[occur.vec>0])-1			#total number of species co-occurrences
occur.freq[sp] <- Np							#total number of plots
}

close (win.pb)

#Output matrix
meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
theta.out <- data.frame(sci.name,local.avgS,tot.cooccur,occur.freq,meanco,meanco.sd,meanco.u,meanco.l,GS,GS.sd)
tkdestroy (base)
write.csv2 (file = 'theta_out.csv', theta.out)

spd <- matrix (unlist (strsplit (as.character(theta.out[,1]), split = '_')), ncol = 2, byrow = T)
spaces <- 50-unlist (lapply ((strsplit(spd[,1], split = '' )), length))
rep.spaces <- function (times) paste (rep (' ', times), collapse = '')
write.table (file = 'theta_import.species.data.txt', paste (spd[,1], unlist (apply (as.matrix (spaces), 1, rep.spaces)), theta.out[,9]), quote = F, row.names = F, col.names = F)

tkmessageBox (type = 'ok', message = paste ('Result files were saved into', getwd (), '\n\nYou can use the file theta_import.species.data.txt to import the species theta values directly to JUICE (File > Import > Species Data).'))
end.end <- T
  }
plot.new ()
}


library (tcltk)

end.end <- F

input.matrix <- convert.table.2c (JUICE.table)


GSmethods <- c ('Multiplicative beta diversity (Zeleny 2009)', 'Additive beta diversity (Fridley et al. 2007)')
base <- tktoplevel()
tkwm.title(base, "Calculation of generalists-specialists according to Fridley et al. 2007")

  spec.frm <- tkframe(base,borderwidth=2)
  frame.a <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.b <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.c <- tkframe(spec.frm, relief="groove", borderwidth=2)
  frame.d <- tkframe(spec.frm, relief="groove", borderwidth=2)
  
  GSmethod <- tclVar ('multiplicative')
  label.radio <- tklabel (frame.a, text = 'Which beta diversity algorithm to use?')
  radio1 <- tkradiobutton (frame.a, text = GSmethods[1], value = 'multiplicative', variable = GSmethod)
  radio2 <- tkradiobutton (frame.a, text = GSmethods[2], value = 'additive', variable = GSmethod)
#  tkselect (radi)
  
  tk.psample <- tclVar (5)
  tk.reps <- tclVar (10)
  label.entry1 <- tklabel (frame.b, text = 'Minimal frequency of species')
  entry1 <- tkentry (frame.b, width = 5, textvariable = tk.psample)
  
  label.entry2 <- tklabel (frame.c, text = 'Number of random subsamples')
  entry2 <- tkentry (frame.c, width = 5, textvariable = tk.reps)
  
  button1 <- tkbutton (frame.d, text = 'Calculate', command = function () end.end <- GS.function (input.matrix, psample = as.numeric (tclvalue (tk.psample)), reps = as.numeric (tkget (entry2)), version = as.character (tclvalue (GSmethod))))
  
  tkpack (label.radio, radio1, radio2)
  tkpack (label.entry1, entry1)
  tkpack (label.entry2, entry2)
  tkpack (button1)
  
  tkpack (frame.a, side = 'top', ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = 'w')
  tkpack (frame.b, frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10)
  tkpack (frame.d, side = 'bottom', pady = 10, padx = 10)
  tkpack (spec.frm)
  repeat  if (dev.cur () > 1) break
  dev.off ()