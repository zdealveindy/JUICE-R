GS.function <- function (input.matrix, reps, psample, version = 'multiplicative') # the function calculating theta value of species specialization
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

calc.spp <- function (sp)
{
#setWinProgressBar (win.pb, sp, label = paste ('Species no. ', sp))
temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == labels (select.spp[sp])]>0,]
temp.matrix <- temp.matrix[,colSums (temp.matrix > 0) > 0]
rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = dim (temp.matrix)[1], byrow = F)
sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))
mc.mat <- array(0,dim=c(psample,dim (temp.matrix)[2],reps))
for(i in 1:reps) {
	mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
}
total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
mean.alpha <- colMeans (apply (mc.mat, c(1,3), sum))
if (version == 'multiplicative') Wbeta.vec <- total.rich/mean.alpha else Wbeta.vec <- total.rich-mean.alpha # here is the option to choose between multiplicative Zeleny's algorithm and additive Fridley's algorithm
GS <- mean(Wbeta.vec)			#mean Whittaker beta value for all reps (= THETA: G-S metric)
GS.sd <- sd(Wbeta.vec)			#s.d. of above
meanco <- mean(total.rich)			#mean # cooccurrences in "psample" plots
meanco.sd <- sd(total.rich)		#s.d. of above
sci.name <- labels (select.spp[sp])	#scientific name
local.avgS <- mean(mean.alpha)				#approximate mean local richness
occur.freq <- as.vector (select.spp[sp])							#total number of plots
c (sci.name,local.avgS,occur.freq,meanco,meanco.sd,GS,GS.sd)
}


input.matrix[input.matrix > 0] <- 1

Nplots <- dim (input.matrix)[1]

plots.per.spp <- colSums (input.matrix)
select.spp <- plots.per.spp[plots.per.spp >= psample]

Nspp <- length (select.spp)

win.pb <- winProgressBar (title = 'Calculation progress', label = paste ('Species no. ', 1), min = 1, max = Nspp, initial = 0, width = 300) 


GS <- foreach (sp=1:Nspp) %dopar% calc.spp(sp)
close (win.pb)

#Output matrix
GS <- t(as.data.frame (GS))
rownames (GS) <- GS[,1]
GS <- GS[,-1]
mode (GS) <- 'numeric'
print (GS)
meanco.u <- qnorm(.975,mean=GS[,3],sd=GS[,4])			#97.5% confidence limit
meanco.l <- qnorm(.025,mean=GS[,3],sd=GS[,4])			#2.5% confidence limit
theta.out <- cbind (GS[,1:4], meanco.u, meanco.l, GS[,5:6])
colnames (theta.out) <- c('local.avgS','occur.freq','meanco','meanco.sd','meanco.u','meanco.l','GS','GS.sd')
tkdestroy (base)
write.csv2 (file = 'theta_out.csv', theta.out)

#spd <- matrix (unlist (strsplit (as.character(theta.out[,1]), split = '_')), ncol = 2, byrow = T)
#spaces <- 50-unlist (lapply ((strsplit(spd[,1], split = '' )), length))
#rep.spaces <- function (times) paste (rep (' ', times), collapse = '')
#write.table (file = 'theta_import.species.data.txt', paste (spd[,1], unlist (apply (as.matrix (spaces), 1, rep#.spaces)), theta.out[,8]), quote = F, row.names = F, col.names = F)

tkmessageBox (type = 'ok', message = paste ('Result files were saved into', getwd (), '\n\nYou can use the file theta_import.species.data.txt to import the species theta values directly to JUICE (File > Import > Species Data).'))
end.end <- T
  }
plot.new ()
}


library (tcltk)

end.end <- F

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
  
  tk.psample <- tclVar (5)
  tk.reps <- tclVar (10)
  label.entry1 <- tklabel (frame.b, text = 'Minimal frequency of species')
  entry1 <- tkentry (frame.b, width = 5, textvariable = tk.psample)
  
  label.entry2 <- tklabel (frame.c, text = 'Number of random subsamples')
  entry2 <- tkentry (frame.c, width = 5, textvariable = tk.reps)
  
  button1 <- tkbutton (frame.d, text = 'Calculate', command = function () end.end <- GS.function (as.matrix (JUICE.table), psample = as.numeric (tclvalue (tk.psample)), reps = as.numeric (tkget (entry2)), version = as.character (tclvalue (GSmethod))))
  
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
  stopWorkers (workers)
