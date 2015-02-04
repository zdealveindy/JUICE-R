# JUICE-R JUICE.table JUICE.short.head JUICE.species.data R.exe no.transf

#setwd("C:/Program Files/R/R-2.15.0/bin/x64/")
#JUICE.table <- read.table ('table.txt', sep = '\t', 
#                check.names = F, head = T, row.names = 1)
#JUICE.short.head <- read.table('ShortHead.txt', sep='\t',
#                head = T)

require (tcltk)
beals.spec.pool <- function (x, grouping, include = TRUE) 
{
	win.pb2 <- winProgressBar (title = 'Beals (the whole matrix)', label = 'start', min = 1, max = nrow (x)+2, initial = 0, width = 300)
	incSp <- include
	x <- as.matrix (x)
	x <- ifelse(x > 0, 1, 0)  # all data are transformed into p/a data
	target <- x[grouping == 1,]  # selects plots for which sp. pool should be calculated (grouping == 1)

	setWinProgressBar (win.pb2, 1, label = 'Crossprod')
	M <- crossprod(x, x)  # calculates number of joint co-occurrences of species
	C <- diag(M)	# diagonal contains number of occurrences for given species

	setWinProgressBar (win.pb2, 1, label = 'First sweep')
	M <- sweep(M, 2, replace(C, C == 0, 1), "/")  # this calculates p[ij]
	if (!incSp) 
		for (i in 1:ncol(x)) M[i, i] <- 0  # if the target species is not included

# the following part calculates beals smoothed values for all data - the only
# purpose of this is to calculate cut-off threshold values for individual species 
	x.S <- rowSums (x)
	x.b <- x

	for (i in 1:nrow (x)) 
		{
		setWinProgressBar (win.pb2, i+1, label = i)
		x.b [i, ] <- rowSums (sweep (M, 2, x[i, ], "*"))
		}

	x.SM<- rep(x.S, ncol(x))	
	if (!incSp) x.SM <- x.SM - x
	x.b <- x.b/replace(x.SM, x.SM == 0, 1)
	close (win.pb2)

	win.pb2 <- winProgressBar (title = 'Beals (target matrix)', label = 'start', min = 1, max = nrow (target)+2, initial = 0, width = 300)
	setWinProgressBar (win.pb2, 1, label = 'Beals threshold')

	outlier.fun <- function (y) 
		{
		beals.occ <- y[y>0]
		stats.temp <- fivenum (beals.occ)
		iqr <- diff (stats.temp [c(2,4)])
		min (beals.occ[!(beals.occ < stats.temp[2] - 1.5* iqr)])
		}

	beals.thresh <- apply (x.b*x, 2, FUN = function (y) outlier.fun(y))  # non-outlier species beals threshold

	S <- rowSums(target)
	b <- target
	for (i in 1:nrow(target))
		{
		setWinProgressBar (win.pb2, i+1, label = i)
		b[i, ] <- rowSums(sweep(M, 2, target[i, ], "*"))
		}
	SM <- rep(S, ncol(target))
	if (!incSp) SM <- SM - target
	b <- b/replace(SM, SM == 0, 1)

	setWinProgressBar (win.pb2, nrow (target)+1, label = 'Calculating species pool matrix')
	spec.pool <- t(apply (b, 1, FUN = function (x)  (x >= beals.thresh)))
	close (win.pb2)
	class (spec.pool) <- 'numeric'
	spec.pool
}

cancel <- tclVar (1)

x <- JUICE.table
sp.p <- beals.spec.pool (x, grouping = JUICE.short.head$Group_Number)
write.table (sp.p, file = 'sp.pool.txt', sep = ';')
save (sp.p, file = 'sp.pool.r')

t.sp.p <- t (sp.p)
file <- 'sp.pool.imp.juice.txt'

cat ('Table from  relevés of the file: result.wct\n', file = file)
cat ('\n', file = file, append = T)
cat (paste ('Number of relevés:', nrow (sp.p), '\n\n'), file = file, append = T)
cat (c(';', colnames (t.sp.p), '\n'), file = file, append = T, sep = ';' )
for (ro in seq (nrow (t.sp.p)))
	cat (c(as.character (JUICE.species.data$Species_Name[ro]), JUICE.species.data$Layer[ro], t.sp.p[ro,], '\n'), file = file, append = T, sep = ';')

cancel <- tkmessageBox (type = 'ok', message = paste ('Result files were saved into', getwd (), '\n\n
The file sp.pool.txt contains the matrix with species pool data (matrix species x samples, presence of species in species pool). You can use the file sp.pool.imp.juice.txt to import the whole matrix with species pool data into JUICE (File > Import > Table > From Spreadsheet File (e.g. Excell table)) and keep clicking Continue - Cover values change into Percentage values).'))
