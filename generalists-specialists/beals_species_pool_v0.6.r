# JUICE-R JUICE.table JUICE.short.head JUICE.species.data R.exe no.transf

# the code above this line is JUICE-R code - it's important to keep at the place if you intend to use this script from JUICE program (Tichy 2002)



### Script for calculation of species pool, using Beals smoothing method 
### Author: David Zeleny (zeleny@sci.muni.cz)
### 23-01-2013
### the function 'beals.spec.pool' is based on function 'beals' (library vegan), written by Miguel De Cáceres and Jari Oksanen



### The whole script can be applied 'as is' from JUICE as a JUICE-R function (see http://www.bit.ly/juice-r). Alternatively, functions 'beals.spec.pool' and
### 'beals.spec.pool.2' can be used independently directly in R. Optionally, this function generates several files, which are stored in working directory:
### 'spec.pool.txt' - tab delimited matrix with species x samples table, presence absence data, species pool of individual samples of target dataset
### 'spec.poo.imp.juice.txt' - file formated for direct import into JUICE program (in JUICE menu go File > Import > Table > From Spreadsheet File (e.g. Excell table)
### and keep clicking Continue - Cover values change into Percentage values))

### NOTE ABOUT USE FROM JUICE - if you want to calculate species pool using target-background method, you need to separate (in JUICE) the matrix into 2 parts - first will become
### target matrix, second will become background matrix. If there is no separation, script will calculate target-target method. If there are more than 2 parts,
### script will ignore third and higher and calculate target-background method from first two

### Function 'beals.spec.pool' is for calculating species pool using method target - background (species pool for target samples is calculated based on
### co-ooccurrence data from background dataset)
### Function 'beals.spec.pool.2'is for calculating species pool using method target - target (co-occurrence data are calculated also from target table, but after
### removing the sample for which the species pool is calculated). This function is much slower (basically it must recalculate the whole matrix many times), so
### there is also parallel version (setup argument parallel = T and no.cores = 8 (if you have 8 cores on your computer)). Parallel calculation require version of
### R >= 2.15.0

### Explanation of arguments of 'beals.spec.pool' function:
###
###	x - species-samples matrix, merging together target releves (sorted at the beginning) and background releves (appended after)
###	(target releves are those the species pool should be calculated for, bacground releves are for calculating co-occurrences)
###
###	grouping - vector of the same length as nrow (x) containing numbers 1 and 2, indicating which releves are target (1) and background (2)
###
###	species.data - matrix or data frame with two columns, first containing species names and second with layer information (1 - trees, 4 - shrubs, 6 - herbs,
###	7 - juveniles, 9 - mosses); species names MUST be sorted in the same way as species in matrix x; this variable is optional - if provided, the file intended to
###	be imported into JUICE is generated ('sp.pool.imp.juice.txt')
###
###	include - should be species, for which probability of occurrence is calculated, included in the releve? (default option is FALSE, which is consistent with 
###	modification of Herben & Munzbergova)
###
###	cooc.matrix - which releves should be used as background to calculated co-occurrences? Default is 'background', i.e. only the releves considered
###		as background releves will be used (indicated by number 2 in 'grouping' vector). Alternative is 'all', meaning that co-occurrences will be
###		calculated from all releves (i.e. the whole table containing also target releves). 'background' method is more straightforward, but it has
###		dissadvantage - the species, which occur in target releves (and are missing in background releves) cannot be included in analysis and will be ignored.
###
###	write.to.txt - should the 'sp.pool.txt' file with species pool data (matrix of species x samples) be generated? Default FALSE
### 	verbal - should be the progress bars and messages shown? This is useful expecially for large datasets, to see the progress of calculation. Default = TRUE
###
### parallel (only in beals.spec.pool.2) - if you have multiple cores on your computer, you may run the calculation in parallel, which will significantly enhance
###   the speed; default is FALSE
###
### no.cores (only in beals.spec.pool.2) - setup number of cores you have on your computer. If the matrix of samples is large, for parallel computing you need also
###   large RAM

beals.spec.pool <- function (x, grouping, species.data = NULL, include = FALSE, cooc.matrix = 'background', write.to.txt = FALSE, verbal = T) 
{
	require (tcltk)

	outlier.fun <- function (y) 
		{
		beals.occ <- y[y>0]
		stats.temp <- fivenum (beals.occ)
		iqr <- diff (stats.temp [c(2,4)])
		min (beals.occ[!(beals.occ < (stats.temp[2] - 1.5* iqr))])
		}

	incSp <- include
	cooc.matrix <- c('background', 'all')[pmatch (cooc.matrix, table = c('background', 'all'))]
	x <- as.matrix (x)
	x <- ifelse(x > 0, 1, 0)  # all data are transformed into p/a data

	target <- x[grouping == 1,,drop = F]  # selects plots for which sp. pool should be calculated (grouping == 1)
	background <- if (cooc.matrix == 'background') x[grouping == 2,, drop = F] else x

	if (verbal) win.pb2 <- winProgressBar (title = 'Beals (background matrix)', label = 'start', min = 1, max = nrow (background)+2, initial = 0, width = 300)
	if (verbal) setWinProgressBar (win.pb2, 1, label = 'Crossprod')
	M <- crossprod(background, background)  # calculates number of joint co-occurrences of species
	C <- diag(M)	# diagonal contains number of occurrences for given species

	if (verbal) setWinProgressBar (win.pb2, 1, label = 'First sweep')
	M <- sweep(M, 2, replace(C, C == 0, 1), "/")  # this calculates p[ij]
	if (!incSp) 
		for (i in 1:ncol(background)) M[i, i] <- 0  # if the target species is not included

# the following part calculates beals smoothed values for background data - the only purpose of this is to calculate cut-off threshold values for individual species 
# which is recommended by Herben & Munzbergova

	x.S <- rowSums (background)
	x.b <- background

	for (i in 1:nrow (background)) 
		{
		if (verbal) setWinProgressBar (win.pb2, i+1, label = i)
		x.b [i, ] <- rowSums (sweep (M, 2, background[i, ], "*"))
		}

	x.SM <- rep(x.S, ncol(background))	
	if (!incSp) x.SM <- x.SM - background
	x.b <- x.b/replace(x.SM, x.SM == 0, 1)
	if (verbal) close (win.pb2)

	if (verbal) win.pb2 <- winProgressBar (title = 'Beals (target matrix)', label = 'start', min = 1, max = nrow (target)+2, initial = 0, width = 300)
	if (verbal) setWinProgressBar (win.pb2, 1, label = 'Beals threshold')

	beals.thresh <- apply (x.b*background, 2, FUN = function (y) outlier.fun(y))  # non-outlier species beals threshold

	S <- rowSums(target)
	b <- target
	for (i in 1:nrow(target))
		{
		if (verbal) setWinProgressBar (win.pb2, i+1, label = i)
		b[i, ] <- rowSums(sweep(M, 2, target[i, ], "*"))
		}
	SM <- rep(S, ncol(target))
	if (!incSp) SM <- SM - target
	b <- b/replace(SM, SM == 0, 1)

	if (verbal) setWinProgressBar (win.pb2, nrow (target)+1, label = 'Calculating species pool matrix')
	spec.pool <- t(apply (b, 1, FUN = function (x)  (x >= beals.thresh)))
	if (verbal) close (win.pb2)
	class (spec.pool) <- 'numeric'

# export of species pool data into txt format (file 'sp.pool.txt', with cells separated by tabulators
	if (write.to.txt) write.table (spec.pool, file = 'sp.pool.txt', sep = '\t')

# export of species pool data into txt format, which could be directly read by JUICE (optional, only if data have been exported from JUICE and are intended to be
# imported into JUICE again; requires JUICE.species.data file to be imported
if(!is.null (species.data))
	{
#	M <- crossprod(background)
	t.sp.p <- t (spec.pool)
	file <- 'sp.pool.imp.juice.txt'

	if (verbal) win.pb3 <- winProgressBar (title = 'Creating JUICE import file', label = 'start', min = 1, max = nrow (t.sp.p), initial = 0, width = 300)
	cat ('Table from  relevés of the file: result.wct\n', file = file)
	cat ('\n', file = file, append = T)
	cat (paste ('Number of relevés:', nrow (spec.pool), '\n\n'), file = file, append = T)
	cat (c(';', colnames (t.sp.p), '\n'), file = file, append = T, sep = ';' )
	for (ro in seq (nrow (t.sp.p)))
		{
		if (verbal) setWinProgressBar (win.pb3, ro, label = ro)
		cat (c(as.character (JUICE.species.data$Species_Name[ro]), JUICE.species.data$Layer[ro], t.sp.p[ro,], '\n'), file = file, append = T, sep = ';')
		}
	if (verbal) close (win.pb3)
	}

# the final message
	if (verbal) cancel <- tkmessageBox (type = 'ok', message = paste ('Result files were saved into', getwd (), '\n\nThe file sp.pool.txt contains the matrix with species pool data (matrix species x samples, presence of species in species pool)',
	if(!is.null (species.data)) '\n\nYou can use the file sp.pool.imp.juice.txt to import the whole matrix with species pool data into JUICE (File > Import > Table > From Spreadsheet File (e.g. Excell table)) and keep clicking Continue - Cover values change into Percentage values).'))

	return (spec.pool)

}

beals.spec.pool.2 <- function (x, species.data = NULL, include = FALSE, cooc.matrix = 'background', write.to.txt = FALSE, verbal = T, parallel = F, no.cores = 4)
{
if (parallel) verbal <- F

if (verbal) win.pb4 <- winProgressBar (title = 'Beals species pool', label = 'start', min = 1, max = nrow (x), initial = 0, width = 300)

if (!parallel)
{
  spec.pool <- x
  for (i in seq (nrow (x)))
	 {
	 if (verbal) setWinProgressBar (win.pb4, i, label = paste ('Sample', i, 'of', nrow (x)))
	 grouping <- rep (2, nrow (x))
	 grouping[i] <- 1
	 spec.pool[i,] <- beals.spec.pool (x, grouping = grouping, species.data = NULL, include = include, cooc.matrix = cooc.matrix, write.to.txt = FALSE, verbal = FALSE)
	}
}

if (parallel)
{
sp <- x
  require (parallel)
  cl <- makeCluster (no.cores)
beals.fun <- function (i, sp, include, cooc.matrix)
  {
	grouping <- rep (2, nrow (sp))
	grouping[i] <- 1
	beals.spec.pool (sp, grouping = grouping, species.data = NULL, include = include, cooc.matrix = cooc.matrix, write.to.txt = FALSE, verbal = FALSE)
   }
#clusterExport (cl, list('x', 'include', 'cooc.matrix', 'beals.spec.pool'))
#clusterExport (cl, 'beals.fun')
clusterExport (cl, 'beals.spec.pool')
spec.pool <- parSapply (cl, seq (nrow (x)), FUN = beals.fun, sp = sp, include = include, cooc.matrix = cooc.matrix)
spec.pool <- as.data.frame (t(spec.pool))
names (spec.pool) <- names (x)
rownames (spec.pool) <- rownames (x)
stopCluster(cl)
}



# export of species pool data into txt format (file 'sp.pool.txt', with cells separated by tabulators
	if (write.to.txt) write.table (spec.pool, file = 'sp.pool.txt', sep = '\t')

# export of species pool data into txt format, which could be directly read by JUICE (optional, only if data have been exported from JUICE and are intended to be
# imported into JUICE again; requires JUICE.species.data file to be imported
if(!is.null (species.data))
	{
#	M <- crossprod(background)
	t.sp.p <- t (spec.pool)
	file <- 'sp.pool.imp.juice.txt'

	if (verbal) win.pb3 <- winProgressBar (title = 'Creating JUICE import file', label = 'start', min = 1, max = nrow (t.sp.p), initial = 0, width = 300)
	cat ('Table from  relevés of the file: result.wct\n', file = file)
	cat ('\n', file = file, append = T)
	cat (paste ('Number of relevés:', nrow (spec.pool), '\n\n'), file = file, append = T)
	cat (c(';', colnames (t.sp.p), '\n'), file = file, append = T, sep = ';' )
	for (ro in seq (nrow (t.sp.p)))
		{
		if (verbal) setWinProgressBar (win.pb3, ro, label = ro)
		cat (c(as.character (JUICE.species.data$Species_Name[ro]), JUICE.species.data$Layer[ro], t.sp.p[ro,], '\n'), file = file, append = T, sep = ';')
		}
	if (verbal) close (win.pb3)
	}

if (verbal) close (win.pb4)
if (verbal) cancel <- tkmessageBox (type = 'ok', message = paste ('Result files were saved into', getwd (), '\n\nThe file sp.pool.txt contains the matrix with species pool data (matrix species x samples, presence of species in species pool)',
if(!is.null (species.data)) '\n\nYou can use the file sp.pool.imp.juice.txt to import the whole matrix with species pool data into JUICE (File > Import > Table > From Spreadsheet File (e.g. Excell table)) and keep clicking Continue - Cover values change into Percentage values).'))

return (spec.pool)
}

# calculation of species pool
calculate.function <- function (parallel, no.cores)
{
if (length (unique (grouping)) == 1)
	sp.p <- beals.spec.pool.2 (x = x, species.data = species.data, include = T, cooc.matrix = 'background', write.to.txt = TRUE, parallel = parallel, no.core = no.cores)
if (length (unique (grouping)) == 2)
	sp.p <- beals.spec.pool (x = x, grouping = grouping, species.data = species.data, include = T, cooc.matrix = 'background', write.to.txt = TRUE )
if (length (unique (grouping)) > 2) {require (tcltk); tkmessageBox (type = 'ok', message = 'Too many groups of samples! You need only one or two!')}
}

### THE FOLLOWING LATER DELETE !!!
# the following is for reading the data (in case that the script is not used within JUICE-R function)
#setwd("C:/Program Files/R/R-2.15.0/bin/x64/")
#JUICE.table <- read.table ('table.txt', sep = '\t', check.names = F, head = T, row.names = 1)
#JUICE.short.head <- read.table('ShortHead.txt', sep='\t',  head = T)
#JUICE.species.data <- read.table ('SpeciesData.txt', sep='\t', head = T, fill = T)

# CALCULATION
# import of data from JUICE
x <- JUICE.table	# matrix of species-sample data
grouping <- JUICE.short.head$Group_Number		# vector of numbers (1,2) denoting target and background releves
species.data <- JUICE.species.data[, c('Species_Name', 'Layer')]

# Tcltk wizard at the beginning
require (tcltk)
cancel <- tclVar (0)
parallel.comp <- tclVar (0)
no.cores <- tclVar (2)
length.target <- sum (grouping == 1)
length.background <- if (length (unique (grouping)) == 1) length.target else sum (grouping == 2)
length.all <- length (grouping)
 
base <- tktoplevel()
tkwm.title(base, 'JUICE-R application' )

spec.frm <- tkframe(base,borderwidth=5)

frame.a <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.b <- tkframe (spec.frm, borderwidth = 0)
frame.c <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
frame.c.1 <- tkframe (frame.c)
tkpack (tklabel (spec.frm, text = "Calculation of species pool based on Beals smoothing"))

tkpack (tklabel (frame.a, text = 'Number of samples used in analysis:'))
tkpack (tklabel (frame.a, text = paste ('target samples: ', length.target, ', background samples: ', length.background, ', all samples: ', length.all, sep = '')))


  parallel.label <- tklabel (frame.c, text = 'Parallel calculation (enable if you have more than one core)')
  parallel.no.cores.label <- tklabel (frame.c.1, text = 'number of cores: ')
  parallel.no.cores.entry <- tkentry (frame.c.1, width = 2, textvariable = no.cores)
  parallel.checkbutton <- tkcheckbutton (frame.c, text = 'use parallel computing', variable = parallel.comp)

tkpack (parallel.no.cores.label, parallel.no.cores.entry, anchor = 'w', side = 'left')
tkpack (parallel.label, parallel.checkbutton, frame.c.1, anchor = 'w')

tkpack (tkbutton (frame.b, text = 'Calculate', command = function () calculate.function (parallel = as.logical (as.numeric (tclvalue (parallel.comp))), no.cores = as.numeric (tclvalue (no.cores)))), side = 'right')



tkpack (frame.a, side = "top", ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w")
tkpack (frame.b, side = 'right', ipady = 2, ipadx = 2, padx = 10, pady = 10, anchor = 'w')
tkpack (frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w")

tkpack (spec.frm)
tkbind (base, "<Destroy>", function() tclvalue(cancel)<-2) 
tkwait.variable (cancel)

