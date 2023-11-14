library(breakpointR)

args = commandArgs(trailingOnly = TRUE)
bamdir = args[1]
outputdir = args[2]
chromosomes = args[3]
pairedEndReads = F

if (chromosomes == 'all'){
    chromosomes = paste0('chr',c(1:22,'X','Y'))
}

breakpointr(inputfolder = bamdir,
            outputfolder = outputdir,
            chromosomes = chromosomes,
            pairedEndReads = pairedEndReads)
