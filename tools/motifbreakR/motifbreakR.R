# setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

## Load required libraries
suppressPackageStartupMessages({
    library("getopt")
    library("motifbreakR")
    library("BSgenome")
    library("MotifDb")
})

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)s

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    'quiet', 'q', 2, "logical",
    'help' , 'h', 0, "logical",
    ,byrow=TRUE, ncol=4)
opt <- getopt(spec)

# print list of snps
cat("verbose = ", opt$quiet,"\n")
cat("preprocess = ",opt$preprocess,"\n")


# 1. Get list of snps
