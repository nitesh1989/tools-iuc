######################
# setup R error handling to go to stderr
######################
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
    'quiet', 'q', 2, "logical",
    'help' , 'h', 0, "logical",
    ,byrow=TRUE, ncol=4)
opt <- getopt(spec)

# If help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}


## Set verbose mode
verbose = if(is.null(opt$quiet)){TRUE}else{FALSE}
if(verbose){
    cat("Verbose mode is ON\n\n")
}


# Enforce the following required arguments
if (is.null(opt$preprocess)) {
    cat("'--preprocess' is required\n")
    q(status=1)
}

# Load required libraries
suppressPackageStartupMessages({
    library("minfi")
    library("ELMER")
    library("doParallel")
})

################################

################################

# Get 450k data
# 1. Needs to be a matrix
# 2. rownames need to be cg-ids
# 3. Remove rownames from the column
# 4. Same column names as RNA data

Meth = read.csv(opt$Meth, header=TRUE)
#colnames(Meth) = gsub(fixed=TRUE,pattern=".",replacement="-",x=colnames(Meth))
#rownames(Meth) = Meth$X
#drops <- c("X");
#Meth = Meth[,!(names(Meth) %in% drops)]
# Convert to matrix
Meth = as.matrix(Meth)

# Get RNA data
# 1. Needs to be a matrix
# 2. rownames need to be Gene ids
# 3. Remove rownames from the column
# 4. Same column names as 450k Data

RNA = read.csv(opt$GeneExp,header=TRUE)
#colnames(RNA) = gsub(fixed=TRUE,pattern=".",replacement="-",x=colnames(RNA))
#drops <- c("X"); 
#RNA = RNA[,!(names(RNA) %in% drops)]
# Convert to Matrix
RNA = as.matrix(RNA)


# TODO: What happens to extra columns in MEth, not in RNA, and vice versa

# TODO: Get Sample information
# Define a CSV file with file names, Sample,Normal information



# mee@sample$TN = set to CSV file





# Data check before proceeding
stopifnot(colnames(RNA) == colnames(Meth))
stopifnot(is(RNA,"matrix"))
stopifnot(is(Meth,"matrix"))


# Get probe annotation
probe <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, what="Locations")
probe <- GRanges(seqnames=probe$chr,ranges=IRanges(probe$pos,width=1,names=rownames(probe)),
                strand=probe$strand,name=rownames(probe))
mee <- fetch.mee(meth=Meth, exp=GeneExp, TCGA=TRUE, probeInfo=probe)

# Get gene information
geneAnnot <- txs()
geneAnnot$GENEID <- paste0("ID",geneAnnot$GENEID)
geneInfo <- promoters(geneAnnot,upstream = 0, downstream = 0)


# Fetch Mee
# TODO: Recommend using the ELMER TCGA pipeline for TCGA data
mee <- fetch.mee(meth=Meth, exp=RNA, TCGA=FALSE, probeInfo=probe, geneInfo=geneInfo)


# Stop if MEE doesn't work
stopifnot(!is.null(mee))



### Results to show after MEE fetch

#Get distal enhancer probes that are 2kb away from TSS and overlap with REMC and FANTOM5
#enhancers on chromosome 1
Probe <- get.feature.probe(rm.chr=paste0("chr",c(2:22,"X","Y")))


# Identify differentially methylated enchancer probes 

# TODO: Get direction of Hypo or Hyper methylation
sig.diff <- get.diff.meth(mee, cores=detectCores()/2, dir.out ="./Result",diff.dir="hypo", pvalue = 0.01,save=TRUE)

sig.probes <- sig.diff[,1]
## Collect nearby 20 genes for Sig.probes
nearGenes <- GetNearGenes(TRange=getProbeInfo(mee,probe=Sig.probes),geneAnnot=getGeneInfo(mee),cores=detectCores()/2)

#Both files contain four columns: probe, GeneID, Symbol, Distance, Sides, Raw.p, Pe.
#1. Probe: the name of probes.
#2. GeneID and Symbol is for the genes which are linked to the probe.
#3. Distance: the distance between the probe and the gene.
#4. Sides: right (R) side or left (L) side of probe and the rank based on distance. For example, L3 means the gene is the
#number 3 closest gene from the left side of the probe.
#5. Raw.p: P value from the Mann-Whitney U test for each pair.
#6. Pe: the empirical P value for each pair.

# Identify significant probe-gene pairs
Hypo.pair <-get.pair(mee=mee,probes=Sig.probes,nearGenes=nearGenes,
                     permu.dir="./ELMER.example/Result/LUSC/permu",permu.size=300,Pe = 0.01,
                     dir.out="Result",cores=detectCores()/2,label= "hypo",save = TRUE)






