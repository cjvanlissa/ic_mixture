#clear console and env
rm(list=ls(all.names = T))
cat("\014")

#load necessary packages
library(tidyLPA)
library(MplusAutomation)
library(sn)

# Load simulation functions from source -----------------------------------
# source('functions.R')

# set conditions for simulation
hyper_parameters <- list(
  reps = 1:100,
  es = c(1, 2.75, 3.88), # Corresponding to entropy of .24, .7 and .9 for high N
  N = c(40, 100, 300),
  maxK = 3,
  dist = c("norm", "sn", "likert")
)


# Create hypergrid with simulation parameters and save it as .RData file extension
summarydata <- expand.grid(hyper_parameters, stringsAsFactors = FALSE)
set.seed(5738)
summarydata$seed <- sample(1:.Machine$integer.max, nrow(summarydata))
summarydata$rownum <- 1:nrow(summarydata)
saveRDS(summarydata, file = "summarydata.RData")
# summarydata<-readRDS("summarydata.RData")

# prepare parallel processing
library(doSNOW)
nclust <- parallel::detectCores() 
cl <- makeCluster(nclust) 
registerDoSNOW(cl) 

# add progression bar
pb <- txtProgressBar(min = 0, max = nrow(summarydata), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# run simulation
tab <- foreach(rownum = 1:nrow(summarydata), .options.snow = opts, .packages = c("MplusAutomation", "bain", "mvtnorm"), .combine = rbind) %dopar% {
  # Set seed
  attach(summarydata[rownum, ])
  set.seed(seed)
  filnm <- paste0("tmp", Sys.getpid())

  # create dataframe for each group
  dataframe <- data.frame(x = switch(dist,
                norm = c(rnorm(floor(N/2)), rnorm(ceiling(N/2), mean = es)),
                sn = c(sn::rsn(floor(N/2), alpha = 2),
                       sn::rsn(ceiling(N/2), xi = es, alpha = 2)),
                likert = as.integer(cut(c(
                  rnorm(floor(N/2)), rnorm(ceiling(N/2), mean = es)), 5))))

  res <- createMixtures(classes = 1:maxK,
                        filename_stem = filnm,
                        rdata = dataframe, run = 1L)
  out <- mixtureSummaryTable(res)
  df <- out$Classes*2 # Je schreef: Hier komt nog de vector classificatie probabilities bij; dit is K-1 (bij K classen). NB Hier klopt 2*K nu eerst wel volgens mij en dan staat het nu goed, maar wilde het wel melden...
  ll <- sapply(res, function(i){i$results$summaries$LL})
  caic <- -2 * ll + N*(N+df) / (N-df-2)
  aicc <- -2 * ll + N*(df+1) / (N-df-2) # It was: ((df-2)+.001) # Why do you use .001 and why this expression?
  out <- data.frame(rownum = rownum, out[c("Classes", "AIC", "BIC", "aBIC", "Entropy", "T11_VLMR_PValue", 
                   "T11_LMR_PValue", "BLRT_PValue", "min_N", "max_N", "min_prob", 
                   "max_prob")
  ], df = df, ll = ll, caic = caic, aicc = aicc)
  write.table(x = out, file = sprintf("./results/results%d.txt" , Sys.getpid()), 
              sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

#Close cluster
stopCluster(cl)


# End of simulation -------------------------------------------------------
stop("End of simulation")
