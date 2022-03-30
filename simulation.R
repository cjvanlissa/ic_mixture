#clear console and env
rm(list=ls(all.names = T))
cat("\014")

#load necessary packages
library(tidyLPA)
library(MplusAutomation)
library(sn)
library(doSNOW)
library(bain)

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

if(dir.exists("results")){
  unlink("results", recursive = T)
}
dir.create("results")

# prepare parallel processing
nclust <- parallel::detectCores()
cl <- makeCluster(nclust) 
registerDoSNOW(cl) 

# add progression bar
pb <- txtProgressBar(min = 0, max = nrow(summarydata), style = 3)
opts <- list(progress = function(n) setTxtProgressBar(pb, n))

# run simulation
tab <- foreach(rownum = 1:nrow(summarydata), .options.snow = opts, .packages = c("MplusAutomation", "bain"), .combine = rbind) %dopar% {
  # Set seed
  attach(summarydata[rownum, ])
  set.seed(seed)
  filnm <- paste0("tmp", Sys.getpid())

  # create dataframe for each group
  dat <- data.frame(x = switch(dist,
                norm = c(rnorm(floor(N/2)), rnorm(ceiling(N/2), mean = es)),
                sn = c(sn::rsn(floor(N/2), alpha = 2),
                       sn::rsn(ceiling(N/2), xi = es, alpha = 2)),
                likert = as.integer(cut(c(
                  rnorm(floor(N/2)), rnorm(ceiling(N/2), mean = es)), 5))))

  res <- createMixtures(classes = 1:maxK,
                        filename_stem = filnm,
                        rdata = dat, run = 1L)
  out <- mixtureSummaryTable(res)
  df <- out$Classes*2 #
  tmp <- rep(NA, maxK)
  ll <- sapply(res, function(i){i$results$summaries$LL})
  ns <- t(sapply(res, function(i){
    tmp[1:length(i$results$class_counts$posteriorProb$count)] <- i$results$class_counts$posteriorProb$count
    tmp
    }))
  
  caic <- -2 * ll + N*(N+df) / (N-df-2)
  aicc <- -2 * ll + N*(df+1) / (N-df-2) # CJ: Aaah ik had een copy-paste fout gemaakt uit jouw email, waardoor er stond / (df-2). Dat werkte niet
  out <- data.frame(rownum = rownum, out[c("Classes", "AIC", "BIC", "aBIC", "Entropy", "T11_VLMR_PValue", 
                   "T11_LMR_PValue", "BLRT_PValue", "min_N", "max_N", "min_prob", 
                   "max_prob")], df = df, ll = ll, caic = caic, aicc = aicc, ns)
  write.table(x = out, file = sprintf("./results/results%d.txt" , Sys.getpid()), 
              sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
}

#Close cluster
stopCluster(cl)
stop("End of simulation")

f <- c(list.files(pattern = "^tmp"), list.files(pattern = "^data_tmp"))
file.remove(f)

# Read files --------------------------------------------------------------


library(data.table)

res <- as.data.table(readRDS(list.files(pattern = "summarydata")))

vars <- c("rownum", "Classes", "AIC", "BIC", "aBIC", "Entropy", "T11_VLMR_PValue", 
  "T11_LMR_PValue", "BLRT_PValue", "min_N", "max_N", "min_prob", 
  "max_prob", "df", "ll", "caic", "aicc", paste0("N", 1:maxK))

f <- list.files("results", full.names = TRUE)
tab <- lapply(f, fread, header = F)
tab <- rbindlist(tab)
setorderv(tab, cols = c("V1", "V2"), order=1L, na.last=FALSE)
if(!(tab$V1[1] == 1 & tail(tab$V1, 1) == nrow(res) & length(unique(tab$V1)) == nrow(res))){
  c(1:nrow(res))[!c(1:nrow(res)) %in% unique(tab$V1)]
  stop()
}
names(tab) <- vars
#tab[, "rownum" := NULL]
merged <- merge(res, tab, by = "rownum")

fwrite(merged, paste0("sim_results_", Sys.Date(), ".csv"))
saveRDS(merged, paste0("sim_results_", Sys.Date(), ".RData"))