setwd("C:/Users/aless/Documents/UniTn/SecondSemester/NetworkDataAnalysis/NitrogenStarvation/")

exp <- exprs(gse[[1]])

rscudo_probes <- read.csv('rscudo_probes.txt')
rscudo_agi <- anno$AGI[rownames(exp) %in% rscudo_probes$x]
write.csv(rscudo_agi, file = 'rscudo_agi.txt', quote = F, row.names = F, col.names = F)
# rscudo_sym <- anno$`Gene Symbol`[rownames(exp) %in% rscudo_probes$x]
# write.csv(rscudo_sym, file = 'rscudo_sym.txt', quote = F, row.names = F, col.names = F)

rf_probes <- read.csv('rf_probes.txt')
rf_agi <- anno$AGI[rownames(exp) %in% rf_probes$x]
write.csv(rf_agi, file = 'rf_agi.txt', quote = F, row.names = F, col.names = F)

lasso_probes <- read.csv('lasso_probes.txt')
lasso_agi <- anno$AGI[rownames(exp) %in% lasso_probes$x]
write.csv(lasso_agi, file = 'lasso_agi.txt', quote = F, row.names = F, col.names = F)

rf_rscudo_probes <- read.csv('rf_rscudo_probes.txt')
rf_rscudo_agi <- anno$AGI[rownames(exp) %in% rf_rscudo_probes$x]
write.csv(rf_rscudo_agi, file = 'rf_rscudo_agi.txt', quote = F, row.names = F, col.names = F)

list_agi <- anno$AGI
write.csv(list_agi, file = 'list_agi.txt', quote = F, row.names = F, col.names = F)
