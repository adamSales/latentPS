pload('~/Box Sync/CT/data/RANDstudyData/MSdata.RData')
datMS <- dat

source('R/prelimStan.r')


totDatMS <- dataPrep(datMS,advanceOrig)
datMS <- totDatMS$dat
advanceMS <- totDatMS$advance
rm(totDat);gc()

sdatMS <- makeStanDat(datMS,advanceMS)


msMod <- stan('R/psmod.stan',data=sdatMS,iter=10000,chains=8)

save(msMod,sdatMS,file="output/msMod.RData")

rm(msMod); gc()
