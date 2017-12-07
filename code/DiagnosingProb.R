library(TMB)


for (i in 1:18) {
load(file.path(paste('2017-04-04_RUN_DIAG_',i,'_250knots/sdreport',i,'.RData',sep = '')));   assign(paste('SD',i,sep =''), SD)
}

Qs <- list(SD1,SD2,SD3,SD4,SD5,SD6,SD7,SD8,SD9,SD10,SD11,SD12,SD13,SD14,SD15,SD16,SD17,SD18)

save(Qs, file = 'QEstimates.RData')



SD1 # cod adu,  ok
SD2 # cod juv,  bad
SD3 # meg adu,  ok
SD4 # meg juv,  ok
SD5 # bud adu,  ok
SD6 # bud juv,  bad 
SD7 # pis adu,  ok
SD8 # pis juv,  ok
SD9 # had adu,  ok
SD10 # had juv,  ok
SD11 # whg adu, ok 
SD12 # whg juv, ok 
SD13 # hke adu, ok
SD14 # hke juv,  ok
SD15 # ple adu, ok 
SD16 # ple juv, ok 
SD17 # sol adu, ok 
SD18 # Sol juv, ok 
