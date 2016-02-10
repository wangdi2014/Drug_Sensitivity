'rank.all.drugs' <- function(){
#  browser()
	all.drugs.res = read.csv('/media/OS/NCIDREAM/DREAM7_DrugSensitivity1_Drug_Response_All.csv', check.names=FALSE)
	row.names=as.character(all.drugs.res[,1])
	
	all.testing.drugs = as.matrix(all.drugs.res[,2:(2+18-1)])
	rownames(all.testing.drugs) = row.names
	storage.mode(all.testing.drugs) = 'double'
	
	all.drugs.res = as.matrix(all.drugs.res[, 2:ncol(all.drugs.res)])
	rownames(all.drugs.res) = row.names
	storage.mode(all.drugs.res) = 'double'
	
	all.testing.drugs.ranking = all.testing.drugs
	for(i in 1:nrow(all.testing.drugs)){
		cur.drug = all.testing.drugs[i,]
		all.testing.drugs.ranking[i,] = rank(-cur.drug, na.last = 'keep', ties.method = 'random')
	}
	write.csv(all.testing.drugs.ranking, file='all.testing.drugs.ranking.csv')
	write.csv(all.testing.drugs, file='all.testing.drugs.res.csv')
	
	all.drugs.res.ranking = all.drugs.res
	for(i in 1:nrow(all.drugs.res)){
		cur.drug = all.drugs.res[i,]
		all.drugs.res.ranking[i,] = rank(-cur.drug, na.last = 'keep', ties.method = 'random')
	}
	write.csv(all.drugs.res.ranking, file='all.drugs.res.ranking.csv')
	write.csv(all.drugs.res, file='all.drugs.res.csv')

}
