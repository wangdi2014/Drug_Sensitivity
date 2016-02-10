'extract.feas.betas.enetobj' <- function(drug.n = 1) {
 browser()
	data.types = c('Gene', 'RNASeq', 'RPPA', 'CNV.GeneLevel', 'Exome.RNA.calls.pvfilter', 'Methylation.pvfilter')
	feas.list = list('matrix', length(data.types))
	max.fea.length = 0
	for(i in 1:length(data.types)) {
		cur.data.type = data.types[i]
		cur.data.type.enetobj.file = paste(cur.data.type, '.elasticnet.obj.drug', drug.n, '.Rdata', sep = '')
		load(cur.data.type.enetobj.file)
		actions = elasticnet.obj$actions
		betas = elasticnet.obj$beta.pure
		steps = NULL
		for(j in 1:length(actions))
			if(!is.null(names(actions[[j]])))
				steps = c(steps, names(actions[[j]]))
		steps = unique(steps)
		idx = match(steps, colnames(betas))
		steps.betas = betas[nrow(betas),idx]
		cur.fea.matrix = cbind(steps, steps.betas)
		colnames(cur.fea.matrix) = c(paste(data.types[i], '.names', sep = ''), paste(data.types[i], '.betas', sep =''))
		feas.list[[i]] = cur.fea.matrix
		if(nrow(cur.fea.matrix) > max.fea.length)
			max.fea.length = nrow(cur.fea.matrix)
	}

	feas.matrix = matrix(' ', max.fea.length, 2 * length(data.types))
	col.idx = 1:2
	col.names = NULL
	for(i in 1:length(feas.list)){
		feas.matrix[1:nrow(feas.list[[i]]), col.idx] = feas.list[[i]]
		col.names = c(col.names, colnames(feas.list[[i]]))
		col.idx = col.idx + 2
	}
	colnames(feas.matrix) = col.names
	write.csv(feas.matrix, file = paste('elasticnet.obj.target.list.drug.', drug.n, '.csv', sep = ''))	
}
