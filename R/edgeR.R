library(edgeR)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--rawcount", help="the raw count file")
argv <- add_argument(argv,"--condition", help="the condition file")
argv <- add_argument(argv,"--comparename", help="the compare name")
argv <- add_argument(argv,"--design", help="the design of experiment")
argv <- add_argument(argv,"--foldchange", help="the foldchange value",type='numeric')
argv <- add_argument(argv,"--pvalue", help="the p value",type='numeric')
argv <- add_argument(argv,"--padj", help="the p adjust value",type='numeric')
argv <- add_argument(argv,"--genename", help="the gene name file")
argv <- add_argument(argv,"--outdir", help="the directory of output file")
argv <- parse_args(argv)

rawcount <- argv$rawcount
condition <- argv$condition
comparename <- argv$comparename
design <- argv$design
fc <- argv$foldchange
p <- argv$pvalue
q <- argv$padj
genename <- argv$genename
outdir <- argv$outdir
setwd(outdir)

if ( !is.na(genename) ) {Genename <- read.delim(genename,header=T,row.names=1)}

compares <- strsplit(comparename,'vs')[[1]]
condition <- read.delim(condition, header=T)
if ( colnames(condition)[-1]=='type' ) { ngroup <- length(colnames(condition))-2}
if ( colnames(condition)[-1]!='type' ) { ngroup <- length(colnames(condition))-1}

if ( ngroup==1 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),]
	colnames(groupdata1)[2] <- 'groups'
	groupdata <- groupdata1}
if ( ngroup==2 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),-3]
	groupdata2 <- condition[which(condition$groups2 %in% compares),-2]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2)}
if ( ngroup==3 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),c(-3,-4)]
	groupdata2 <- condition[which(condition$groups2 %in% compares),c(-2,-4)]
	groupdata3 <- condition[which(condition$groups3 %in% compares),c(-2,-3)]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	colnames(groupdata3)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2,groupdata3)}
if (ngroup==4 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),c(-3,-4,-5)]
	groupdata2 <- condition[which(condition$groups2 %in% compares),c(-2,-4,-5)]
	groupdata3 <- condition[which(condition$groups3 %in% compares),c(-2,-3,-5)]
	groupdata4 <- condition[which(condition$groups4 %in% compares),c(-2,-3,-4)]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	colnames(groupdata3)[2] <- 'groups'
	colnames(groupdata4)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2,groupdata3,groupdata4)}
if ( ngroup==5 ) {
	groupdata1 <- condition[which(condition$groups1 %in% compares),c(-3,-4,-5,-6)]
	groupdata2 <- condition[which(condition$groups2 %in% compares),c(-2,-4,-5,-6)]
	groupdata3 <- condition[which(condition$groups3 %in% compares),c(-2,-3,-5,-6)]
	groupdata4 <- condition[which(condition$groups4 %in% compares),c(-2,-3,-4,-6)]
	groupdata5 <- condition[which(condition$groups5 %in% compares),c(-2,-3,-4,-5)]
	colnames(groupdata1)[2] <- 'groups'
	colnames(groupdata2)[2] <- 'groups'
	colnames(groupdata3)[2] <- 'groups'
	colnames(groupdata4)[2] <- 'groups'
	colnames(groupdata5)[2] <- 'groups'
	groupdata <- rbind(groupdata1,groupdata2,groupdata3,groupdata4,groupdata5)}

rownames(groupdata) <- groupdata[,1]
group <- factor(groupdata$groups)

if ( design=='normal' | design=='multi' ) {
	designs <- model.matrix(~group)}
if ( design=='batch' | design=='pair' | design=='blocking' ) {
	type <- factor(groupdata$type)
	designs <- model.matrix(~type+group)}


rawcount <- read.delim(rawcount, header=T, row.names=1)
counts <- subset(rawcount, select=rownames(groupdata))
counts <- na.omit(counts)
counts <- round(counts)
counts <- counts[rowSums(counts)>1,]

y <- DGEList(counts=counts, group=group)
y <- calcNormFactors(y)
rownames(designs) <- colnames(y)

if ( design=='normal' ) {

	if ( table(groupdata$groups)[[compares[1]]] + table(groupdata$groups)[[compares[2]]] > 2 ) {
		y <- estimateCommonDisp(y)
		y <- estimateTagwiseDisp(y)
		exact <- exactTest(y, pair=c(compares[2],compares[1]))
		result <- topTags(exact,n=NULL,sort.by='none')}

	if ( table(groupdata$groups)[[compares[1]]] + table(groupdata$groups)[[compares[2]]] == 2 ) {
		exact <- exactTest(y, pair=c(compares[2],compares[1]),dispersion=0.04)
		result <- topTags(exact,n=NULL,sort.by='none')
		y <- equalizeLibSizes(y)}}

if ( design=='batch' | design=='pair' | design=='multi' | design=='blocking' ) {
	y <- equalizeLibSizes(y)
	y <- estimateGLMCommonDisp(y,designs)
	y <- estimateGLMTrendedDisp(y,designs)
	y <- estimateGLMTagwiseDisp(y,designs)

	fit <- glmFit(y, designs)
	if ( design=='batch' | design=='pair' ) {lrt <- glmLRT(fit)}
	if ( design=='multi' ) {lrt <- glmLRT(fit,coef=2:ncol(designs))}
	if ( design=='blocking' ) {lrt <- glmLRT(fit,coef=(ncol(designs)-1):ncol(designs))}
	result <- topTags(lrt,n=NULL,sort.by='none') }

ID <- rownames(result)
if ( !is.na(genename) ) {Anno <- Genename[ID,]}

Count <- y$pseudo.counts
Count <- Count[ID,]

if ( is.na(genename) ) {
	result <- cbind(ID,as.data.frame(Count),as.data.frame(result))}
if ( !is.na(genename) ) {
	result <- cbind(ID,as.data.frame(Count),as.data.frame(result),Anno)}

if ( design=='normal' ) {
	result <- subset(result,select=-logCPM)}
if ( design!='normal' ) {
	result <- subset(result,select=-c(logCPM,LR))}

if ( design!='multi' & design !='blocking' ) {
	names(result)[names(result)=="logFC"] <- 'log2FoldChange'}

names(result)[names(result)=="PValue"] <- 'pvalue'
names(result)[names(result)=="FDR"] <- 'padj'
result$padj[is.na(result$padj)]  <- 1
result <- result[order(result$pvalue),]

if ( design!='multi' & design !='blocking' ) {
	if ( is.na(p) ) {
		ALL <- subset(result,padj <= q & abs(log2FoldChange) >= log(fc,2))
		UP <- subset(ALL,log2FoldChange > log(fc,2))
		DOWN <- subset(ALL,log2FoldChange < -log(fc,2))}
	if ( is.na(q) ) {
		ALL <- subset(result,pvalue <= p & abs(log2FoldChange) >= log(fc,2))
		UP <- subset(ALL,log2FoldChange > log(fc,2))
		DOWN <- subset(ALL,log2FoldChange < -log(fc,2))}}

if ( design=='multi' | design=='blocking' ) {
	if ( is.na(p) ) {
		ALL <- subset(result,padj <= q)}
	if ( is.na(q) ) {
		ALL <- subset(result,pvalue <= p)}}


write.table(result,file=paste(outdir,'/',comparename,'_DEG.xls',sep=''),sep='\t',quote=F,row.names=F)
write.table(ALL,file=paste(outdir,'/',comparename,'_DEG_ALL.xls',sep=''),sep='\t',quote=F,row.names=F)
write.table(UP,file=paste(outdir,'/',comparename,'_DEG_UP.xls',sep=''),sep='\t',quote=F,row.names=F)
write.table(DOWN,file=paste(outdir,'/',comparename,'_DEG_DOWN.xls',sep=''),sep='\t',quote=F,row.names=F)
