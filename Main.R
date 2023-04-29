#   Differential expression analysis with limma
#   Load the packages (make sure to install them via console before loading them)
library(GEOquery)
library(limma)
library(umap)

# This code loads a gene expression dataset from the Gene Expression Omnibus (GEO) 
#database using the getGEO() function from the GEOquery R package. 
#The dataset being loaded is identified by the accession code "GSE21962". 

# The if statement checks if the loaded dataset contains multiple platforms (i.e., microarray chips) or not.
# If there are multiple platforms, it selects the one with platform ID "GPL5175" using the grep() function 
# and assigns the corresponding index to idx.

gset <- getGEO("GSE21962", GSEMatrix =TRUE, AnnotGPL=FALSE)
idx <- 1
gset <- gset[[idx]]

# fvarLabels -> feature_variable_labels
print(fvarLabels(gset))

#In this line, the function fvarLabels() is used to retrieve and modify the feature variable labels (fvarLabels) of a GEOquery object gset.
#The make.names() function is used to ensure that the column names are syntactically valid in R, as column names in R have certain naming rules, such as not allowing spaces or special characters. This function replaces invalid characters with a dot (.) and adds a prefix if the first character is not a letter.
#Therefore, the purpose of the line fvarLabels(gset) <- make.names(fvarLabels(gset)) is to make sure that the fvarLabels of the gset object contain valid column names in R syntax.
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples (0 corresponds to case and 1 corresponds to normal)
gsms <- paste0("01010101010101010101010101010101010101010101010101",
               "010101010101010101010101010101010101010101")

#Split the string into individual characters
sml <- strsplit(gsms, split="")[[1]]

# Extract the expression matrix from the gset object
ex <- exprs(gset)
# If the vlaue are negative assign them as NaN
ex[which(ex <= 0)] <- NaN


library(impute)

# For the values which are missing, handle them using the KNN impute
# perform KNN imputation on the filtered expression matrix
imputed_exprs <- impute.knn(ex, k = 10)

#Applying the log2 transformation 
exprs(gset) <- log2(imputed_exprs[["data"]]) # log2 transform

# Normalise the data using the quantile based normalisation 
exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data


# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("case","control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)



print(n_sig) # The value of n_sig is 6397

# Once we got count of the DEG's, copy them to another table
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000)

# Since all the columns are not revelant to us we only extract the ones 
# which are concerned to us
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_LIST","SPOT_ID","RANGE_GB","RANGE_STRAND","RANGE_START"))

# Calculate how many genes are there that satisfy the conditions of p.value and fold change value range
n_sig <- sum(tT$P.Value < 0.01 & (tT$logFC < 1/1.5 | tT$logFC > 1.5))
print(n_sig)
tT <- head(tT,n=n_sig)

# We will also write this table to a csv file to do further analysis
write.csv(tT, file = "DEGs.csv", row.names = TRUE)



# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))




