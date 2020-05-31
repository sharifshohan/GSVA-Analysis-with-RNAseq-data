# The following code gives you the main patheays that are involved in your sample types based on gene set variation analysis

#At first you need to set the directory
#run the code code in where your "mytable.csv" file was saved
##setwd("/Users/mohammadumersharifshohan/Desktop/Sample projects online/")

mydata <- read.csv(file= "T_ALL_rlog2 copyBIS.csv",header=TRUE) #Load your data
mydata #Check your data
length(mydata$X) #check the length of the mydata file. basically how many rows are there
#As the GSV package works on based on gene name and sample. But the gene names needs to be in Entrez ID format. So the next few lines we will 
#collect the name and split the name sothat we can get the basic ENSEMBL ID data. the data given looks to have ENSG00000268663.1_4 and we are 
#going to get rid of the ".1_4". So anything after "." we will get rid of to run using BiomaRt and get ht Entrez ID

rowname <- mydata$X #Store the ensembl gene name into rowname
View(rowname) #Check if you have got the names

#mydata$X <- NULL
x<- strsplit(rowname, "\\.") #Slit the name based on where "." is and this code divides into two parts the name. 
x# check x
z<- unlist(x) #We need to unlist the data

y<- matrix(z, ncol=2, byrow=TRUE) # This code makes it a matrix with two columns. 1st coulmn before  "." and 2nd column after "."
y
#The follong two codes are to check whether we have missed any data till now
length(mydata$X)
length(y[,1]) #This line means all the data in the first column and then we check the length
y[,1] #You can see what is in y[,1]. This is everything in the first row of the matrix
new_rowname <- y[,1] #Save the name in a new variable called new_rowname
length(new_rowname) #check whether we have the same legth as before

#As we now have the ID of ENSEMBL in proper format. We will now convert it to Entrez ID and But before we will do some checks
#whether our converted id and the ensembl id in the main file is in the proper position. IF not then we will have problem 

#Now, we are just going to give the mydata file id for each gene. Starting from one to the number of rows present.
mydata$id  <- 1:nrow(mydata)
mydata #you can see a new id in the end of the file in your console
# In the following code we have the new_rownames that we created that does not have "." in it and added it to main mydata file as new column

mydata$X1 <- new_rowname ## the reason we have given this column name "X1" is because it will help us merge the names with entrez id. Remember the line name
mydata #Check if it was added
length(mydata$X1) # check if the length is still the same or whether we have missed any file or not. Looks good 62299

#------
##BiocManager::install("biomaRt") install this package if not installed
library(biomaRt) #load the package

#The following code is inbuilt code for biomaRt to convert ensembl id to entrez id. Although you will see that all the ids could not be converted
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "http://www.ensembl.org")
#Sometimes shows that they will try some other hosting service. Ignore the error
#the following code will check the ensembl gene id from your file "||value = new_rowname||" and give you a new file with two columns that have
#ensembl_gene_id and entrezgene_id
genes <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","entrezgene_id"),
               values = new_rowname, 
               mart = mart)
genes #Check what is inside genes. You can see there are NA in the file. That means it could not find the entrez id for that ensembl 
colnames(genes) <- c("X1","V1") #we Are changing the column name to for merging as stated in line 37
genes

#Now we the want to merge the genes$X1 column to mydata$X1 column and get the gene$V1 or the entrez id in the mydata file
#The following code is for ordering purpose. This is a functionc created for ordering
keeping.order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
}
genes
#In the following code we created a mydata2 a new file. here we merge genes file and mydata file by column X1. So we try to match them and
#create a new mydata2 file with all the contents from mydata present
mydata2<- keeping.order(mydata, merge, y=genes, by = "X1", all.x = TRUE) 
mydata2
##you can see mydata2 file with lots of NA in the V1 which is the entrez ID. You can also see the file by writing it out in your desktope
###write.csv(mydata2, "newfile.csv") //although not necessary
##The following codes we are getting rid of the un necessary columns in the mydata2 file
mydata2$X1 <- NULL
mydata2$X <- NULL
mydata2$id <- NULL 
mydata2
library(tidyr)
mydata2<- mydata2 %>% drop_na() ##We dropped the files with no Entrez IDs and got 20,000 genes with entrez id to work with
length(mydata2$ALLSIL) #Checked the length 19452 genes
mydata2

mydata2 <- mydata2[!duplicated(mydata2$V1), ]  ## Removing all the duplicates from the code

rownames(mydata2) <- mydata2$V1 #here we are giving the rownames the entrez id name. so column v1 is added to the rownames
mydata2
#As you can see there were few duplicates. the code got rid of them
mydata2$V1 <- NULL # we get rid of the V1 column now
mydata2
length(mydata2$ALLSIL)  #check the length of mydata2 that is created for final
length(mydata$ALLSIL) # The main file 

########## GSVA ############
#Now we are removing rows with no values or 0
keep <- apply(mydata2, 1, function(x) var(x, na.rm = TRUE) > 0)
keep
mydata2 <- mydata2[keep,]
length(mydata2$ALLSIL) ##Around 200 more genes were removed
mydata2 <- data.matrix(mydata2) #Convert it to data matrix

#we need a set of genes for comparison what our file has genes expression more from
library(GSEABase)
library(GSVAdata)
data(c2BroadSets)

canonicalC2BroadSets <-c2BroadSets[c(grep("^KEGG",names(c2BroadSets)),
                                     grep("^REACTOME", names(c2BroadSets)),
                                     grep("^BIOCARTA", names(c2BroadSets)))]

esrnaseq <- gsva(mydata2, canonicalC2BroadSets, min.sz=5, max.sz=500,
                 kcdf="Poisson", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

write.csv(esrnaseq, "esrnaseq_new.csv")

