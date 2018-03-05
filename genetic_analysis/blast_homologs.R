## CS640 Bioinformatics - Fall 2017
## Author: Melanie Baybay
## Last Modified: Sept. 21, 2017


# Adapted from page bottom of:
# http://www.molgen.ua.ac.be/bioinfo/acourse/Bioconductor-by-examples.txt

# see also:
#http://bioconductor.fhcrc.org/packages/release/bioc/vignettes/Biostrings/inst/doc/MultipleAlignments.R

# Some web services based packages
# ================================

#some useful packages
install.packages("RCurl")
library(RCurl)
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biostrings")
library(Biostrings)
install.packages("rentrez")
library(rentrez)

setwd("/Users/MelanieBaybay/Documents/USF/2017-2018Fall/Bioinformatics/R-Assignments")


#Here is the FASTA sequence provided in the BLAST assignment:
seq <- "HLYPGEVCPGMDIRNNLTRLHELENCSVIEGHLQILLMFKTRPEDFRDLSFPKLIMITDYLLLFRVYGLESLKDLFPNLTVIRGSRLFFNYALVIFEMVHLKELGLYNLMNITRGSVRIEKNNELCYLATIDWSRILDSVEDNHIVLNKDDNEECGDICPGTAKGKTNCPATVINGQFVERCWTHSHCQKVCPTICKSHGCTAEGLCCHSECLGNCSQPDDPTKCVACRNFYLDGRCVETCPPPYYHFQDWRCVNFSFCQDLHHKCKNSRRQGCHQYVIHNNKCIPECPSGYTMNSSNLLCTPCLGPCPKVCHLLEGEKTIDSVTSAQELRGCTVINGSLIINIRGGNNLAAELEANLGLIEEISGYLKIRRSYALVSLSFFRKLRLIRGETLEIGNYSFYALDNQNLRQLWDWSKHNLTITQGKLFFHYNPKLCLSEIHKMEEVSGTKGRQERNDIALKTNGDKASCENELLKFSYIRTSFDKIS"

#After finding the best hit that is actual, not predicted sequence, we can use its accession number
#seq <- "NP_000199.2"

# blast a sequence using the NCBI REST interface
# enable low-complexity filtering with FILTER="L", whereas FILTER="F" turns filltering off
# FORMAT_TYPE="Text" vs "HTML"
# use NCBI GI numbers in the output page with NCBI_GI="on"  
# May no longer be supported (not listed in current API):
#    CLIENT="web", SERVICE="plain"
ncbi_blast = 
  function(seq,database="nr",expect=10,program="blastp",hitlistsize=100, filter="F", matrix="BLOSUM90", word_size=6, n_alignments=10, format_type="HTML") {
  # use of Rcurl,getForm(): employs NCBI REST interface (https://ncbi.github.io/blast-cloud/dev/api.html)
  job = as.character(getForm("https://blast.ncbi.nlm.nih.gov/Blast.cgi", QUERY=seq, DATABASE=database, MATRIX=matrix, 
                             HITLIST_SIZE=hitlistsize, ALIGNMENTS=n_alignments, FILTER=filter, EXPECT=expect, FORMAT_TYPE=format_type, PROGRAM=program, CLIENT="web",
                             WORD_SIZE=word_size, SERVICE="plain", NCBI_GI="on", CMD="Put"))
  m = regexpr("RID = ([^\n]+)",job)
  id = substring(job,m[1]+6,m[1]+attributes(m)$match.length-1)
  m = regexpr("RTOE = ([^\n]+)",job)
  rtoe = substring(job,m[1]+7,m[1]+attributes(m)$match.length-1)
  result = "Status=WAITING"
  while (length(grep("Status=WAITING",result))) {
    #Sys.sleep(1)
    Sys.sleep(60)
    result = as.character(getForm("https://blast.ncbi.nlm.nih.gov/Blast.cgi", RID=id, FORMAT_TYPE="Text", CMD="Get"))
  }
  getForm("https://blast.ncbi.nlm.nih.gov/Blast.cgi", RID=id, CMD="Delete")
  return(result)
}
b = ncbi_blast(seq,database="refseq_protein",program="blastp") 

write.table(b, "blast_data0_humanProteinSeq.html", sep="\t")

#1
# Starting with the accession number of the protein that was your best hit of actual
# (not predicted) sequences in the above search
# perform a new search for distantly related homologs
# Save in a file named "blast_data#1_distantHomologs.html"
seq <- "NP_000199.2"
b1 <- ncbi_blast(seq, database="refseq_protein", matrix="PAM250")
write.table(b1, "blast_data#1_distantHomologs.html", sep="\t")

#2
# Use a BLAST search to find the accession number of the nucleotide sequence that corresponds with the
# protein sequence that you used as query above. Retrieve the pairwise alignments
# of the best candidates as aligned with your query. 
# Save in file named "blast_data#2_nucleotideSeq.html
b2 <- ncbi_blast(seq, database="nr", program="tblastn", hitlistsize=20, n_alignments=5)
write.table(b2, "blast_data#2_nucleotideSeq.html", sep="\t")

#3
#Find homologs to the SLC16A13 protein (paralog of the query protein) 
# that you used in part 1  #5
# Save in a file named "blast_data#3_SLC16A13Homologs.html"
seq <- "NP_963860.1" # accession number of SLC16A13
b3 <- ncbi_blast(seq, database="refseq_protein", expect=30, program="blastp", matrix="PAM250")
write.table(b3, "blast_data#3_SLC16A13Homologs.html", sep="\t")

