
##just need to change working directory, file names, and file paths and run in directory with files in it. 
##input is SAM format resulting from Bismark alignment (with no header)

library(dplyr)
library(stringr)
library(readr)


##set unique file ID's
unique_fileID <- c("SRR13012805",
                   "SRR13012806",
                   "SRR13012807",
                   "SRR13012808",
                   "SRR13012809",
                   "SRR13012810",
                   "SRR13012811",
                   "SRR13012812",
                   "SRR13012813",
                   "SRR13012814",
                   "SRR13012815",
                   "SRR13012816",
                   "SRR13012817",
                   "SRR13012818",
                   "SRR13012819",
                   "SRR13012820",
                   "SRR13012821",
                   "SRR13012822",
                   "SRR13012823",
                   "SRR13012824",
                   "SRR13012825",
                   "SRR13012826",
                   "SRR13012827",
                   "SRR13012828",
                   "SRR13012829",
                   "SRR13012830",
                   "SRR13012831",
                   "SRR13012832",
                   "SRR13012833",
                   "SRR13012834",
                   "SRR13012835",
                   "SRR13012836",
                   "SRR13012837",
                   "SRR13012838",
                   "SRR13012839",
                   "SRR13012840",
                   "SRR13012841",
                   "SRR13012842",
                   "SRR13012843",
                   "SRR13012844",
                   "SRR13012845",
                   "SRR13012846",
                   "SRR13012847",
                   "SRR13012848",
                   "SRR13012849",
                   "SRR13012850",
                   "SRR13012851",
                   "SRR13012852",
                   "SRR13012853",
                   "SRR13012854",
                   "SRR13012855",
                   "SRR13012856",
                   "SRR13012857",
                   "SRR13012858",
                   "SRR13012859",
                   "SRR13012860",
                   "SRR13012861",
                   "SRR13012862",
                   "SRR13012863",
                   "SRR13012864",
                   "SRR13012865",
                   "SRR13012866",
                   "SRR13012867",
                   "SRR13012868",
                   "SRR13012869",
                   "SRR13012870",
                   "SRR13012871",
                   "SRR13012872",
                   "SRR13012873",
                   "SRR13012874",
                   "SRR13012875",
                   "SRR13012876",
                   "SRR13012877",
                   "SRR13012878",
                   "SRR13012879",
                   "SRR13012880",
                   "SRR13012881",
                   "SRR13012882",
                   "SRR13012883",
                   "SRR13012884",
                   "SRR13012885",
                   "SRR13012886",
                   "SRR13012887",
                   "SRR13012888",
                   "SRR13012889",
                   "SRR13012890",
                   "SRR13012891",
                   "SRR13012892",
                   "SRR13012893",
                   "SRR13012894",
                   "SRR13012895",
                   "SRR13012896",
                   "SRR13012897",
                   "SRR13012898",
                   "SRR13012899",
                   "SRR13012900",
                   "SRR13012901",
                   "SRR13012902",
                   "SRR13012903",
                   "SRR13012904",
                   "SRR13012905",
                   "SRR13012906",
                   "SRR13012907",
                   "SRR13012908",
                   "SRR13012909",
                   "SRR13012910",
                   "SRR13012911",
                   "SRR13012912",
                   "SRR13012913",
                   "SRR13012914",
                   "SRR13012915",
                   "SRR13012916",
                   "SRR13012917",
                   "SRR13012918",
                   "SRR13012919",
                   "SRR13012920",
                   "SRR13012921",
                   "SRR13012922",
                   "SRR13012923",
                   "SRR13012924",
                   "SRR13012925",
                   "SRR13012926",
                   "SRR13012927",
                   "SRR13012928",
                   "SRR13012929",
                   "SRR13012930",
                   "SRR13012931",
                   "SRR13012932",
                   "SRR13012933",
                   "SRR13012934",
                   "SRR13012935",
                   "SRR13012936",
                   "SRR13012937",
                   "SRR13012938")

##set file paths
filePaths <- paste0("/scratch/emb19132/GENOME_PDR/RAT/sam_no_header/", unique_fileID, "_sorted.txt")

##make for loop to calculate PDN
for(i in 1:length(unique_fileID)){
  ##read in SAM file written as a plain text file with headers removed
  SAM <- read.table(filePaths[i], header = FALSE, fill = TRUE)
  ##remove XM:Z: string tag from methylation string (column V14)
  SAM$plain.string <- str_remove_all(SAM$V14, "XM:Z:")
  ##count number of CpGs in string
  SAM$number.of.unmeth.cpg <- str_count(SAM$plain.string, "z")
  SAM$number.of.meth.cpg <- str_count(SAM$plain.string, "Z")
  SAM$number.of.cpg <- (SAM$number.of.meth.cpg + SAM$number.of.unmeth.cpg)
  SAM$actual.read.length <- str_count(SAM$plain.string, "") ##should be same as CIGAR specification if no mismatches
  ##extract just the CpGs in the methylation string
  SAM$cpg.meth.string <- gsub("[^zZ]", "", SAM$plain.string)
  ##calculate number of neighbors (i.e. observations of pairs that could be concordant/discordant)
  SAM$num.neighbors <- (((SAM$number.of.cpg - 2)*2)+(2))
  ##set negative values (to cytosines with no cytosines) to zero
  SAM$num.neighbors[SAM$num.neighbors < 0] <- 0
  ##count discordant occcurances of a string
  SAM$zZ.occur <- str_count(SAM$cpg.meth.string,"zZ")
  SAM$Zz.occur <- str_count(SAM$cpg.meth.string,"Zz")
  ##for each occurance of zZ or Zz, read loses 2 concordance points from max concordance
  SAM$discordant.occur <- (SAM$zZ.occur + SAM$Zz.occur)
  SAM$discordance.deduction <- (SAM$discordant.occur*2)
  SAM$concordance.score <- (SAM$num.neighbors - SAM$discordance.deduction)
  SAM$concordance.perc <- (SAM$concordance.score/SAM$num.neighbors)
  ##print only read information as a text file for each sample
  SAM$zero.based.start <- (SAM$V4 - 1)
  SAM$zero.based.end <- (SAM$zero.based.start + SAM$actual.read.length)
  SAM$POS <- SAM$V4
  SAM$CIGAR <- SAM$V6
  SAM$CHROM <- SAM$V3
  cleaned.SAM <- na.omit(SAM) ##remove NAs
  write.table(cleaned.SAM, file=paste0("/scratch/emb19132/GENOME_PDR/RAT/PDN_output/", unique_fileID[i], "_cleaned_SAM.txt"))
  ##now make txt file in bed format for each with the 0 based start and end positions after removing cytosines without a neighbor on the nearest read (0 or 1 CpG)
  perRead_bed <- cleaned.SAM[,c(34,30,31,1,29)]
  #Write the result table in an individual text file
  write.table(perRead_bed, file=paste0("/scratch/emb19132/GENOME_PDR/RAT/PDN_output/", unique_fileID[i], "_PDN_per_read.txt"))
}



