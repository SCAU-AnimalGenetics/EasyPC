############################################################################################
#                       Easy Pedigree Checking (EasyPC) V1.0                               #
#                                                                                          #                                                                      
#  A software for pedigree checking with high density SNP data.                            #                                                                      
#                                           Designed by Zhe Zhang, Yuanyu Luo & Jinlong He #
#                                                           Email: zhezhang@scau.edu.cn    #
############################################################################################

# Input file: a genotype file and a pedigree file
#   The first column in the genotype file is ID that corressponding to the ID in pedigree file.
#   The rest columns in the genotype file are SNP genotypes coded as0, 1, 2, and 3 for 
#             unknown, aa, Aa, and AA, respectively.   
#   The pedigree file includes two columns with the offsprings in the first and a maternal or
#             paternal candidate in the second column. 
# Output file??
#   'Mendelian_error_rate.txt' shows the results for mendelian error checking between all individual pairs. 
#        Four colcumns in this file include two IDs for the first and second invidual, mendelian 
#              error count and mendlian error rate.  
#  'ped_right.txt' shows the correct pedigree appeared in 'pedigree file' that identifyed by EasyPC.
#  'ped_error.txt' shows the wrong pedigree appeared in 'pedigree file' that identifyed by EasyPC.
#  'ped_no_gene.txt' shows the pedigree that include ungenotyped individual(s). 
#  'ped_right_inferred.txt' shows the inferred correct invidual paris that did not appear in 'pedigree file'.

##============================ functions starts here ==================
mend.error <- function(gen){
  dir.create("EasyPCres", showWarnings=F); setwd("./EasyPCres")
  mend_err <- matrix(0, nrow(gen)**2, 4)
  colnames(mend_err) <- c("off","par","Mendelian_error","error_rate")
  mend_err[,1] <- rep(gen[,1], each=nrow(gen))
  mend_err[,2] <- rep(gen[,1], times=nrow(gen))
  row.names(gen[,1]); gen <- gen[,-1]
  for(i in 1:nrow(gen)){ # off
    for(j in 1:nrow(gen)){ # par
      err <- sum((gen[i, ] == 3 & gen[j, ] ==1) | (gen[i, ] == 1 & gen[j,] ==3))
      if(sum(gen[j, ] == 2) != ncol(gen)){
        errate <- err / sum((gen[i, ]*gen[j, ] != 0) & (gen[j, ] !=2)) 
      }else{ errate <- 999}
      mend_err[(i-1)*nrow(gen)+j, 3] <- err
      mend_err[(i-1)*nrow(gen)+j, 4] <- errate
    }
  }
  mend_err <- mend_err[mend_err[,1] != mend_err[,2], ]
  errate_plot <- hist(as.numeric(mend_err[,4]),main="",xlab="Mendelian error rate",ylab="Frequency")
  write.table(mend_err,"Mendelian_error_rate.txt", row.names=F, col.names=T, quote=F)
  setwd("../")
  return(mend_err=mend_err)
}
ped.infer <- function(ped, mend_err, thre=0.01){
  dir.create("EasyPCres", showWarnings=F); setwd("./EasyPCres")
  mend_right <- mend_err[as.numeric(mend_err[,4])<thre,]
  mend_error <- mend_err[as.numeric(mend_err[,4])>=thre,]
  tmp1 <- paste(ped[,1], ped[,2], sep="")
  tmp2 <- paste(mend_err[,1], mend_err[,2], sep="")
  tmp3 <- paste(mend_right[,1], mend_right[,2], sep="")
  tmp4 <- paste(mend_error[,1], mend_error[,2], sep="")
  no.gene <- ped[!tmp1 %in% tmp2,]
  write.table(no.gene, "ped_no_gene.txt", row.names=F, col.names=T, quote=F)
  ped <- ped[tmp1 %in% tmp2,]
  right <- tmp3 %in% tmp1
  error <- tmp4 %in% tmp1
  Num.of.ped.no.gene <- nrow(no.gene)
  Num.of.ped.right <- sum(right)
  Num.of.ped.error <- nrow(ped) - sum(right)
  Errate.of.ped <-   Num.of.ped.error/(Num.of.ped.right+Num.of.ped.error)*100
  ped.right <- mend_right[right,]
  ped.error <- mend_error[error,]
  ped.right.infer <- mend_right[!tmp3 %in% paste(ped[,1], ped[,2], sep=""),]
  ped.right.infer <- ped.right.infer[!paste(ped.right.infer[,2],ped.right.infer[,1],sep="") %in% paste(ped[,1], ped[,2], sep=""),]
  Num.of.ped.right.infer <- nrow(ped.right.infer)
  res <- matrix(0,5,2)
  res[,1] <- c("Ungenotyped_Ped_Num","Ped_Correct_Num", "Ped_Error_Num","Inferred_Ped_Num","Ped_Error_Rate")
  res[,2] <- c(Num.of.ped.no.gene,Num.of.ped.right, Num.of.ped.error,Num.of.ped.right.infer,Errate.of.ped)
  write.table(ped.right, "ped_right.txt", quote=F, row.names=F)
  write.table(ped.error,"ped_error.txt", quote=F, row.names=F)
  write.table(ped.right.infer,"ped_right_infer.txt", quote=F, row.names=F)
  write.table(res,"result.txt", quote=F, row.names=F, col.names=F)
  setwd("../")
}

##============================ program starts here ==================
setwd("c:/")			# set your work directory
gen <- as.matrix(read.table("gen.txt", header=F)) 	# read genotype file
ped <- as.matrix(read.table("ped.txt", header=T)) 	# read pedigree fil
mend_err <- mend.error(gen)				# calculte mendelian error rate within population
ped.infer(ped, mend_err, thre=0.01)			# pedigree cheacking, thre= threshold value.







