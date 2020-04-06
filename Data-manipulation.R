library(tidyverse)


#import dataset
dab <- read_delim("VDJBdatabase.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


#combining TRA and TRB in a single dataframe with controles 
for (i in dab$complex.id){ #for loop over all ids
  if (i %in% done != TRUE){ #if id not in list with previously done ids, continue
    done <- c(done, i) #add id to list with previously done ids
    pos <- (which(dab$complex.id == i)) #find all positions of the id
    for (j in pos){ #for loop over all positions of id
      r <- dab[j,]  #assign row to new dataframa
      if (r$Gene == "TRA"){ #filter on TRA
        zonder <- str_remove_all(r$Meta,"\"") #removal off "\" in Meta
        stuks <- strsplit(zonder,",") #split meta in categories
        for (z in stuks[[1]]){ #loop over categories
          if (isTRUE(grepl('subject.id:', z)) == TRUE){ #if categorie contains subject.id, continue
            #rename data for distinquishing between TRA and TRB
            CDR3ai <- r$CDR3
            TRAVi <- r$V
            TRAJi <- r$J
            if(isTRUE(is.null(PCi)) == TRUE){ #if patientcount (PCi) is empty
              PCi <- str_remove(z," subject.id: ") #assign patientcount to PCi
            }else if(PCi != str_remove(z," subject.id: ")){ #else if patientcount (PCi) is not containing the same patient id 
              print(c("ERROR",i)) #give an error
            }
          }
        }
      } else if(r$Gene == "TRB"){ #filter on TRB
        zonder <- str_remove_all(r$Meta,"\"") #removal off "\" in Meta
        stuks <- strsplit(zonder,",") #split meta in categories
        for (z in stuks[[1]]){ #loop over categories
          if (isTRUE(grepl('subject.id:', z)) == TRUE){ #if categorie contains subject.id, continue
            #rename data for distinquishing between TRA and TRB
            CDR3bi <- r$CDR3
            TRBVi <- r$V
            TRBJi <- r$J
            if(isTRUE(is.null(PCi)) == TRUE){ #if patientcount (PCi) is empty
              PCi <- str_remove(z," subject.id: ") #assign patientcount to PCi
            }else if(PCi != str_remove(z," subject.id: ")){ #else if patientcount (PCi) is not containing the same patient id 
              print(c("ERROR",i)) #give an error
            }
          }
        }
      }
    }
    #add previously defined data to new dataframe and empty PCi
    nie <- add_row(nie, complex.id = i,CDR3b = CDR3bi, TRBV = TRBVi, TRBJ = TRBJi, CDR3a = CDR3ai, TRAV = TRAVi, TRAJ = TRAJi, PatientCounts = PCi)
    PCi <- NULL
  }
}
#write data to DATAZ.txt without column complex.id
DATA <- nie
DATA$complex.id <- NULL
write.table(DATA, "D:/gliph-1.0/gliph/bin/DATAz.txt", sep="\t", row.names = FALSE, quote = FALSE)


#make a new dataframe with a new column containing patient ids
for(e in dab$Meta){ #loop over meta
  verw <- str_remove_all(e,"\"")#remove "\"
  gesplit <- strsplit(verw,",")#split in categories
  for (x in gesplit[[1]]){ #loop over categories
    if (isTRUE(grepl('subject.id:', x)) == TRUE){#if categorie contains subjectid
      subid <- str_remove(x," subject.id: ") #name subject id as subid
      niecol <-c(niecol, subid) #add subid to new column
    }
  }
}
#add subids to new dataframe containing the all data
dabsubid <- dab
dabsubid$subid <- niecol


#make a table with HLA information for GLIPHs HLA association test
for(e in dabsubid$subid){ #for every subid
  if (e %in% klaar != TRUE){#check if subid is in list with already done subids
    klaar <- c(klaar, e)#add subid to list with already done subids
    pla <- (which(dabsubid$subid == e)) #on which positions is this subid located
    nierow <- c(e) #add this subid to the beginning of a new list (nierow)
    for(p in pla){ #for every position of the subid
      if(dabsubid[p,]$`MHC A` %in% nierow != TRUE){ #if MHC not yet in nierow
        nierow <- c(nierow, dabsubid[p,]$`MHC A`) #add to nierow
      }
      if(dabsubid[p,]$`MHC B` %in% nierow != TRUE){ #see above, only for MHC B instead of MHC A
        nierow <- c(nierow, dabsubid[p,]$`MHC B`)
      }
    }
    #write nierow and a linebreak to file containing the HLA tabel and empty nierow
    write_lines(nierow, "D:/gliph-1.0/gliph/bin/HLA_tabel.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/HLA_tabel.txt" ,append = TRUE, sep = "")
    nierow <- NULL
  }
}


#Make an HLA_data txt file for the visualisation
for(q in dabsubid$`MHC A`){ #for loop over every MHC A
  if (q %in% qlaar != TRUE){#check if MHC is in list with already done MHCs 
    qlaar <- c(qlaar, q)#add MHC to list with already done MHCs
    qpla <- (which(dabsubid$`MHC A` == q)) #get positions of every occurunce of this MHC
    nieqrow <- c(q) #make a new list with the MHC as first entry
    for(qp in qpla){ #loop over all positions
      if(dabsubid[qp,]$CDR3 %in% nieqrow != TRUE){ #if CDR3 on this location not yet in list 
        nieqrow <- c(nieqrow, dabsubid[qp,]$CDR3) #add CDR3 to list
      }
    }
    #write nieqrow and a linebreak to file containing the HLA tabel and empty qnierow
    write_lines(nieqrow, "D:/gliph-1.0/gliph/bin/HLA_DATA.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/HLA_DATA.txt" ,append = TRUE, sep = "")
    nieqrow <- NULL
  }
}


#getting the individual cluster can be achieved by changing the variable nr, in the line below and the name of the file in line 126
nr <- 2534 #change this number to the cluster you want the information of
afdr <- NULL #empty afdr
clus <- DATAz.convergence.groups #import the clusters
sclus <- strsplit(as.character(clus[nr,3])," ")#split up the target sequences in individual sequences
for (s in sclus){
  afdr <- c(afdr, s) #add every sequence to afdr
}
afdr <- c(nr,str_remove(clus[nr,2],"CRG-"),afdr) #add the number of the cluster, the source sequences and all target sequences of the cluster
write_lines(afdr, "D:/gliph-1.0/gliph/bin/cluster2534.txt" ,append = TRUE, sep="\t") #write this line to a text file of the cluster to be used in visualisation
#list of clusters i have used in visualisation: 18,2,139,111,33,113,591,457,1370,2534

