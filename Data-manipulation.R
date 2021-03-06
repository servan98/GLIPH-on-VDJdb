library(tidyverse)
library(stringr)

#import dataset
dab <- read_delim("ABsam.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


#combining TRA and TRB in a single dataframe with controles
done <- c()
PCi <- NULL
nie <- tibble(complex.id = NA,CDR3b = NA, TRBV = NA, TRBJ = NA, CDR3a = NA, TRAV = NA, TRAJ = NA, PatientCounts = NA)

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
niecol <- c()

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
klaar <-c()

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
qlaar <- c()
for(q in dabsubid$`MHC A`){ #for loop over every MHC A
  if (q %in% qlaar != TRUE){#check if MHC is in list with already done MHCs 
    qlaar <- c(qlaar, q)#add MHC to list with already done MHCs
    qpla <- (which(dabsubid$`MHC A` == q)) #get positions of every occurunce of this MHC
    sp <- strsplit(q,":")
    if(length(sp[[1]])>1){
      typ <- str_c(sp[[1]][[1]],':',sp[[1]][[2]])
    }
    else{
      typ <- sp[[1]][[1]]
    }
    nieqrow <- c(typ) #make a new list with the MHC as first entry
    for(qp in qpla){ #loop over all positions
      if(dabsubid[qp,]$CDR3 %in% nieqrow != TRUE){ #if CDR3 on this location not yet in list 
        nieqrow <- c(nieqrow, dabsubid[qp,]$CDR3) #add CDR3 to list
      }
    }
    #write nieqrow and a linebreak to file containing the HLA tabel and empty nieqrow
    write_lines(nieqrow, "D:/gliph-1.0/gliph/bin/HLA_DATAres3.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/HLA_DATAres3.txt" ,append = TRUE, sep = "")
    nieqrow <- NULL
  }
}


#Make an Patientid_data txt file for the visualisation
slaar <- c()
for(s in dabsubid$subid){ #for loop over every patient id
  if (s %in% slaar != TRUE){#check if pID is in list with already done pIDs 
    slaar <- c(slaar, s)#add pID to list with already done pIDs
    spla <- (which(dabsubid$subid == s)) #get positions of every occurunce of this pID
    niesrow <- c(s) #make a new list with the pID as first entry
    for(sp in spla){ #loop over all positions
      if(dabsubid[sp,]$CDR3 %in% niesrow != TRUE){ #if CDR3 on this location not yet in list 
        niesrow <- c(niesrow, dabsubid[sp,]$CDR3) #add CDR3 to list
      }
    }
    #write niesrow and a linebreak to file containing the Patienid tabel and empty niesrow
    write_lines(niesrow, "D:/gliph-1.0/gliph/bin/Patientid_DATA.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/Patientid_DATA.txt" ,append = TRUE, sep = "")
    niesrow <- NULL
  }
}

#Make an VgeneUsage_data txt file for the visualisation
vlaar <- c()
for(v in dabsubid$V){ #for loop over every V gene
  if (v %in% vlaar != TRUE){#check if V gene is in list with already done V genes 
    vlaar <- c(vlaar, v)#add V gene to list with already done V genes
    vpla <- (which(dabsubid$V == v)) #get positions of every occurunce of this V gene
    nievrow <- c(v) #make a new list with the V gene as first entry
    for(vp in vpla){ #loop over all positions
      if(dabsubid[vp,]$CDR3 %in% nievrow != TRUE){ #if CDR3 on this location not yet in list 
        nievrow <- c(nievrow, dabsubid[vp,]$CDR3) #add CDR3 to list
      }
    }
    #write nievrow and a linebreak to file containing the Vgeneusage tabel and empty nievrow
    write_lines(nievrow, "D:/gliph-1.0/gliph/bin/VgeneUsage_DATA.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/VgeneUsage_DATA.txt" ,append = TRUE, sep = "")
    nievrow <- NULL
  }
}

#Make an JgeneUsage_data txt file for the visualisation
jlaar <- c()
for(j in dabsubid$J){ #for loop over every J gene
  if (j %in% jlaar != TRUE){#check if J gene is in list with already done J genes 
    jlaar <- c(jlaar, j)#add J gene to list with already done J genes
    jpla <- (which(dabsubid$J == j)) #get positions of every occurunce of this V gene
    niejrow <- c(j) #make a new list with the J gene as first entry
    for(jp in jpla){ #loop over all positions
      if(dabsubid[jp,]$CDR3 %in% niejrow != TRUE){ #if CDR3 on this location not yet in list 
        niejrow <- c(niejrow, dabsubid[jp,]$CDR3) #add CDR3 to list
      }
    }
    #write niejrow and a linebreak to file containing the Jgeneusage tabel and empty niejrow
    write_lines(niejrow, "D:/gliph-1.0/gliph/bin/JgeneUsage_DATA.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/JgeneUsage_DATA.txt" ,append = TRUE, sep = "")
    niejrow <- NULL
  }
}

#Make an Epitope_data txt file for the visualisation
elaar <- c()
for(e in dabsubid$Epitope){ #for loop over every Epitope
  if (e %in% elaar != TRUE){#check if Epitope is in list with already done Epitopes 
    elaar <- c(elaar, e)#add Epitope to list with already done Epitopes
    epla <- (which(dabsubid$Epitope == e)) #get positions of every occurunce of this Epitope
    nieerow <- c(e) #make a new list with the Epitope as first entry
    for(ep in epla){ #loop over all positions
      if(dabsubid[ep,]$CDR3 %in% nieerow != TRUE){ #if CDR3 on this location not yet in list 
        nieerow <- c(nieerow, dabsubid[ep,]$CDR3) #add CDR3 to list
      }
    }
    #write nieerow and a linebreak to file containing the Eptipe tabel and empty niejrow
    write_lines(nieerow, "D:/gliph-1.0/gliph/bin/Epitope_DATA.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/Epitope_DATA.txt" ,append = TRUE, sep = "")
    nieerow <- NULL
  }
}


#getting the individual cluster can be achieved by changing the variable nr in the line below, the " ".convergence.groups and the name of the file in line 220
nr <- 2534 #change this number to the cluster you want the information of
afdr <- NULL #empty afdr
clus <- read.delim("D:/gliph-1.0/gliph/bin/DATAz-convergence-groups.txt", header=FALSE) #import the clusters
sclus <- strsplit(as.character(clus[nr,3])," ")#split up the target sequences in individual sequences
for (s in sclus){
  afdr <- c(afdr, s) #add every sequence to afdr
}
afdr <- c(nr,str_remove(clus[nr,2],"CRG-"),afdr) #add the number of the cluster, the source sequences and all target sequences of the cluster
write_lines(afdr, "D:/gliph-1.0/gliph/bin/cluster2534.txt" ,append = TRUE, sep="\t") #write this line to a text file of the cluster to be used in visualisation
#list of clusters i have used in visualisation: 18,2,139,111,33,113,591,457,1370,2534
