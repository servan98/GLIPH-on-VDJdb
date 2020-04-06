#dab <- read_delim("ABsam.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
library(tidyverse)
ptm <- proc.time()
for (i in dab$complex.id){# for loop over alle id's
  if (i < 100000 ){#verklein database tot eerste 1000 complexen
    if (i %in% done != TRUE){#kijk of id in de lijst met gedane ids staat
      done <- c(done, i)#id aan lijst met gedane ids toevoegen
      pos <- (which(dab$complex.id == i))#op welke posities staat de id
      for (j in pos){ #for loop over alle posities || straks hier een loop over complete dab van maken
        r <- dab[j,]  #assign de hele regel naar n nieuwe dataframe
        if (r$Gene == "TRA"){ #filteren op TRA
          zonder <- str_remove_all(r$Meta,"\"")#verwijderen van \
          stuks <- strsplit(zonder,",")#splitten in categorien
          for (z in stuks[[1]]){ #ga over catgerien heen
            if (isTRUE(grepl('subject.id:', z)) == TRUE){#is het de categorie met subject.id
              complexidi <- i
              CDR3ai <- r$CDR3
              TRAVi <- r$V
              TRAJi <- r$J
              if(isTRUE(is.null(PCi)) == TRUE){
                PCi <- str_remove(z," subject.id: ")
              }else if(PCi != str_remove(z," subject.id: ")){
                print(c("ERROR",i))
              }
            }
          }
        } else if(r$Gene == "TRB"){#filteren op TRB
          zonder <- str_remove_all(r$Meta,"\"")#verwijderen van \
          stuks <- strsplit(zonder,",")#splitten in categorien
          for (z in stuks[[1]]){ #ga over catgerien heen
            if (isTRUE(grepl('subject.id:', z)) == TRUE){#is het de categorie met subject.id
              CDR3bi <- r$CDR3
              TRBVi <- r$V
              TRBJi <- r$J
              if(isTRUE(is.null(PCi)) == TRUE){
                PCi <- str_remove(z," subject.id: ")
              }else if(PCi != str_remove(z," subject.id: ")){
                print(c("ERROR",i))
              }
            }
          }
        }
      }
      nie <- add_row(nie, complex.id = i,CDR3b = CDR3bi, TRBV = TRBVi, TRBJ = TRBJi, CDR3a = CDR3ai, TRAV = TRAVi, TRAJ = TRAJi, PatientCounts = PCi)
      PCi <- NULL
    }
  }
}
proc.time() - ptm

library(tidyverse)
#maak een lijst met posities van elke subid
for(e in dabsubid$subid){
  if (e %in% klaar != TRUE){#kijk of sub in de lijst met gedane subss staat
    klaar <- c(klaar, e)#sub aan lijst met gedane subs toevoegen
    pla <- (which(dabsubid$subid == e))
    nierow <- c(e)
    for(p in pla){
      if(dabsubid[p,]$`MHC A` %in% nierow != TRUE){
        nierow <- c(nierow, dabsubid[p,]$`MHC A`)
      }
      if(dabsubid[p,]$`MHC B` %in% nierow != TRUE){
        nierow <- c(nierow, dabsubid[p,]$`MHC B`)
      }
    }
    write_lines(nierow, "D:/gliph-1.0/gliph/bin/HLA_tabel.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/HLA_tabel.txt" ,append = TRUE, sep = "")
    nierow <- NULL
  }
}
print(dabsubid[p,]$`MHC A`)

for(q in dabsubid$CDR3){
  if (q %in% qlaar != TRUE){#kijk of sub in de lijst met gedane subss staat
    qlaar <- c(qlaar, q)#sub aan lijst met gedane subs toevoegen
    qpla <- (which(dabsubid$CDR3 == q))
    nieqrow <- c(q)
    for(qp in qpla){
      if(dabsubid[qp,]$`MHC A` %in% nieqrow != TRUE){
        nieqrow <- c(nieqrow, dabsubid[qp,]$`MHC A`)
      }
    }
    write_lines(nieqrow, "D:/gliph-1.0/gliph/bin/HLA_DATAvv.txt" ,append = TRUE, sep="\t")
    write_lines("\n", "D:/gliph-1.0/gliph/bin/HLA_DATAvv.txt" ,append = TRUE, sep = "")
    nieqrow <- NULL
  }
}
qlaar <- NULL

nr <- 2534
afdr <- NULL
clus <- DATAz.convergence.groups
sclus <- strsplit(as.character(clus[nr,3])," ")#splitten in categorien
for (s in sclus){
  afdr <- c(afdr, s)
}

afdr <- c(nr,str_remove(clus[nr,2],"CRG-"),afdr)

write_lines(afdr, "D:/gliph-1.0/gliph/bin/cluster2534.txt" ,append = TRUE, sep="\t")
#18,2,139,111,33,113,591,457,1370,2534









#maak een tsv for HLA type
nierow <- c(e)
for(p in pla){
  nierow <- c(nierow, p)
}

write_lines(nierow, "D:/gliph-1.0/gliph/bin/HLA_tabel.txt" ,append = TRUE, sep="\t")
write_lines("\n", "D:/gliph-1.0/gliph/bin/HLA_tabel.txt" ,append = TRUE, sep = "")


write.table(nierow, file = "nieurow.txt", row.names = F , sep ="/t")

#maak een nieuwe column met subids
for(e in dab$Meta){
  verw <- str_remove_all(e,"\"")#verwijderen van \
  gesplit <- strsplit(verw,",")#splitten in categorien
  for (x in gesplit[[1]]){ #ga over catgerien heen
    if (isTRUE(grepl('subject.id:', x)) == TRUE){#is het de categorie met sub
      subid <- str_remove(x," subject.id: ")
      niecol <-c(niecol, subid)
    }
  }
}

dabsubid <- dab
dabsubid$subid <- niecol

nierow <-NULL
klaar <- NULL # maak een lege lijst voor gedane METAs/subs
niecol <- NULL

library(TCellPack)
PlotTCellPack(gliph = "D:/gliph-1.0/gliph/bin/DATAz-convergence-groups.txt", cell.data = cdb, legend = TRUE)
PlotTCellPack(gliph = gliph.example, cell.data = cdb, legend = TRUE)
PlotTCellPack(gliph = gliph.example, cell.data = cell.data.continuous.example, legend = TRUE)
cdb <- cell.data.discrete.example
cdb$cell <- 'AAAGATGAGGACTGGT-17'
gliphex <- "D:/gliph-1.0/gliph/bin/DATAz-convergence-groups.txt"



write.table(DATA, "D:/gliph-1.0/gliph/bin/DATAz.txt", sep="\t", row.names = FALSE, quote = FALSE)

done <- NULL # leeg de lijst met gedan ids
library(tidyverse) #import library
nie <- tibble(complex.id = NA,CDR3b = NA, TRBV = NA, TRBJ = NA, CDR3a = NA, TRAV = NA, TRAJ = NA, PatientCounts = NA)#maak t dataframe  
nie <- nie[FALSE,] #leeg het dataframe

#write.table(mydata, "c:/mydata.txt", sep="\t")

#nue <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA) #maak een dataframe
#names(nue) <- c("complex.id","CDR3b", "TRBV", "TRBJ", "CDR3a", "TRAV", "TRAJ", "PatientCounts") #geef kolommen titels
#nue <- nue[FALSE,] #leeg het dataframe

#print(str_remove(z," subject.id: "))#verwijder titel
#nie <- add_row(nie, complex.id = i,CDR3b = NA, TRBV = NA, TRBJ = NA, CDR3a = r$CDR3, TRAV = r$V, TRAJ = r$J, PatientCounts = z)
#nue <- rbind(nue, c(i, NA, NA, NA, r$CDR3, r$V, r$J, z))
