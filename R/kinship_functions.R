# algorithm for calculating kinship coefficients

# pedigree example from Lange
# 1-----2
#    |
# -------
# |     |
# 3-----4
#    |
# -------
# |     |
# 5     6


#' Scanped Function
#' 
#' Reads in a pedigree from a specified input file. 
#' If parameter header=TRUE the first line should contain a header line of the form:
#' Indiv  Father  Mother  Sex
#' Each subsequent line should contain 3 integer indexes (individual, father, mother) and a letter M or F (sex)
#' All individuals must have an integer label larger than that of any ancestor and the k founders should have labels 1..k
#' @param filename
#' @param header. Defaults to TRUE
#' @export
#' @examples
#' scanped()


scanped <- function(filename)
{
  dat <- read.table(filename,sep="",header=TRUE)
  msg <- "test"
#  msg <- paste("Scanned",dt,"data for",nrow(fcsv),"individuals.")
  print(msg)
  return(dat)
}

#' Kinship function
#' 
#' Generates a n*n matrix of kinship coefficients for n individuals of a pedigree
#' @param pedigree
#' @export
#' @examples
#' kinship()
#' 
 
kinship <- function(pedigree)
{
  father=2
  mother=3
  noind=length(pedigree$ind)
  kmat <- matrix(data=rep(0,noind^2),nrows<-noind,ncols<-noind)
  for(i in 1:noind)
    for(j in i:noind)
    {
      if(i==j)
      {
        if(pedigree[i,father]==0)
          kmat[i,i]=0.5
        else
          if(pedigree[i,father] < pedigree[i,mother])
            kmat[i,i]=0.5+0.5*kmat[pedigree[i,father],pedigree[i,mother]]
          else
            kmat[i,i]=0.5+0.5*kmat[pedigree[i,mother],pedigree[i,father]]
      }
      else
      {
          if(pedigree[j,father]>0)
          {
            if((pedigree[j,father]>=i)&&(pedigree[j,mother]>=i))
              kmat[i,j]=0.5*kmat[i,pedigree[j,father]]+0.5*kmat[i,pedigree[j,mother]]
            else
              if((pedigree[j,father]>=i)&&(pedigree[j,mother]<i))
                kmat[i,j]=0.5*kmat[i,pedigree[j,father]]+0.5*kmat[pedigree[j,mother],i]
              else
                if((pedigree[j,father]<i)&&(pedigree[j,mother]>=i))
                  kmat[i,j]=0.5*kmat[pedigree[j,father],i]+0.5*kmat[i,pedigree[j,mother]]
                else
                  if((pedigree[j,father]<i)&&(pedigree[j,mother]<i))
                    kmat[i,j]=0.5*kmat[pedigree[j,father],i]+0.5*kmat[pedigree[j,mother],i]
          }
      }
      kmat[j,i]=kmat[i,j]
    }
  
  return(kmat)  
}

#' CoeffInbr function
#' 
#' Generates a vector of length n of inbreeding coefficients for n individuals of a pedigree
#' @param pedigree
#' @export
#' @examples
#' coeffInbr()
#' 

coeffInbr <- function(pedigree)
{
  father=2
  mother=3
  CoI <- rep(0,length(pedigree$ind))
  kmat <- kinship(pedigree)
  for(i in 1:length(pedigree$ind))
    if(pedigree[i,father]==0 || pedigree[i,mother]==0)
      CoI[i]=0
    else
      CoI[i] = kmat[pedigree[i,father],pedigree[i,mother]]
  return(CoI)
}

