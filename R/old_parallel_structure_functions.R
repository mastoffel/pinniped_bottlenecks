# this is the old function from the parallelstructure package to run the code without breaking.

runsToDfStructure <- function(files = NULL)
{
    if (is.null(files) | (length(files) == 0)) stop("runsToDfStructure: No input files.\n")
    #number of files selected
    number <- length(files)
    #cat(paste("Number of files selected: ", number, "\n", sep = ""))
    
    #check file
    if (any(checkRuns(files)$type != "STRUCTURE")) stop("runsToDfStructure: Input contains one or more non-STRUCTURE files/Incorrect input format.\n")
    
    i <- 1
    dlist <- vector("list",length = number)
    len1 <- length(files)
    for (i in 1:len1)
    {
        name <- basename(files[i]) 
        file1 <- readLines(as.character(files[i]), warn = FALSE)
        
        #find individuals and get number of individuals
        ind <- as.numeric(as.character(base::gsub("\\D", "", grep("\\d individuals", file1, perl = TRUE, ignore.case = TRUE, value = TRUE)[1])))
        if (is.na(ind)) cat(paste0("Number of individuals is NA in file: ", name))
        
        #get value of k & error check
        k <- as.numeric(as.character(base::gsub("\\D", "", grep("\\d populations assumed", file1, perl = TRUE, ignore.case = TRUE, value = TRUE)[1])))
        if (is.na(k)) cat(paste0("Value of K is NA in file: ", name))
        
        file1 <- file1[grep(".+\\(\\d+\\).+\\:.+",file1)]
        if(length(file1) == 0)
        {
            cstart <- base::charmatch("Inferred ancestry of individuals", file1)
            cend <- base::charmatch("Estimated Allele Frequencies in each", file1)
            file1 <- file1[(cstart+2):(cend-1)]
        }
        
        file_a <- file1[file1 != ""]
        rm(file1)
        file_a <- base::gsub("\\([0-9.,]+\\)","",file_a)
        file_b <- base::gsub(":  ", "", substr(file_a, base::regexpr(":\\W+\\d\\.\\d+", file_a), base::nchar(file_a)-1))
        file_b <- base::sub("\\s+$","",base::sub("^\\s+","",file_b))
        rm(file_a)
        file_c <- as.vector(as.numeric(as.character(unlist(base::strsplit(file_b, " ")))))
        rm(file_b)
        dframe <- as.data.frame(matrix(file_c, nrow = ind, byrow = TRUE),stringsAsFactors = FALSE)
        dframe <- as.data.frame(sapply(dframe, as.numeric),stringsAsFactors = FALSE)
        colnames(dframe) <- paste0("Cluster", 1:k)
        dlist[[i]] <- dframe
        #names(dlist[[i]]) <- as.character(name)
    }
    if (number>1) {return(dlist)} else{return(dframe)}
}


checkRuns <- function(files=NULL, warn=FALSE)
{
    if (is.null(files)) stop("checkRuns: Input is empty.\n")
    len1 <- length(files)
    
    
    checkvec <- rep("UNIDENTIFIED",length=len1)
    subtype <- rep(NA,length=len1)
    for(i in 1:len1)
    {
        chk <- FALSE
        read1 <- readLines(files[i], n=7, warn = FALSE)
        
        #read TESS file
        chk <- grepl("ESTIMATED CLUSTERING PROBABILITIES", toupper(read1)[1])
        if (chk)
        {
            checkvec[i] <- "TESS"
        }
        
        #read STRUCTURE file
        if (!chk)
        {
            chk <- grepl("STRUCTURE BY PRITCHARD", toupper(read1)[4])
            if (chk)
            {
                checkvec[i] <- "STRUCTURE"
            }
        }
        rm(read1)
        
        #read MATRIX/TAB files
        if (!chk)
        {
            seps <- c("","\t",",")
            subtypes <- c("SPACE","TAB","COMMA")
            k=1
            while(!chk)
            {
                if(class(try(suppressWarnings(read.table(files[i],header=FALSE,sep=seps[k],nrows=1,quote="",stringsAsFactors = FALSE))))!="try-error")
                {
                    df <- read.table(files[i],header=FALSE,sep=seps[k],nrows=1,quote="",stringsAsFactors = FALSE)
                    if(all(sapply(df, is.numeric))) {
                        checkvec[i] <- "MATRIX"
                        subtype[i] <- subtypes[k]
                        chk <- TRUE
                    }else{
                        if((ncol(df) > 2) && (is.character(df[,1])))
                        {
                            checkvec[i] <- "TAB"
                            chk <- TRUE
                        }
                    }
                }
                k=k+1
                if(k>3)
                {
                    break
                }
            }
        }
        if((!chk) && warn) warning(paste0("checkRuns: ",files[i]," is not a STRUCTURE, TESS, MATRIX or TAB file.\n"))
    }
    return(list(type=checkvec,subtype=subtype))
}
