
# Calculate Monoisotopic mass ---------------------------------------------


# This function was adapted from BRAIN package


calculate_monoisotopic_mass <- function(seq, IAA = TRUE, AcK = FALSE) {
  
  # Split sequence by amino acids
  seq <- toupper(seq)
  seq <- gsub(pattern = "[[:space:]]+", replacement = "", x = seq)
  seq <- strsplit(x = seq, split = "")
  
  # Create the weight scale
  weight <- c(
    A = 71.037113805,
    R = 156.101111050,
    N = 114.042927470,
    D = 115.026943065,
    C = ifelse(IAA == FALSE, 103.009184505, 103.009184505 + 57.021464),
    E = 129.042593135,
    Q = 128.058577540,
    G = 57.021463735,
    H = 137.058911875,
    I = 113.084064015,
    L = 113.084064015,
    K = ifelse(AcK == FALSE, 128.094963050, 128.094963050 + 42.010565),
    M = 131.040484645,
    F = 147.068413945,
    P = 97.052763875,
    S = 87.032028435,
    T = 101.047678505,
    W = 186.079312980,
    Y = 163.063328575,
    V = 99.068413945,
    H2O = 18.01056)
  
  # Sum the weight of each amino acid and add H2O weight
  unlist(lapply(seq, function(seq){
    sum(weight[c(seq, "H2O")])
  }))
}



# digest protein sequence -------------------------------------------------

# This function was adapted from BRAIN package


digest_aa_sequence <- function(sequence, 
                               enzyme = "trypsin", 
                               missed = 0){
  
  
  ## determine cleavage sites according to enzyme specific rules
  
  seq_vector <- strsplit(sequence, split = "")[[1]]
  end_position <- length(seq_vector)
  
  if(enzyme == "trypsin") {                                
    if(seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    } else seq_string <- sequence    
    seq_string <- gsub("KP", "!P", seq_string)          # prevent cleavage at K and R if followed by P  
    seq_string <- gsub("RP", "!P", seq_string)     
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1    
  }
  
  if(enzyme == "trypsin.strict") {                                
    if(seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    } else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1    
  }
  
  if(enzyme == "pepsin") {
    if(seq_vector[end_position] == "F" | seq_vector[end_position] == "L" |
       seq_vector[end_position] == "W" | seq_vector[end_position] == "Y" | 
       seq_vector[end_position] == "A" | seq_vector[end_position] == "E" | 
       seq_vector[end_position] == "Q") {
      seq_vector[end_position] <- "!"
    } 
    stop <- grep("F|L|W|Y|A|E|Q", seq_vector)       
    start <- stop + 1
  }
  
  if(enzyme == "arg.c") {                                
    if(seq_vector[end_position] == "R"){
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    } else seq_string <- sequence    
    seq_string <- gsub("RP", "!P", seq_string)          # prevent cleavage at R if followed by P  
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("R", seq_vector)
    start <- stop + 1    
  }
  
  if(enzyme == "arg.c, glu.c") {                                
    if(seq_vector[end_position] == "R" | seq_vector[end_position] == "E") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    } else seq_string <- sequence    
    seq_string <- gsub("RP", "!P", seq_string)          # prevent cleavage at R if followed by P  
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("R|E", seq_vector)
    start <- stop + 1    
  }
  ## error checking
  
  if(enzyme != "trypsin" & enzyme != "trypsin.strict" & enzyme != "pepsin" & enzyme != "arg.c" & enzyme != "arg.c, glu.c") 
    stop("undefined enzyme, defined enzymes are trypsin, trypsin.strict, and pepsin and arg.c, and glu.c")                    
  if(length(stop) == 0) warning("sequence does not contain cleavage sites")
  if(missed > length(stop)) stop("number of specified missed cleavages is greater than the maximum possible")
  
  
  ## cleave sequence
  
  cleave <- function(sequence, start, stop, misses) {
    peptide <- substring(sequence, start, stop)
    mc <- rep(misses, times = length(peptide))
    result <- data.frame(peptide, start, stop, mc, stringsAsFactors = FALSE)
    return(result) 
  }
  
  start <- c(1, start)                           # peptides if 0 missed cleavages
  stop <- c(stop, end_position)
  results <- cleave(sequence, start, stop, 0)
  
  if(missed > 0) {                               # peptides if missed cleavages > 0
    for(i in 1:missed) {
      start_tmp <- start[1:(length(start) - i)]
      stop_tmp <- stop[(1 + i):length(stop)]
      peptide <- cleave(sequence, start_tmp, stop_tmp, i)
      results <- rbind(results, peptide) 
    } 
  }
  return(results)
}



# import fasta as dataframe -----------------------------------------------

# Takes a fasta file input and generates a dataframe.
# The column names are:
  # ProteinAccession
  # ProteinDescription
  # GeneName
  # Organism
  # ProteinSequence
# The names of the column are formatted for using with EncyclopeDIA data output



import_fasta_as_df <- function(fname){
  
  
  # Error handling ----------------------------------------------------------
  
  
  # File must be a fasta
  if(!grepl(".fasta",fname)){
    stop(paste0("Input file must be a fasta"))
  }
  
  
  # Meat and potatoes -------------------------------------------------------
  
  
  
  
  require(seqinr)
  # Reading in the fasta file
  fasta <- read.fasta(fname,
                      seqtype = "AA",
                      set.attributes = FALSE, 
                      whole.header = TRUE)
  
  # Creating list for Protein Accession
  # This ID matches what is exported from EncyclopeDIA
  ProteinAccession <- unlist(lapply(names(fasta), function(x) unlist(strsplit(x, split = " "))[1]))
  
  # Parsing the Protein Description
  ProteinDescription <- substr(names(fasta), start = regexpr(pattern = " ", names(fasta)) + 1, stop = nchar(names(fasta)))
  ProteinDescription <- unlist(lapply(ProteinDescription, function(x) unlist(strsplit(x, split = " OS="))[1]))
  
  # Parsing the Gene name
  GeneName <- unlist(lapply(names(fasta), function(x) unlist(strsplit(x, split = "GN="))[2]))
  GeneName <- unlist(lapply(GeneName, function(x) unlist(strsplit(x, split = " PE="))[1]))
  
  # Parsing Organism
  Organism <- substr(ProteinAccession, start = regexpr("_", ProteinAccession) + 1, stop = nchar(ProteinAccession))
  
  # Creating a dataframe from the extracted information
  fasta_df <- data.frame(cbind(ProteinAccession, ProteinDescription, GeneName, Organism))
  
  # Protein Sequence
  fasta_df$ProteinSequence <- sapply(fasta, FUN = paste, collapse = "")
  
  return(fasta_df)
  
}

