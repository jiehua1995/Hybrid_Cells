# Read fasta to 2bit
folder <- "H:\\HybridCells_Data\\FigS2_Assembly of Ras3 genome\\"
fasta <- paste0(folder,"Ras3.fasta")
destfile <- paste0(folder,"Ras3.2bit")


# Build package
library(BSgenomeForge)
fastaTo2bit(fasta, destfile, assembly_accession = NA)


forgeBSgenomeDataPkgFromTwobitFile(destfile, organism="Drosophila melanogaster", provider="JieHua", genome="Ras3",pkg_maintainer="Jie Hua, <Jie.Hua@lmu.de>", pkg_version="1.0.0", pkg_license="Artistic-2.0",destdir="C:/temp/BSgenome_Ras3")