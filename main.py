import os
from helpers import fastq
from helpers import bio_files

def filter_fastq(input_fastq: str, 
                 output_fastq: str = "filtered_fastq.fastq", 
                 gc_bounds: float | int | tuple[float, float] = (0, 100),
                 length_bounds: float | int | tuple[float, float] = (0, 2**32), 
                 quality_threshold: float | int = 0) -> None: 
    """   
    Filters reads in a fastq file by GC composition, read length and average quality. 
    Writes filtered read to a output FASTQ file

    Arguments:
    input_fastq: str - path to the input FASTQ file
    output_fastq: str - path to the output FASTQ file (optional)
    gc_bounds: float | int | tuple - GC content interval in percent. A single number is treated as the upper bound
    length_bounds: float | int | tuple - read length interval. A single number is treated as the upper bound
    quality_threshold: float | int - minimal average read quality (phred33)

    Returns None
    """ 
     
    os.makedirs("filtered", exist_ok=True)

    gc_bounds, length_bounds, quality_threshold = fastq.uniform_input(gc_bounds=gc_bounds, 
                                                                      length_bounds=length_bounds, 
                                                                      quality_threshold=quality_threshold)
    
    with open(input_fastq, mode = "r") as fastq_file, open(os.path.join("filtered", output_fastq), mode = "w") as out_fastq:
        read = []
        for line in fastq_file:
            read.append(line.rstrip("\n"))
            if len(read) == 4:
                if fastq.check_read(read = read, 
                                    gc_bounds=gc_bounds, 
                                    length_bounds=length_bounds, 
                                    quality_threshold=quality_threshold):
                    out_fastq.write("\n".join(read) + "\n")
                  
                read = []


def convert_multiline_fasta_to_oneline(input_fasta: str, 
                                       output_fasta: str = "oneline_fasta_out.fasta") -> None:
    """
    Rewrites a fasta file so that each sequence occupies a single line

    Arguments:
    input_fasta: str - path to the input fasta file
    output_fasta: str - path to the output file (optional)

    Returns None
    """
  
    with open(input_fasta, mode = 'r') as inp_fasta, open(output_fasta, mode = "w") as out_fasta:
        header = None
        part = []
      
        for line in inp_fasta:
            if line.startswith(">"):
                if header is not None:
                    bio_files.write_fasta(output_fasta_file = out_fasta, 
                                          header = header, 
                                          part = part)
                    part = []
                header = line.rstrip("\n")
            else:
                part.append(line.strip())

        if header is not None:
            bio_files.write_fasta(output_fasta_file = out_fasta, header = header, part = part)



def parse_blast_output(input_file: str, 
                       output_file: str = "top_matched_proteins.txt") -> None:
    """
    Selects the first protein's name from the "Description" column for each QUERY for a BLAST text report. 
    The function sorts the list of selected names in alphabetical order and writes it to a new file.
    
    Arguments: 
    input_file: str - path to the input file with blast results 
    output_file: str - path to the output file (optional)

    Returns None
    """
    with open(input_file, mode = "r") as inp_file:
        proteins = []
        description = False
        description_width = None
        for line in inp_file:
            if description:
                proteins.append(line[:description_width].strip())
                description = False
            if line.lower().startswith("description"):
                description = True
                description_width = line.index("Name")
        list_of_proteins = sorted(proteins, key = str.lower)
        bio_files.write_proteins(output_file = output_file, list_of_proteins = list_of_proteins)



def select_genes_from_gbk_to_fasta(input_gbk: str, 
                                   genes: str | list[str], 
                                   output_fasta: str = "flanking_genes.fasta", 
                                   n_before: int = 1, 
                                   n_after: int = 1):
    """
    Finds genes of interest in a GenBank file and writes protein sequences of their neighbours to a FASTA file. 
    The genes of interest themselves are normally not written, but if you write neighboring genes in genes argument,
    they will be written.

    Arguments:
    input_gbk: str - path to the input GenBank file
    genes: str | list[str] - gene name or collection of gene names of interest
    output_fasta: str - path to the output FASTA file (optional)
    n_before: int - number of genes before each gene of interest that need to be writen to a FASTA file
    n_after: int - number of genes after each gene of interest that need to be written to a FASTA file

    Returns None
    """
    cds_list = bio_files.parse_gbk_file(input_file = input_gbk)
    genes_indx = bio_files.find_genes_indx(genes = genes, cds_list = cds_list)

    with open(output_fasta, mode = "w") as out_fasta:
        for indx in sorted(genes_indx):
            flanking_genes = cds_list[max(indx - n_before, 0):indx] + cds_list[indx + 1:indx + n_after + 1]
            for gene in flanking_genes:
                if not gene.get("/translation"):
                    continue
                gene_name = gene.get("/gene") or gene.get("/locus_tag")
                out_fasta.write(">" + gene_name + "\n" + gene["/translation"] + "\n")


  