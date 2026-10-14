import os
from helpers import fastq
from helpers import bio_files

def filter_fastq(input_fastq: str, 
                 output_fastq: str = "filtered_fastq.fastq", 
                 gc_bounds: float | int | tuple[float, float] = (0, 100),
                 length_bounds: float | int | tuple[float, float] = (0, 2**32), 
                 quality_threshold: float | int = 0) -> None: 
    """   
    Функция для фильтрации ридов в fastq файле по GC составу, длине рида, среднему   значению качества прочтения
    """ 
     
    os.makedirs("filtered", exist_ok=True)

    gc_bounds, length_bounds, quality_threshold = fastq.uniform_input(gc_bounds=gc_bounds, length_bounds=length_bounds, quality_threshold=quality_threshold)
    
    with open(input_fastq, mode = "r") as fastq_file, open(os.path.join("filtered", output_fastq), mode = "w") as out_fastq:
        read = []
        for line in fastq_file:
            read.append(line.rstrip("\n"))
            if len(read) == 4:
                if fastq.check_read(read = read, gc_bounds=gc_bounds, length_bounds=length_bounds, quality_threshold=quality_threshold):
                    out_fastq.write("\n".join(read) + "\n")
                  
                read = []


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = "oneline_fasta_out.fasta"):
    with open(input_fasta, mode = 'r') as inp_fasta, open(output_fasta, mode = "w") as out_fasta:
        header = None
        part = []
      
        for line in inp_fasta:
            if line.startswith(">"):
                if header is not None:
                    bio_files.write_fasta(output_fasta_file = out_fasta, header = header, part = part)
                    part = []
                header = line.rstrip("\n")
            else:
                part.append(line.strip())

        if header is not None:
            bio_files.write_fasta(output_fasta_file = out_fasta, header = header, part = part)



def parse_blast_output(input_file, output_file):
    with open(input_file, mode = "r") as inp_file:
        proteins = []
        description = False
        for line in inp_file:
            if description:
                proteins.append(line.split("  ")[0])
                description = False
            if line.lower().startswith("description"):
                description = True
        list_of_proteins = sorted(proteins, key = str.lower)
        bio_files.write_proteins(output_file = output_file, list_of_proteins = list_of_proteins)



def select_genes_from_gbk_to_fasta(input_gbk, genes, output_fasta = "flanking_genes.fasta", n_before = 1, n_after = 1):
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


  