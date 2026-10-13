import os
from helpers import fastq
from helpers import bio_files

def filter_fastq(input_fastq: str, 
                 output_fastq: str = "fitered_fastq.fastq", 
                 gc_bounds: float | int | tuple[float, int] = (0, 100),
                 length_bounds: float | int | tuple[float, int] = (0, 2**32), 
                 quality_threshold: float | int = 0) -> None: 
    """   
    Функция для фильтрации ридов в fastq файле по GC составу, длине рида, среднему   значению качества прочтения
    """ 
     
    os.makedirs("filtered", exist_ok=True)

    gc_bounds, length_bounds, quality_threshold = fastq.uniform_input(gc_bounds=gc_bounds, length_bounds=length_bounds, quality_threshold=quality_threshold)
    
    with open(input_fastq, mode = "r") as fastq_file, open(os.path.join("filtered", output_fastq), mode = "w") as out_fastq:
        read = []
        for line in fastq_file:
            read.append(line.strip())
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
                
            if not line.startswith(">"):
                part.append(line.rstrip("\n"))

        if header is not None:
            bio_files.write_fasta(output_fasta_file = out_fasta, header = header, part = part)
