def write_fasta(output_fasta_file, 
                header: str, 
                part: list[str]) -> None:
    """
    Writes one FASTA record with the sequence on a single line

    Arguments:
    output_fasta_file - opened file for writing
    header: str - header line starting with ">"
    part: list[str] - lines of the sequence 

    Returns None
    """
  
    output_fasta_file.write(header + "\n" + "".join(part) + "\n")
  

def write_proteins(output_file: str, 
                   list_of_proteins: list[str]) -> None:
    """
    Writes sorted alphabetically protein names to the output file, one per line

    Arguments:
    output_file: str - path to the output file
    list_of_proteins: list[str] - list of proteins names

    Returns None
    """
  
    with open(output_file, mode = "w") as out_file:
        out_file.write("\n".join(list_of_proteins) + "\n")
      

def parse_gbk_file(input_file: str) -> list[dict[str, str]]:
    """
    Reads a GenBank file (.gbk) and collects qualifiers of every CDS 

    Arguments: 
    input_file: str - path to the GenBank file

    Returns list[dict[str, str]]. A list of dictionaries, where keys are qualifier names with the "/"
    """
  
    with open(input_file, mode = "r") as inp_file:
        is_cds = False
        cds = {}
        all_cds = []
        last_key = None
        header_indent = " " * 5
        body_indent   = " " * 21
      
        for line in inp_file:
            if is_cds:
                if line.startswith(body_indent):
                    if line.strip().startswith("/"):
                        line_type, line_info = line.strip().split("=", 1)
                        last_key = line_type.strip('"')
                        cds[last_key] = line_info.strip('"')
                    else:
                        cds[last_key] += line.strip().strip('"')
                    
            if line.startswith(header_indent) and not line.startswith(body_indent):
                if "CDS" in line:
                    is_cds = True
                    if len(cds) != 0:
                        all_cds.append(cds)
                        cds = {}
                        last_key = None
                else:
                    is_cds = False
        all_cds.append(cds)
        return all_cds
  

def find_genes_indx(genes: str | list[str], 
                    cds_list: list[dict[str, str]]) -> set[int]:
    """
    Finds position of the given genes of interest in the list of all CDS

    Arguments:
    genes: str | list[str] - gene name or list of gene names
    cds_list: list[dict[str, str]] - a list containing all CDS, each CDS is a dict in which the keys are qualifiers.

    Returns set[int] - set of indexes 
    """
    if isinstance(genes, str):
        genes = [genes]
      
    genes_indx = []
    for gene in genes:
        for i, cds in enumerate(cds_list):
            if gene == cds.get("/gene"):
                genes_indx.append(i)
    return set(genes_indx)

                
    
    
          
                    