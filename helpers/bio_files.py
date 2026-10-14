def write_fasta(output_fasta_file, header, part):
  
    output_fasta_file.write(header + "\n" + "".join(part) + "\n")
  

def write_proteins(output_file, list_of_proteins):
  
    with open(output_file, mode = "w") as out_file:
      out_file.write("\n".join(list_of_proteins) + "\n")
      

def parse_gbk_file(input_file) -> list:
  
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
  

def find_genes_indx(genes, cds_list):
  
    if isinstance(genes, str):
      genes = [genes]
      
    genes_indx = []
    for gene in genes:
        for i, cds in enumerate(cds_list):
            if gene == cds.get("/gene"):
              genes_indx.append(i)
    return set(genes_indx)

                
    
    
          
                    