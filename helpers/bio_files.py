def write_fasta(output_fasta_file, header, part):
    output_fasta_file.write(header + "\n" + "".join(part) + "\n")