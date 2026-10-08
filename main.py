import os
from helpers import fastq

def filter_fastq(input_fastq: str, 
                 output_fastq: str = "fitered_fastq.fastq", 
                 gc_bounds: float | int | tuple[float, int] = (0, 100),
                 length_bounds: float | int | tuple[float, int] = (0, 2**32), 
                 quality_threshold: float | int = 0) -> None: 
  """
  Функция для фильтрации ридов в fastq файле по GC составу, длине рида, среднему значению качества прочтения
  """
  
    os.makedirs("filtered", exist_ok=True)
    with open(input_fastq, mode = "r") as fastq_file, open(os.path.join("filtered", output_fastq), mode = "w") as out_fastq:
        read = []
        for line in fastq_file:
            read.append(line.strip())
            if len(read) == 4:
                if (fastq.check_len(read[1], length_bounds) and 
                    fastq.check_qt(read[3], quality_threshold) and 
                    fastq.check_gc(read[1], gc_bounds)):
                  
                  out_fastq.write("\n".join(read) + "\n")
                  
                read = []
