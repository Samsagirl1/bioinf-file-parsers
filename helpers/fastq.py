import math

def uniform_input(gc_bounds: float | int | tuple[float, float], 
                  length_bounds: float | int | tuple[float, float], 
                  quality_threshold) -> tuple:
    """
    Converts all interval arguments to the form (min, max)

    Arguments:
    gc_bounds: float | int | tuple - GC content interval in percent. A single number is treated as the upper bound
    length_bounds: float | int | tuple - read length interval. A single number is treated as the upper bound
    quality_threshold: float | int - minimal average read quality (phred33)

    Returns tuple (gc_bounds, length_bounds, quality_threshold)
    """
    
    result = []
    if isinstance(gc_bounds, (tuple, list)):
        result.append(gc_bounds)
    else:
        result.append((0, gc_bounds))

    if isinstance(length_bounds, (tuple, list)):
        result.append(length_bounds)
    else:
        result.append((0, length_bounds))

    result.append(quality_threshold)
    return tuple(result)

    
def average_qt(quality_seq: str) -> float:
    """
    Computes the average quality (phred33) value of a read 

    Arguments:
    quality_seq: str - quality sequence of a read from FASTQ file

    Returns mean quality: float. For an empty line returns 0
    """
    
    if not quality_seq:
        return 0
        
    quality_sum = 0
    for char in quality_seq:
        q_score = ord(char) - 33
        p = 10 ** (-q_score / 10)
        quality_sum += p
    return -10 * math.log10(quality_sum / len(quality_seq))


def count_gc(seq : str) -> float:
    """
    Computes GC content in sequence

    Arguments:
    seq: str - RNA/DNA sequence

    Returns GC content: float. For an empty sequence returns 0
    """
    if not seq:
        return 0
        
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) * 100 / len(seq)


def check_gc(seq: str, 
             gc_bounds: tuple[float, float]) -> bool:
    """
    Checks that GC content of a sequence is within the bounds

    Arguments:
    seq: str - RNA/DNA sequence
    gc_bounds: tuple[float, float] - (min, max) GC content in percent

    Return bool. True if the read passes the check, False if not
    """
    
    gc_content = count_gc(seq)
    return gc_content >= gc_bounds[0] and gc_content <= gc_bounds[1]


def check_qt(quality_seq: str, 
             quality_threshold: float | int) -> bool:
    """
    Checks that mean read quality is not lower than the threshold

    Arguments:
    quality_seq: str - quality line of the read from FASTQ file
    quality_threshold: float | int - miniman mean quality (phred33)

    Returns bool. True if the read passes the check, False if not
    """
    
    avg_qt = average_qt(quality_seq)
    return avg_qt >= quality_threshold


def check_len(seq: str, 
              length_bounds: tuple[float, float]) -> bool:
    """
    Checks that sequence length is within the bounds

    Arguments:
    seq: str - RNA/DNA sequence
    lentgh_bounds: tuple[float, float] - (min, max) length 

    Returns bool. True of the read passes the check
    """
    return len(seq) >= length_bounds[0] and len(seq) <= length_bounds[1]


def check_read(read: list[str], 
               gc_bounds: tuple[float, float], 
               length_bounds: tuple[float, float], 
               quality_threshold: float | int) -> bool:
    """
    Applies all filters to ont FASTQ read

    Arguments:
    read: list[str] - four line of the read (header, sequence, "+", quality)
    gc_bounds: tuple[float, float] - GC content bounds
    length_bounds: tuple[float, float] - length bounds
    quality_threshold: float | int - minimal mean quality

    Returns True if the read passes all checks
    """
    
    return (check_len(read[1], length_bounds) and 
            check_qt(read[3], quality_threshold) and
            check_gc(read[1], gc_bounds))