import math

def uniform_input(gc_bounds, length_bounds, quality_threshold) -> tuple:
    """
    Приводит все границы к виду (min, max)
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
    Функция для подсчета значения среднего качества рида
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
    Функция для подсчета GC состава
    """
    if not seq:
        return 0
        
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) * 100 / len(seq)


def check_gc(seq: str, gc_bounds: float | int | tuple[float, float]) -> bool:
    """
    Функция для фильтрации по GC составу по заданному порогу
    """
    
    gc_content = count_gc(seq)
    return gc_content >= gc_bounds[0] and gc_content <= gc_bounds[1]


def check_qt(quality_seq: str, quality_threshold: float | int) -> bool:
    """
    Функция для фильтрации по значению среднего качества прочтения 
    """
    
    avg_qt = average_qt(quality_seq)
    return avg_qt >= quality_threshold


def check_len(seq: str, length_bounds: float | int | tuple[float, float]) -> bool:
    """
    Функция для фильтрации ридов по заданному интервалу длины 
    """
    return len(seq) >= length_bounds[0] and len(seq) <= length_bounds[1]


def check_read(read, gc_bounds, length_bounds, quality_threshold) -> bool:
    """
    Принимает список из 4 элементов (строк рида) и делает проверку
    """
    
    return (check_len(read[1], length_bounds) and 
            check_qt(read[3], quality_threshold) and
            check_gc(read[1], gc_bounds))