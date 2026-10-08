import math


def average_qt(quality_seq: str) -> float:
  """
  Функция для подсчета значения среднего качества рида
  """
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
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) * 100 / len(seq)


def check_gc(seq: str, gc_bounds: float | int | tuple[int, float]) -> bool:
  """
  Функция для фильтрации по GC составу по заданному порогу
  """
    gc_content = count_gc(seq)
    if isinstance(gc_bounds, (tuple, list)):
        return gc_content >= gc_bounds[0] and gc_content <= gc_bounds[1]
    return gc_content <= gc_bounds


def check_qt(quality_seq: str, quality_threshold: float | int) -> bool:
  """
  Функция для фильтрации по значению среднего качества прочтения 
  """
    avg_qt = average_qt(quality_seq)
    return avg_qt >= quality_threshold


def check_len(seq: str, length_bounds: float | int | tuple[float, int]) -> bool:
  """
  Функция для фильтрации ридов по заданному интервалу длины 
  """
    if isinstance(length_bounds, (tuple, list)):
        return len(seq) >= length_bounds[0] and len(seq) <= length_bounds[1]
    return len(seq) <= length_bounds