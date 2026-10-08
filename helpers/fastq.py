import math

# подсчет значения среднего качества рида
def average_qt(quality_seq) -> float:
    quality_sum = 0
    for char in quality_seq:
        q_score = ord(char) - 33
        p = 10 ** (-q_score / 10)
        quality_sum += p
    return -10 * math.log10(quality_sum / len(quality_seq))


# подсчет gc состава
def count_gc(seq) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) * 100 / len(seq)


# True, если процент GC в риде входит в интервал или меньше введенного значения
def check_gc(seq, gc_bounds) -> bool:
    gc_content = count_gc(seq)
    if isinstance(gc_bounds, (tuple, list)):
        return gc_content >= gc_bounds[0] and gc_content <= gc_bounds[1]
    return gc_content <= gc_bounds


# True, если значение среднего качества не ниже порога
def check_qt(quality_seq, quality_threshold) -> bool:
    avg_qt = average_qt(quality_seq)
    return avg_qt >= quality_threshold


# True, если длина рида входит в интервал или меньше введенного значения
def check_len(seq, length_bounds) -> bool:
    if isinstance(length_bounds, (tuple, list)):
        return len(seq) >= length_bounds[0] and len(seq) <= length_bounds[1]
    return len(seq) <= length_bounds