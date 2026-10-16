# bioinformatics-file-parsers

Небольшой пакет для работы с биоинформатическими файлами (FASTQ, FASTA, BLAST, GenBank). Сделан в рамках практикума по Python в биоинформатике. Все функции используют только стандартную библиотку Python.

## Установка

Нужен Python 3.9 или новее. Внешние библиотеки не потребуются.

```bash
git clone https://github.com/Samsagirl1/bioinf-file-parsers.git
cd bioinf-file-parsers
```

## Структура пакета

```
-/
|- README.md
|- главный скрипт.py # утилиты (импорты, filter_fastq и функции из части B)
|- helpers/
    |- fastq.py       # вспомогательные функции для filter_fastq
    |- bio_files.py   # вспомогательные функции для части B
```

## Основные функции

### `filter_fastq`
Фильтрует риды/fastq-последовательности по GC составу, длине нуклеотидной последовательности и по среднему качеству рида (шкала phred 33). Фильтрация проходит для каждого рида в момент прочтения входного файла. Результат сохраняется в папку "filtered". ч

```python
from main import filter_fastq

filter_fastq(input_fastq = "example_fastq.fastq", output_fastq = "filtered_fastq.fastq", gc_bounds = (40, 60), length_bounds = (50, 150), quality_threshold = 20)
```

- `gc_bounds` и `length_bounds` принимают интервал `(min, max)` или одно число, которое будет считаться верхней границей с нижней по-умолчанию равной 0.
- `quality_threshold` принимает число - минимальное среднее качество рида
- convert_multiline_fasta_to_oneline

### `convert_multiline_fasta_to_oneline` 
Принимает на вход fasta-файл и преобразует каждую последовательность таким образом, чтобы каждая последовательность (ДНК/РНК/Белка/...) умещалась в одной строке. В результате работы функции все преобразованные последовательности будут сохранены в новый fasta-файл

```python
from main import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline(input_fasta = "example_fasta.fasta", output_fasta = "output_fasta.fasta")
```

- Аргумент `output_fasta` необязятелен

### `parse_blast_output`
Берет из отчета работы алгоритма BLAST первое название из столбца `Description` для каждого QUERY запроса и сохраняет уникальные названия по алфавиту, по одному на строку.

```python
from main import parse_blast_output

parse_blast_output(input_file = "example_blast_result.txt", output_file = "parsed_blast_result.txt")
```

- Аргумент `output_file` необязателен

### `select_genes_from_gbk_to_fasta
Находит в GenBank-файле (.gbk) гены интереса и сохраняет в FASTA-файл аминокислотные последовательности (квалификатор /translation) их соседей (фланкирующих генов). Последовательности генов интереса не записываются в выходной файл.

```python
from main import select_genes_from_gbk_to_fasta

select_genes_from_gbk_to_fasta(input_gbk= "example_gbk.gbk", genes = ["sdhA", "sucC"], output_fasta = "output_fasta.fasta", n_before = 1, n_after = 1)
```
- `n_before`, `n_after` принимают количество генов до гена интереса и после. По-умолчанию равны 1
- `genes` принимает ген интереса или список генов интереса
- `output_fasta` необязателен

## Ограничения

- Работают только файлы в ожидаемых форматах (см. папку `example_data`)
- `parse_blast_output` читает только признаки `CDS`

## Автор

Валерия Киренская
Telegram: @samsagirl
