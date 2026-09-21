# ДЗ №4 — анализ Hi-C

[Отчёт](report/HW4_HiC_research_report.md)

## Данные

In situ Hi-C, MboI, GRCh38; мононуклеарные клетки костного мозга двух пациентов с B-ALL.

| Образец | Файл 4DN | Эксперимент |
|---|---|---|
| MCG037 | [4DNFIJ2JKO7D](https://data.4dnucleome.org/files-processed/4DNFIJ2JKO7D/) | [4DNES8BK27SC](https://data.4dnucleome.org/experiment-set-replicates/4DNES8BK27SC/) |
| MCG011 | [4DNFI8I9LN74](https://data.4dnucleome.org/files-processed/4DNFI8I9LN74/) | [4DNES2DB4F2D](https://data.4dnucleome.org/experiment-set-replicates/4DNES2DB4F2D/) |

Локус: **chr17:10,000,000–11,000,000**, разрешение **25 кб**. Insulation: окна 100, 200, 300, 500 и 1000 кб; основное окно — 200 кб. Альтернативный инструмент — OnTAD.

## Запуск

Python 3.12, `curl`, `git`, C++11-компилятор, libcurl и zlib.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python download_data.py
bash scripts/setup_ontad.sh
python scripts/analyze_opt.py
python scripts/alternative_tad.py
python scripts/verify.py
```

URL оригиналов и MD5 находятся в [data/sources.json](data/sources.json). Матрицы скачиваются отдельно и не включены в Git.

## Результаты

- [HiGlass: MCG037](report/higlass/MCG037_HiGlass.png), [MCG011](report/higlass/MCG011_HiGlass.png).
- [Таблицы и графики](results/): атрибуты Cooler, bins, raw/balanced контакты, CLI dump, P(s), insulation, ТАДы и сравнение с OnTAD.
- Границы: [MCG037](results/MCG037_boundaries.bed), [MCG011](results/MCG011_boundaries.bed).

BED6: `chrom start end name score strand`, координаты 0-based, полуоткрытые. `score = round(1000 × strength / (1 + strength))`; исходная сила границы — `boundary_strength_200000` в `*_insulation_locus.tsv`.
