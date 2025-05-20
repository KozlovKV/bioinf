# *Козлов Кирилл*, Биоинформатика-1, Задание 3 - Построение пайплайна получения генетических вариантов
## Стэк
Minimap + Toil

Работаю на WSL 2

---

Взял геном человека WES - https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR14949306&display=download :
```bash
sudo apt install sra-toolkit -y # Для удобного скачивания fastq-файлов
fastq-dump ERR14949306 # ID из NCBI SRA
```
либо
```bash
wget -O ERR14949306.fastq.gz https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq\?acc\=ERR14949306
```

---

Скачиваем референсный геном:
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

---

Устанавливаем утилиты:
```bash
sudo apt update
sudo apt install fastqc samtools minimap2 -y
```

## Пайплайн в bash-скрипте
Файл [`pipeline.sh`](./pipeline.sh):
```bash
./pipeline.sh <seqence_file.fasta.gz> <reference_file.fa>
```

Вывод `fastqc` и `samtools flagstat` можно найти в [этой папке](./results/bash/)

## Toil
### Найстройка окружения
```bash
sudo apt update && sudo apt install graphviz -y # Для визуализации пайплайна

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Hello world
Код в файле [`./hello_world.py`](./hello_world.py)
```bash
source .venv/bin/activate
python3 hello_world.py file:hello-world-store
```

### Пайплайн
Код в файле [`./pipeline.py`](./pipeline.py)
```bash
source .venv/bin/activate
python3 pipeline.py file:<dir_for_tmp_files> <seqence_file.fasta.gz> <reference_file.fa> <output_dir>
# например
# python3 pipeline.py file:file_store seq.fastq.gz hg38.fa ./toil_pipeline_outputs --defaultMemory=8Gi --defaultCores=8 --defaultDisk=20Gi
```

FastQC-репорт, оценку выравнивания и логи пайплайна можно найти в [этой папке](./results/toil/)

Встроенное средство для построения графа пайплайна в Toil оказалаось нерабочим, поэтому не вижу смысла рисовать картинку вручную, так как по последовательности выполнения работы нет никаких отличий
- Единтсвенное - анализ исходного генома через **FastQC запускается параллельно** основному пайплайну, так как его результат нигде далее не задействуется