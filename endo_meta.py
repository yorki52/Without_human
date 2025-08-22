#!/usr/bin/env python
# coding: utf-8

import os
import glob
import subprocess
import pandas as pd
import matplotlib.pyplot as plt


N = get_ipython().getoutput('pwd ## имя рабочей машины')

kraken2_db = N[0] + "/Downloads/k2_standard_08gb_20250402/"
##установленная база данных для KRAKEN2


# Пути
fastq_dir = N[0]  + '/Downloads/250730/new/'
output_dir = fastq_dir + "/results"
# !mkdir output_dir

# kraken2_db = "/home/paul/Downloads/k2_standard_08gb_20250402/"

# Создаем папки для результатов
os.makedirs(os.path.join(output_dir, "kraken_reports"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "kraken_outputs"), exist_ok=True)

# Находим все файлы R1
r1_files = glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz"))

for r1_path in r1_files:
    # Получаем имя файла без пути
    r1_filename = os.path.basename(r1_path)
    # Формируем имя R2 файла по имени R1
    r2_filename = r1_filename.replace('_R1', '_R2')
    r2_path = os.path.join(fastq_dir, r2_filename)

    # Проверяем, существует ли R2 файл
    if not os.path.exists(r2_path):
        print(f"Файл пары не найден для {r1_filename}")
        continue

    # Имя образца (можно изменить по необходимости)
    sample_name = r1_filename.split('_')[0]

    print(f"Обработка: {r1_filename} и {r2_filename}")

    report_path = os.path.join(output_dir, "kraken_reports", f"{sample_name}_report.txt")
    output_path = os.path.join(output_dir, "kraken_outputs", f"{sample_name}_output.txt")

    cmd = [
        "kraken2",
        "--db", kraken2_db,
        "--paired",
        r1_path,
        r2_path,
        "--report", report_path,
        "--output", output_path
    ]

    # Запуск Kraken2
    subprocess.run(cmd)

print("Обработка завершена.")


fold_rep =  output_dir + "/kraken_reports"
os.chdir(fold_rep)


#для групп надо комбинировать репорты c помщью combine_kreports
comb = N[0]  + "/KrakenTools/combine_kreports.py"
get_ipython().system('python $comb -r *report.txt -o COMBINED.KREPORT')


files = get_ipython().getoutput('ls *COMBINED.KREPORT')
files



##BRACKEN

# Настройки
krona_db = kraken2_db
read_len = 150
threshold = 10

# Ваши файлы
files = get_ipython().getoutput('ls *COMBINED.KREPORT')
files

# Уровни таксонов, которые можно использовать
tax_levels = ["K", "P", "C", "O", "F", "G", "S"]

for file in files:
    sample_name = os.path.splitext(os.path.basename(file))[0]
    for level in tax_levels:
        # Запускаем для каждого уровня
        cmd_bracken = [
            "/home/paul/Bracken/bracken",
            "-d", krona_db,
            "-i", file,
            "-o", f"{sample_name}_{level}.bracken",
            "-r", str(read_len),
            "-l", level,
            "-t", str(threshold)
        ]
        print(f"Запуск: {' '.join(cmd_bracken)}")
        subprocess.run(cmd_bracken, check=True)
        print(f"Обработка файла {file} для уровня {level} завершена.")


# ##ПОИСК И УДАЛЕНИЕ ДАННЫх человеческих айди "9606","9443", "7711","33208","9605","9604","40674"!!!

# Путь к вашему скрипту фильтрации
script_path = N[0] + "/KrakenTools/filter_bracken.out.py"

files = get_ipython().getoutput('ls *.bracken                         # Ваши файлы с разными уровнями таксонов')
files

# Обработка каждого файла
for file in files:
    output_file = file.replace('.bracken', '_filtered.bracken')

    cmd = [
        "python3", script_path,
        "-i", file,
        "-o", output_file,
        "--exclude", "9606","9443", "7711","33208","9605","9604","40674",  # исключение человеческого ID
    ]

    print(f"Запуск: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"Файл {file} успешно обработан, результат: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Ошибка при обработке файла {file}: {e}")


# In[12]:


get_ipython().system('mkdir filt')


# In[13]:


cp *filtered.bracken filt/


# In[14]:


F = get_ipython().getoutput('pwd')
F = F[0] + "/filt"
os.chdir(F)


# In[ ]:





# In[15]:


bukvs = ["O","P", "F", "S", "G"]
bukvs

for bukva in bukvs:
    print(bukva)
    files = glob.glob("*"+bukva+"_filtered.bracken")  # или укажите конкретные имена

    all_samples = pd.DataFrame()

    for file in files:
        sample_name = file.replace(bukva + '_filtered.bracken', '')
        df = pd.read_csv(file, sep='\t')
        series = df.set_index('name')['fraction_total_reads']
        series.name = sample_name
        all_samples = pd.concat([all_samples, series], axis=1)

    # Заполняем пропуски нулями (если есть)
    all_samples.fillna(0, inplace=True)
    # Сохраняем итоговую таблицу
    all_samples.to_csv(bukva + '_combined_taxonomic_profiles.csv')


# In[16]:


# Список файлов для обработки
files = get_ipython().getoutput('ls *combined_taxonomic_profiles.csv')
files

# Функция для обработки каждого файла
def process_file(filename):
    print(f"Обработка файла: {filename}")
    df = pd.read_csv(filename, index_col=0)

    # Определяем образцы (столбцы)
    samples = list(df.columns)

    # Попытка определить контрольные образцы по названию
    control_samples = [s for s in samples if 'ctrl' in s.lower()]
    other_samples = [s for s in samples if s not in control_samples]

    # Порядок: контрольные в начале
    ordered_samples = control_samples + other_samples

    # Перестраиваем DataFrame по новому порядку
    df_ordered = df[ordered_samples]

    # Нормализация по каждому образцу
    df_normalized = df_ordered.div(df_ordered.sum(axis=0), axis=1)

    # Выбираем топ-10 видов по сумме долей
    species_sums = df_normalized.sum(axis=1)
    top10_species = species_sums.sort_values(ascending=False).head(10).index

    # Вырезаем только топ-10 видов
    df_top10 = df_normalized.loc[top10_species]

    # Транспонируем для построения графика (образцы - строки, виды - столбцы)
    df_top10_transposed = df_top10.T

    # Построение графика
    plt.figure(figsize=(15,10))
    ax = df_top10_transposed.plot(kind='bar', stacked=True)

    # plt.ylabel('Proportion')
    # plt.xlabel('Samples')

    title_name = os.path.splitext(os.path.basename(filename))[0]
    # plt.title(f'Top 10 in {title_name} (stacked bar 100%)')

    plt.legend( bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()

    # Сохраняем график или показываем
    output_image = f"{title_name}_top10_stacked_bar.png"
    plt.savefig(output_image)
    print(f"График сохранен как {output_image}")
    plt.show()
    plt.close()

# Обработка всех файлов
for file in files:
    process_file(file)

