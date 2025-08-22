import os
import requests
from tqdm import tqdm
import tarfile

# --- Установка необходимых инструментов ---

# Клонирование репозиториев Kraken2, KrakenTools, SPAdes и Bracken
os.system('git clone https://github.com/DerrickWood/kraken2.git')
os.system('git clone https://github.com/jenniferlu717/KrakenTools.git')
#os.system('git clone https://github.com/ablab/spades.git')
os.system('git clone https://github.com/jenniferlu717/Bracken.git')

# Установка через conda (предполагается, что conda уже установлен)
# Обратите внимание: команда 'conda install' внутри скрипта может требовать активации окружения.
# Лучше создать отдельное окружение и активировать его перед запуском скрипта.
os.system('conda install -c bioconda bracken -y')
os.system('conda install -c bioconda kraken2 -y')

# Установка Python-библиотек requests и tqdm
os.system('pip install requests tqdm')


# --- Скачивание базы данных Kraken ---

# URL файла базы данных
url = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20250714.tar.gz"

# Имя файла для сохранения
local_filename = url.split('/')[-1]

# Функция для скачивания с прогрессбаром
def download_with_progress(url, filename):
    response = requests.get(url, stream=True)
    total_size_in_bytes= int(response.headers.get('content-length', 0))
    block_size = 1024 # 1 Kibibyte
    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
    with open(filename, 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()
    if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
        print("Ошибка при скачивании файла")
    else:
        print(f"Файл сохранен как {filename}")

# Скачиваем файл базы данных
print("Скачивание базы данных Kraken...")
download_with_progress(url, local_filename)

# Распаковка архива
print("Распаковка архива...")
with tarfile.open(local_filename, 'r:gz') as tar:
    tar.extractall()
print("Распаковка завершена.")

# Удаление исходного архива
os.remove(local_filename)
print(f"Удален архив {local_filename}")
