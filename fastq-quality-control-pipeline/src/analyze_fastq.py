# Проект: Анализ качества FASTQ и фильтрация ридов перед выравниванием

import gzip
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# Настройка matplotlib для корректного отображения русского текста
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

# ============ НАСТРОЙКА ПУТЕЙ ============
# Получаем путь к корневой директории проекта (где находится папка src)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)  # поднимаемся на уровень выше

# Пути к папкам
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')

# Создаём папки, если их нет
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Пути к файлам
INPUT_FILE = os.path.join(DATA_DIR, 'sample.fastq')
OUTPUT_FILE = os.path.join(DATA_DIR, 'filtered_reads.fastq')
QUALITY_PLOT = os.path.join(RESULTS_DIR, 'quality_plot.png')
FILTERING_PLOT = os.path.join(RESULTS_DIR, 'filtering_stats.png')

# ============ ОСНОВНЫЕ ФУНКЦИИ ============

def parse_fastq(filepath):
    """
    Парсит FASTQ файл (обычный или сжатый .gz)
    Возвращает список ридов, где каждый рид = (заголовок, последовательность, качества)
    """
    reads = []
    
    # Проверяем, существует ли файл
    if not os.path.exists(filepath):
        print(f"[ОШИБКА] Файл не найден: {filepath}")
        return reads
    
    # Открываем файл (с поддержкой gzip, если нужно)
    if filepath.endswith('.gz'):
        f = gzip.open(filepath, 'rt')
    else:
        f = open(filepath, 'r')
    
    with f:
        lines = [line.strip() for line in f if line.strip()]
    
    # FASTQ: 4 строки на рид
    for i in range(0, len(lines), 4):
        if i+3 < len(lines):
            header = lines[i]
            sequence = lines[i+1]
            quality = lines[i+3]
            reads.append((header, sequence, quality))
    
    return reads

def filter_reads(reads, min_quality=25, min_length=50, min_percent_good=0.85):
    """
    Фильтрует риды по качеству и длине
    
    Параметры:
    - min_quality: минимальный Phred score (25 = 99.7% точность)
    - min_length: минимальная длина рида в нуклеотидах
    - min_percent_good: какой процент позиций должен быть выше min_quality
    """
    passed_reads = []
    stats = {
        'total': len(reads),
        'filtered_by_length': 0,
        'filtered_by_quality': 0,
        'passed': 0
    }
    
    for header, seq, qual in reads:
        # Фильтр по длине
        if len(seq) < min_length:
            stats['filtered_by_length'] += 1
            continue
        
        # Фильтр по качеству: считаем, сколько позиций с качеством >= порога
        good_positions = 0
        for char in qual:
            score = ord(char) - 33
            if score >= min_quality:
                good_positions += 1
        
        good_ratio = good_positions / len(qual) if len(qual) > 0 else 0
        
        if good_ratio >= min_percent_good:
            passed_reads.append((header, seq, qual))
            stats['passed'] += 1
        else:
            stats['filtered_by_quality'] += 1
    
    return passed_reads, stats

def analyze_position_quality(reads, max_length=100):
    """
    Анализирует среднее качество по каждой позиции в риде (от 1 до max_length)
    """
    position_scores = []
    
    # Собираем все последовательности качеств
    for _, _, qual in reads:
        scores = [ord(c) - 33 for c in qual[:max_length]]
        position_scores.append(scores)
    
    if not position_scores:
        return [], {}
    
    # Транспонируем и считаем среднее и перцентили
    avg_scores = []
    percentiles = {25: [], 50: [], 75: []}
    
    for pos in range(max_length):
        pos_values = [scores[pos] for scores in position_scores if pos < len(scores)]
        if pos_values:
            avg_scores.append(sum(pos_values) / len(pos_values))
            for p in percentiles:
                percentiles[p].append(np.percentile(pos_values, p))
        else:
            avg_scores.append(0)
            for p in percentiles:
                percentiles[p].append(0)
    
    return avg_scores, percentiles

def visualize_quality(avg_scores, percentiles):
    """
    Визуализация качества по позициям (график с перцентилями и порогами)
    """
    positions = range(1, len(avg_scores) + 1)
    
    plt.figure(figsize=(12, 6))
    
    # Извлекаем данные для перцентилей
    q25 = percentiles[25]
    q50 = percentiles[50]
    q75 = percentiles[75]
    
    # Зона между 25 и 75 перцентилями
    plt.fill_between(positions, q25, q75, alpha=0.3, color='blue', label='25-75 перцентиль')
    # Медиана
    plt.plot(positions, q50, 'b-', linewidth=2, label='Медиана (Q50)')
    # Среднее значение
    plt.plot(positions, avg_scores, 'r--', linewidth=1.5, label='Среднее', alpha=0.7)
    
    # Пороги качества
    plt.axhline(y=30, color='g', linestyle='--', label='Отлично (Q30)', alpha=0.7)
    plt.axhline(y=25, color='y', linestyle='--', label='Хорошо (Q25)', alpha=0.7)
    plt.axhline(y=20, color='orange', linestyle='--', label='Приемлемо (Q20)', alpha=0.7)
    
    # Подписи осей и заголовок
    plt.xlabel('Позиция в риде (нуклеотид)', fontsize=12)
    plt.ylabel('Phred Quality Score', fontsize=12)
    plt.title('Качество ридов по позициям', fontsize=14, fontweight='bold')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 45)
    plt.tight_layout()
    
    # Сохраняем график в папку results
    plt.savefig(QUALITY_PLOT, dpi=150, bbox_inches='tight')
    print(f"\n[OK] График качества сохранён как '{QUALITY_PLOT}'")
    plt.show()

def visualize_filtering_stats(stats):
    """
    Визуализация статистики фильтрации (столбчатая диаграмма)
    """
    plt.figure(figsize=(10, 5))
    
    # Данные для столбцов
    categories = ['Всего ридов', 'Прошли фильтр', 'Отсев по длине', 'Отсев по качеству']
    values = [stats['total'], stats['passed'], 
              stats['filtered_by_length'], stats['filtered_by_quality']]
    
    # Цвета для столбцов
    colors = ['#3498db', '#2ecc71', '#e74c3c', '#f39c12']
    bars = plt.bar(categories, values, color=colors)
    
    # Подписи осей и заголовок
    plt.ylabel('Количество ридов', fontsize=12)
    plt.title('Результаты фильтрации ридов', fontsize=14, fontweight='bold')
    
    # Добавляем значения на столбцы
    for bar, value in zip(bars, values):
        if value > 0:
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(value), ha='center', va='bottom', fontweight='bold')
    
    # Поворачиваем подписи для лучшей читаемости
    plt.xticks(rotation=15)
    plt.tight_layout()
    
    # Сохраняем график в папку results
    plt.savefig(FILTERING_PLOT, dpi=150, bbox_inches='tight')
    print(f"[OK] График статистики сохранён как '{FILTERING_PLOT}'")
    plt.show()

def save_filtered_reads(passed_reads, output_path):
    """Сохраняет отфильтрованные риды в FASTQ формате"""
    with open(output_path, 'w') as f:
        for header, seq, qual in passed_reads:
            f.write(f"{header}\n{seq}\n+\n{qual}\n")
    
    print(f"[OK] Сохранено {len(passed_reads)} ридов в {output_path}")

def generate_test_fastq(filepath, num_reads=100):
    """Генерирует тестовый FASTQ файл, если он не существует"""
    print(f"[ПРЕДУПРЕЖДЕНИЕ] {filepath} не найден, создаю тестовый файл...")
    
    with open(filepath, 'w') as f:
        for i in range(num_reads):
            f.write(f"@read_{i}\n")
            f.write("AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n")
            f.write("+\n")
            qualities = []
            for j in range(50):
                if j < 35:
                    q = 38
                else:
                    q = 20
                qualities.append(chr(q + 33))
            f.write(''.join(qualities) + "\n")
    
    print(f"[OK] Создан тестовый файл: {filepath}")

# ============ ОСНОВНАЯ ПРОГРАММА ============

if __name__ == "__main__":
    # Проверка наличия необходимых библиотек
    try:
        import matplotlib
        import numpy
    except ImportError:
        print("[ОШИБКА] Требуются библиотеки matplotlib и numpy")
        print("Установите их командой: pip install matplotlib numpy")
        sys.exit(1)
    
    print("=" * 60)
    print("АНАЛИЗ КАЧЕСТВА FASTQ ФАЙЛА")
    print("=" * 60)
    print(f"Корневая папка проекта: {PROJECT_ROOT}")
    print(f"Папка с данными: {DATA_DIR}")
    print(f"Папка с результатами: {RESULTS_DIR}")
    print()
    
    # 1. Проверяем/создаём входной файл
    if not os.path.exists(INPUT_FILE):
        generate_test_fastq(INPUT_FILE, num_reads=100)
    
    print(f"Загрузка {INPUT_FILE}...")
    reads = parse_fastq(INPUT_FILE)
    
    if not reads:
        print("[ОШИБКА] Не удалось загрузить риды из файла")
        sys.exit(1)
    
    print(f"[OK] Загружено ридов: {len(reads)}")
    
    # 2. Анализ качества по позициям
    avg_quality, percentiles = analyze_position_quality(reads, max_length=100)
    
    if avg_quality:
        print("\n" + "=" * 60)
        print("КАЧЕСТВО ПО ПОЗИЦИЯМ (средний Phred score)")
        print("=" * 60)
        
        # Показываем первые 20 позиций
        print("Первые 20 позиций:")
        for pos, q in enumerate(avg_quality[:20]):
            print(f"Позиция {pos+1:3d}: {q:.1f}")
        
        # Находим проблемные позиции с низким качеством
        low_quality_positions = [i+1 for i, q in enumerate(avg_quality) if q < 25]
        if low_quality_positions:
            print(f"\n[ПРЕДУПРЕЖДЕНИЕ] Позиции с низким качеством (<25): {low_quality_positions[:10]}")
            if len(low_quality_positions) > 10:
                print(f"         ... и ещё {len(low_quality_positions)-10} позиций")
        else:
            print("\n[OK] Все позиции имеют качество выше порога 25")
        
        # ГРАФИК 1: Качество по позициям
        visualize_quality(avg_quality, percentiles)
    
    # 3. Фильтрация ридов
    print("\n" + "=" * 60)
    print("ФИЛЬТРАЦИЯ РИДОВ")
    print("=" * 60)
    passed_reads, stats = filter_reads(reads, min_quality=25, min_length=50, min_percent_good=0.85)
    
    # Выводим статистику фильтрации
    print(f"Всего ридов:              {stats['total']}")
    print(f"Отфильтровано по длине:   {stats['filtered_by_length']}")
    print(f"Отфильтровано по качеству: {stats['filtered_by_quality']}")
    print(f"[OK] Прошло фильтрацию:   {stats['passed']}")
    print(f"Процент сохранённых данных: {stats['passed']/stats['total']*100:.1f}%")
    
    # ГРАФИК 2: Статистика фильтрации
    visualize_filtering_stats(stats)
    
    # 4. Сохранение очищенных данных
    if passed_reads:
        save_filtered_reads(passed_reads, OUTPUT_FILE)
    
    # 5. Детальный отчёт по первым ридам
    print("\n" + "=" * 60)
    print("ДЕТАЛЬНЫЙ ОТЧЁТ ПО РИДАМ (первые 10)")
    print("=" * 60)
    for i, (header, seq, qual) in enumerate(reads[:10]):
        good_positions = sum(1 for c in qual if (ord(c)-33) >= 25)
        good_ratio = good_positions / len(qual) if len(qual) > 0 else 0
        status = "ПРОШЁЛ" if (len(seq) >= 50 and good_ratio >= 0.85) else "ОТСЕВ"
        status_symbol = "[OK]" if status == "ПРОШЁЛ" else "[X]"
        print(f"{status_symbol} Рид {i+1:2d}: {status:7s} | длина={len(seq):3d} | хороших позиций={good_ratio*100:5.1f}%")
    
    # 6. Вывод рекомендаций для дальнейшего анализа
    print("\n" + "=" * 60)
    print("РЕКОМЕНДАЦИИ")
    print("=" * 60)
    retention_rate = stats['passed'] / stats['total'] * 100
    if retention_rate > 80:
        print("[OK] Качество данных хорошее. Можно выполнять выравнивание на референсный геном.")
    elif retention_rate > 50:
        print("[ПРЕДУПРЕЖДЕНИЕ] Качество данных среднее. Рекомендуется тримминг концов ридов.")
    else:
        print("[ОШИБКА] Качество данных низкое. Рекомендуется пересмотреть протокол секвенирования.")
    
    # Анализ падения качества к концу рида
    if len(avg_quality) > 20:
        start_quality = np.mean(avg_quality[:20])
        end_quality = np.mean(avg_quality[-20:])
        quality_drop = start_quality - end_quality
        if quality_drop > 10:
            print(f"[ИНФО] Качество падает на {quality_drop:.1f} баллов от начала к концу.")
            print("      Рекомендуется обрезка 3'-концов ридов.")
    
    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)
    print("\nСозданные файлы:")
    print(f"  [1] {QUALITY_PLOT}")
    print(f"  [2] {FILTERING_PLOT}")
    print(f"  [3] {OUTPUT_FILE}")