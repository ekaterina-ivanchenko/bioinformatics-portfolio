
# Проект: Анализ FASTA файлов: GC-состав и поиск генов (ORF)

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

# Настройка matplotlib для русского текста
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ============ НАСТРОЙКА ПУТЕЙ ============
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')

# Создаём папки, если их нет
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# Пути к файлам
INPUT_FILE = os.path.join(DATA_DIR, 'sequence.fasta')
OUTPUT_TARGET = os.path.join(DATA_DIR, 'potential_target.fasta')
OUTPUT_PROTEIN = os.path.join(DATA_DIR, 'translated_protein.fasta')
GC_PLOT = os.path.join(RESULTS_DIR, 'gc_content_plot.png')

# ============ ГЕНЕТИЧЕСКИЙ КОД ============
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# ============ ОСНОВНЫЕ ФУНКЦИИ ============

def read_fasta(filepath):
    """
    Читает FASTA файл и возвращает последовательность и заголовок
    """
    if not os.path.exists(filepath):
        print(f"[ОШИБКА] Файл не найден: {filepath}")
        return "", ""
    
    sequence = ""
    header = ""
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    break
                header = line[1:]
            elif line:
                sequence += line.upper()
    
    return sequence, header

def gc_content(seq):
    """Рассчитывает процент G и C нуклеотидов"""
    if not seq:
        return 0.0
    
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    return (g + c) / total * 100

def gc_content_sliding(seq, window_size=100, step=50):
    """
    Рассчитывает GC-состав скользящим окном
    """
    if len(seq) < window_size:
        return [], []
    
    gc_values = []
    positions = []
    
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        gc = (window.count('G') + window.count('C')) / window_size * 100
        gc_values.append(gc)
        positions.append(i + window_size // 2)
    
    return positions, gc_values

def find_open_reading_frames(seq, min_len=300):
    """
    Ищет открытые рамки считывания (ORF)
    Старт: ATG, Стоп: TAA, TAG, TGA
    """
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []
    
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i+3]
            
            if codon == start_codon:
                for j in range(i + 3, len(seq) - 2, 3):
                    if seq[j:j+3] in stop_codons:
                        orf_len = j + 3 - i
                        if orf_len >= min_len:
                            orfs.append((i, j + 3, seq[i:j+3]))
                        break
                i = j + 3 if 'j' in locals() else i + 3
                if 'j' in locals():
                    del j
            else:
                i += 3
    
    orfs.sort(key=lambda x: len(x[2]), reverse=True)
    return orfs

def reverse_complement(seq):
    """Возвращает обратную комплементарную цепь ДНК"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def translate_dna_to_protein(dna_seq):
    """Переводит ДНК в белок"""
    protein = ""
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        protein += aa
        if aa == '*':
            break
    return protein

def visualize_gc_content(positions, gc_values):
    """Визуализация GC-состава"""
    plt.figure(figsize=(12, 6))
    
    plt.plot(positions, gc_values, 'g-', linewidth=1.5, label='GC-состав')
    plt.axhline(y=50, color='r', linestyle='--', linewidth=1, label='Средний GC (50%)', alpha=0.7)
    plt.fill_between(positions, gc_values, 50, where=(np.array(gc_values) > 50),
                     color='green', alpha=0.2, label='GC > 50%')
    plt.fill_between(positions, gc_values, 50, where=(np.array(gc_values) < 50),
                     color='red', alpha=0.2, label='GC < 50%')
    
    plt.xlabel('Позиция в геноме (нуклеотиды)', fontsize=12)
    plt.ylabel('GC-состав (%)', fontsize=12)
    plt.title('Распределение GC-состава по геному', fontsize=14, fontweight='bold')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 100)
    plt.tight_layout()
    
    plt.savefig(GC_PLOT, dpi=150, bbox_inches='tight')
    print(f"\n[OK] График сохранён: {GC_PLOT}")
    plt.show()

def save_fasta(header, sequence, output_path, line_width=60):
    """Сохраняет последовательность в FASTA формате"""
    with open(output_path, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + "\n")
    print(f"[OK] Сохранено: {output_path}")

def generate_test_fasta(filepath, length=2000):
    """Генерирует тестовый FASTA файл"""
    import random
    
    print(f"[ПРЕДУПРЕЖДЕНИЕ] Файл не найден, создаю тестовый: {filepath}")
    
    bases = ['A', 'T', 'G', 'C']
    sequence = ''.join(random.choice(bases) for _ in range(length))
    
    # Добавляем ATG для ORF
    sequence = sequence[:500] + 'ATG' + sequence[503:1000] + 'TAA' + sequence[1003:]
    
    with open(filepath, 'w') as f:
        f.write(">test_sequence random genome\n")
        for i in range(0, len(sequence), 70):
            f.write(sequence[i:i+70] + "\n")
    
    print(f"[OK] Создан тестовый файл, длина: {length} нт")

# ============ ОСНОВНАЯ ПРОГРАММА ============

if __name__ == "__main__":
    print("=" * 60)
    print("АНАЛИЗ FASTA: GC-СОСТАВ И ПОИСК ГЕНОВ (ORF)")
    print("=" * 60)
    
    # 1. Проверяем входной файл
    if not os.path.exists(INPUT_FILE):
        generate_test_fasta(INPUT_FILE)
    
    # 2. Загрузка последовательности
    print(f"\nЗагрузка: {INPUT_FILE}")
    sequence, header = read_fasta(INPUT_FILE)
    
    if not sequence:
        print("[ОШИБКА] Не удалось загрузить последовательность")
        sys.exit(1)
    
    print(f"[OK] Заголовок: {header if header else 'без заголовка'}")
    print(f"[OK] Длина генома: {len(sequence)} нт")
    
    # 3. GC-состав
    print("\n" + "=" * 60)
    print("GC-СОСТАВ")
    print("=" * 60)
    
    gc_global = gc_content(sequence)
    print(f"Общий GC-состав: {gc_global:.2f}%")
    
    if gc_global > 60:
        print("[ИНФО] Высокий GC (>60%) — геном стабилен")
    elif gc_global > 40:
        print("[ИНФО] Нормальный GC (40-60%) — оптимально")
    else:
        print("[ИНФО] Низкий GC (<40%) — геном менее стабилен")
    
    # 4. Визуализация
    if len(sequence) >= 100:
        positions, gc_sliding = gc_content_sliding(sequence, window_size=min(100, len(sequence)//5), step=min(50, len(sequence)//10))
        if gc_sliding:
            visualize_gc_content(positions, gc_sliding)
    
    # 5. Поиск ORF
    print("\n" + "=" * 60)
    print("ПОИСК ГЕНОВ (ORF)")
    print("=" * 60)
    
    orfs_forward = find_open_reading_frames(sequence, min_len=150)
    print(f"На прямой цепи: {len(orfs_forward)} ORF")
    
    rev_seq = reverse_complement(sequence)
    orfs_reverse = find_open_reading_frames(rev_seq, min_len=150)
    print(f"На обратной цепи: {len(orfs_reverse)} ORF")
    
    all_orfs = orfs_forward + orfs_reverse
    all_orfs.sort(key=lambda x: len(x[2]), reverse=True)
    
    print(f"\n[OK] Всего найдено потенциальных генов: {len(all_orfs)}")
    
    for i, (start, end, orf_seq) in enumerate(all_orfs[:5]):
        print(f"  {i+1}. Длина: {len(orf_seq)} нт ({len(orf_seq)//3} а.к.) | позиции: {start}-{end}")
    
    # 6. Трансляция и сохранение
    if all_orfs:
        longest_orf = all_orfs[0][2]
        protein = translate_dna_to_protein(longest_orf)
        
        print("\n" + "=" * 60)
        print("ТРАНСЛЯЦИЯ В БЕЛОК")
        print("=" * 60)
        print(f"Самый длинный ORF: {len(longest_orf)} нт")
        print(f"Длина белка: {len(protein)} аминокислот")
        print(f"Начало белка: {protein[:50]}..." if len(protein) > 50 else f"Белок: {protein}")
        
        save_fasta(f"translated_protein_{len(protein)}aa", protein, OUTPUT_PROTEIN)
        save_fasta(f"potential_target_gene_{len(longest_orf)}nt", longest_orf, OUTPUT_TARGET)
    
    # 7. Итог
    print("\n" + "=" * 60)
    print("РЕЗУЛЬТАТЫ")
    print("=" * 60)
    print(f"  [1] {GC_PLOT}")
    print(f"  [2] {OUTPUT_TARGET}")
    print(f"  [3] {OUTPUT_PROTEIN}")
    
    print("\n" + "=" * 60)
    print("АНАЛИЗ ЗАВЕРШЁН")
    print("=" * 60)