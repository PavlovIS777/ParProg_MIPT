import random

# Параметры генерации
NUM_ELEMENTS = 5000000  # 1 миллион элементов
MIN_VAL = -1000000
MAX_VAL = 1000000

# Генерация случайного массива
random_data = [random.randint(MIN_VAL, MAX_VAL) for _ in range(NUM_ELEMENTS)]

# Запись в файл
with open('input_data.txt', 'w') as f:
    f.write(f"{NUM_ELEMENTS} ")  # Первая строка - количество элементов
    for num in random_data:
        f.write(f"{num} ")

print(f"Сгенерирован файл input_data.txt с {NUM_ELEMENTS} случайными числами")