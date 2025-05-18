#!/bin/bash

# Конфигурация
EPS="1e-6"
THREADS_LIST=(1 2 4 8 12)
EXECUTABLE="./a.out"
LOGFILE="integration_benchmark.log"

# Очистка лог-файла
echo "Benchmark started at $(date)" > "$LOGFILE"
echo "EPSILON = $EPS" >> "$LOGFILE"
echo "Threads | Time (s) | Speedup | Efficiency" >> "$LOGFILE"
echo "--------|----------|---------|-----------" >> "$LOGFILE"

# Временное хранилище для T1
T1=0

for THREADS in "${THREADS_LIST[@]}"; do
    echo -n "Running with $THREADS thread(s)... "

    # Замер времени выполнения
    START=$(date +%s.%N)
    $EXECUTABLE "$THREADS" "$EPS" "0.02" "6.0" > /dev/null
    END=$(date +%s.%N)

    # Расчёт времени
    TIME=$(echo "$END - $START" | bc)

    # Для одного потока — сохраняем базовое время
    if [ "$THREADS" -eq 1 ]; then
        T1=$TIME
        SPEEDUP=1
        EFFICIENCY=1
    else
        SPEEDUP=$(echo "scale=4; $T1 / $TIME" | bc)
        EFFICIENCY=$(echo "scale=4; $SPEEDUP / $THREADS" | bc)
    fi

    # Запись в лог
    printf "%7d | %8.4f | %7.4f | %9.4f\n" "$THREADS" "$TIME" "$SPEEDUP" "$EFFICIENCY" >> "$LOGFILE"

    echo "done (time: $TIME s)"
done

echo "Results saved to $LOGFILE"
