# Get number of logical cores
num_cores=$(nproc)
echo "Logical cores: $num_cores"
iter=1000

# Get processor name (remove spaces and slashes for filename safety)
cpu_name=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^[ \t]*//' | tr ' /' '__')
results_file="results_${cpu_name}_${num_cores}cores.csv"

# Create results.csv header if file does not exist
header="datetime,test_name,exec_time,seed,all_equal"
for ((i=0; i<num_cores; i++)); do
    header+=",thread_$i"
done
if [ ! -f "$results_file" ]; then
    echo "$header" > "$results_file"
fi

# List of tests (add/remove as needed)
# tests=("USE_INT64")
# tests=("USE_FLOAT")
# tests=("USE_DOUBLE")
# tests=("USE_MATMUL")
# tests=("USE_ROOT")
tests=("USE_INT64" "USE_FLOAT" "USE_DOUBLE" "USE_MATMUL" "USE_ROOT")

# Compile and run for each test
for test in "${tests[@]}"; do
    exe="main_${test,,}" # e.g. main_use_float
    echo "Compiling main.cpp with -D$test -> $exe"
    g++ -Werror -Wall -Wextra -pedantic -ffinite-math-only -fno-signed-zeros --std=c++17 -O0 main.cpp -o "$exe" -lpthread -D$test
    if [ $? -ne 0 ]; then
        echo "Compilation failed for $test"
        continue
    fi

    for ((run=1; run<=iter; run++)); do
        datetime=$(date '+%Y/%m/%d %H:%M:%S')
        echo "[$datetime] Running test: $exe (try $run/$iter)"
        "./$exe" "$num_cores" "$results_file"
        status=$?
        if [ $status -ne 0 ]; then
            echo "$exe exited with failure (exit code $status)."
            break
        fi
    done
done