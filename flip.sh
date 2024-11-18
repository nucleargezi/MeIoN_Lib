g++ -std=c++20 -DMeIoN test.cpp -o flip_test
g++ -std=c++20 -DMeIoN tt.cpp -o flip_std
g++ -std=c++20 -DMeIoN gen.cpp -o flip_gen

cnt = 1

while true; do
    ./flip_gen > flip_data
    ./flip_test < flip_data > flip_my_out
    ./flip_std < flip_data > flip_std_out
    echo "times: $cnt"
    if ! diff -b -B flip_my_out flip_std_out > /dev/null; then
        echo "wrong answer"
        break
    else 
        echo "accept"
    fi
    ((cnt++))
done