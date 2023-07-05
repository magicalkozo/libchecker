#!/bin/bash

set -xe

gcc main.c -o main -std=c11 -Wall -Wextra -lm -Wno-strict-aliasing -Wno-unknown-pragmas
# g++ main.cpp -o main -std=c++17 -Wall -Wextra -I .

./main

rm main