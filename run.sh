#!/bin/bash

set -xe

gcc main.c -o main -std=c11 -Wall -Wextra -Wno-strict-aliasing -Wno-unknown-pragmas -lm

./main

rm main