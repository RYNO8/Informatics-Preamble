#!/bin/bash

for FILE in *.test.cpp;
do g++ $FILE -o temp && ./temp;
done