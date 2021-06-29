#!/bin/bash
g++ -o Bin/gato.o -c gato.cpp
 
g++ -o Bin/sed.o -c sed.cpp
g++ -o Bin/sed Bin/gato.o Bin/sed.o 
