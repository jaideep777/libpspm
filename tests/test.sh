FILE=species_init
g++ -I../include -c ../src/species.cpp $FILE.cpp
g++ -o 1 species.o $FILE.o
./1 && \
	printf "%b" "\033[0;32m[PASS]\033[m" ": $FILE \n"  || \
	printf "%b" "\033[1;31m[FAIL]\033[m" ": $FILE \n"
