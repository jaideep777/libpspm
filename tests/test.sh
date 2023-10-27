g++ -g -pg -o 1 species.o rkck45.o lsoda.o solver.o ebt.o iebt.o $FILE.o
g++ -g -pg -O3 -I../include -c ../src/species.cpp ../src/rkck45.cpp ../src/lsoda.cpp ../src/solver.cpp ../src/ebt.cpp ../src/iebt.cpp  $FILE.cpp
FILE=iebt_equil
./1 && \
	printf "%b" "\033[0;32m[PASS]\033[m" ": $FILE \n"  || \
	printf "%b" "\033[1;31m[FAIL]\033[m" ": $FILE \n"
