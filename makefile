all: sp sfc srr sbc cmp

sp: 
	g++ -O4 -o sp4 src/splitPairs.cpp -std=c++11

sfc: 
	gcc -O4 -o sfc src/split_columns.c

srr: 
	gcc -O4 -o srr src/split_read_rsw.c

sbc: 
	gcc -O4 -o sbc src/split_on_chrom.c


cmp:
	g++ -O4 -o rsr_compare src/compare.cpp -std=c++11
	

clean:
	rm -f sp4 sfc srr sbc rsr_compare


