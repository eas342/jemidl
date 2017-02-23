all:
	make mult2
	make fzroots2
	make drizzle1d
	make clean
	make jem_avsigclip

mult2: mult2.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c mult2.c
	ld mult2.so.1 -o mult2.so.1.0 -lc mult2.o

fzroots2: fzroots2.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c fzroots2.c
	ld fzroots2.so.1 -o fzroots2.so.1.0 -lc fzroots2.o

drizzle1d: drizzle1d.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c drizzle1d.c
	ld drizzle1d.so.1 -o drizzle1d.so.1.0 -lc drizzle1d.o
	cp drizzle1d.so.1.0 drizzle1d.so

clean: clean.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c clean.c
	ld clean.so.1 -o clean.so.1.0 -lc clean.o
	cp clean.so.1.0 clean.so

clean2: clean2.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c clean2.c
	ld clean2.so.1 -o clean2.so.1.0 -lc clean2.o
	cp clean2.so.1.0 clean2.so

clean3: clean3.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c clean3.c
	ld clean3.so.1 -o clean3.so.1.0 -lc clean3.o
	cp clean3.so.1.0 clean3.so

jem_avsigclip: jem_avsigclip.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c jem_avsigclip.c
	ld jem_avsigclip.so.1 -o jem_avsigclip.so.1.0 -lc jem_avsigclip.o
	cp jem_avsigclip.so.1.0 jem_avsigclip.so


keyword_demo: keyword_demo.c
	gcc -I/Applications/itt/idl/idl81/external/include -fPIC -c keyword_demo.c
	ld keyword_demo.so.1 -o keyword_demo.so.1.0 -lc keyword_demo.o
