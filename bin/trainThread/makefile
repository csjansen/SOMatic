CC = g++
FILES = main.cpp
OUT_EXE = trainsom

build: $(FILES)
		$(CC) $(FILES) -O3 -fpermissive -std=c++11 -pthread -o $(OUT_EXE)

clean:
	rm -f trainsom
