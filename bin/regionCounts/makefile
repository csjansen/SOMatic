CC = g++
FILES = main.cpp
OUT_EXE = regionCounts

build: $(FILES)
		$(CC) $(FILES) -O3 -fpermissive -o $(OUT_EXE) -std=c++11 -pthread

clean:
	rm -f regionCounts
