# Simple DNA Profiler Makefile

build:
	mkdir -p bin
	g++ -std=c++17 -Wall -Wextra -Iheader src/checkSTR.cpp src/dnaprofiler.cpp -o bin/dnaprofiler

test-small: build
	@echo "\nTesting small datasets..."
	@for i in 3 4; do \
		echo "Testing data/small/$$i.txt..."; \
		bin/dnaprofiler -d data/small/database.csv -s data/small/$$i.txt; \
		echo ""; \
	done

test-large: build
	@echo "\nTesting large datasets..."
	@for i in 11 17; do \
		echo "Testing data/large/$$i.txt..."; \
		bin/dnaprofiler -d data/large/database.csv -s data/large/$$i.txt; \
		echo ""; \
	done

clean:
	rm -rf bin