
FILES = cluster cutree diffmap genes GO map partition regionCounts score trainQuick trainThread units metasomgene metasomThread LinkUnits LinkMeta LinkAnalyze 

build:
	for dir in $(FILES); do \
		$(MAKE) -C $$dir; \
	done
clean:
	for dir in $(FILES); do \
                $(MAKE) -C $$dir clean; \
        done

