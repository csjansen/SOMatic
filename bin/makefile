
FILES = cluster cutree diffmap genes GO map partition regionCounts score train trainQuick trainThread units 

build:
	for dir in $(FILES); do \
		$(MAKE) -C $$dir; \
	done
clean:
	for dir in $(FILES); do \
                $(MAKE) -C $$dir clean; \
        done

