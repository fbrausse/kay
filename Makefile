
DESTDIR ?= /usr/local
INSTALL  = install

HEADERS = \
	kay/bits.hh \
	kay/gmpxx.hh \
	kay/flintxx.hh \
	kay/numbers.hh \

.PHONY: install uninstall clean

install:
	mkdir -p $(DESTDIR)/include/kay && \
	$(INSTALL) -m 0644 $(addprefix include/,$(HEADERS)) $(DESTDIR)/include/kay/

uninstall:
	$(RM) $(addprefix $(DESTDIR)/include/,$(HEADERS))

clean:
