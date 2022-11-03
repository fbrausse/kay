
DESTDIR ?= /usr/local
INSTALL  = install

HEADERS = \
	kay/bits.hh \
	kay/compiletime.hh \
	kay/gmpxx.hh \
	kay/flintxx.hh \
	kay/numbers.hh \
	kay/numbits.hh \
	kay/dbl-ival.hh \

.PHONY: install uninstall clean

$(DESTDIR)/%/:
	mkdir -p $@

$(DESTDIR)/include/%: include/% | $(DESTDIR)/include/kay/
	$(INSTALL) -m 0644 $< $@

install: $(addprefix $(DESTDIR)/include/,$(HEADERS))

uninstall:
	$(RM) $(addprefix $(DESTDIR)/include/,$(HEADERS))

clean:
