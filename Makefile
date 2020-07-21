
DESTDIR ?= /usr/local
INSTALL  = install

HEADERS = \
	kay/bits.hh \
	kay/compiletime.hh \
	kay/gmpxx.hh \
	kay/flintxx.hh \
	kay/numbers.hh \
	kay/numbits.hh \

.PHONY: install uninstall clean

$(DESTDIR)/include/%: include/% | $(DESTDIR)/include/
	$(INSTALL) -m 0644 $< $@

install: $(addprefix $(DESTDIR)/include/,$(HEADERS))

$(DESTDIR)/%/:
	mkdir -p $@

uninstall:
	$(RM) $(addprefix $(DESTDIR)/include/,$(HEADERS))

clean:
