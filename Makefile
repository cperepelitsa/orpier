PREFIX=/usr/local
SYSCONFDIR=$(PREFIX)/etc
MANDIR=$(PREFIX)/man

orpie: install.ml
	ocamlbuild -use-ocamlfind -pp camlp4o.opt main.native
	mv main.native orpie

DOCS=doc/orpie.1 doc/orpierc.5 doc/orpie-curses-keys.1
$(DOCS):
	make -C doc orpie.1 orpierc.5 orpie-curses-keys.1

install.ml:
	sed -n \
		-e "1s!@prefix@!$(PREFIX)!p" \
		-e "2s!@sysconfdir@!$(SYSCONFDIR)!p" \
		install.ml.in > install.ml

.PHONY: install-bin install-doc install
install-bin: orpie
	install -d "$(PREFIX)/bin" "$(SYSCONFDIR)"
	install orpie "$(PREFIX)/bin"
	install -m 644 orpierc "$(SYSCONFDIR)"

install-doc: $(DOCS)
	install -d "$(MANDIR)/man1" "$(MANDIR)/man5"
	install -m 644 doc/orpie.1 "$(MANDIR)/man1"
	install -m 644 doc/orpie-curses-keys.1 "$(MANDIR)/man1"
	install -m 644 doc/orpierc.5 "$(MANDIR)/man5"

install: install-bin install-doc

clean:
	ocamlbuild -clean
	rm -f install.ml
