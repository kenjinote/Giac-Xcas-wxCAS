SUBDIRS = 	src check intl po debianold debian doc examples 
DIST_SUBDIRS = 	src check intl po debianold debian doc examples 
EXTRA_DIST = path.h.in setup inst TROUBLES Init.scm compte INSTALACION LEAME.WIN LEAME.ipaq giac.specs giac.specs32 fltk.1.1.changes README.MACOSX pariinl.h changeversion mkosx mkdebian mkdebian2 win2unix.cc addu.cc xcas.nsi xcas.desktop xcas.xpm xcas.applications README.visualcpp config.h.win mkiphone build_for_iphoneos path.iphone simupath.iphone config.h.iphone compile.mingw

# Rule to build tar-bzipped distribution package
$(PACKAGE)-$(VERSION).tar.gz: dist

## FIXME: installing autoconf-generated config.h is just PLAIN WRONG.

install-exec-hook:
	cp config.h $(includedir)/giac
	./mkinstalldirs $(prefix)/share/
	./mkinstalldirs $(prefix)/share/pixmaps
	cp xcas.xpm $(prefix)/share/pixmaps/xcas.xpm
	./mkinstalldirs $(prefix)/share/applications
	cp xcas.desktop $(prefix)/share/applications/xcas.desktop
	./mkinstalldirs $(prefix)/share/application-registry
	cp xcas.applications $(prefix)/share/application-registry/xcas.applications

ACLOCAL_AMFLAGS = -I m4

