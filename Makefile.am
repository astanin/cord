# AUTOMAKE_OPTIONS = foreign
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I./iniparser/src

SUBDIRS = iniparser .

# iniParser
cord_LDADD = -L./iniparser -liniparser

confdir = $(datadir)/$(PACKAGE)
docdir = $(datadir)/doc/$(PACKAGE)

bin_PROGRAMS = cord
bin_SCRIPTS = scripts/cord_h5togp scripts/cord_h5todx scripts/cord_gp2pdf scripts/cord_pdf2jpg scripts/cord_plots.rb scripts/cord_gp2eps scripts/cord_eps2jpg
cord_SOURCES = cord.cc dmesh.cc dmeshops.cc pradi.cc solver.cc tumeval.cc growth.cc nutrient.cc spmat.cc utils.cc rdsolve.cc slesolve.cc
conf_DATA = cord.ini
doc_DATA = README NEWS COPYING AUTHORS

# headers
cord_SOURCES += amesh2d.h boundary.h dmesh.h dmeshops.h function.h gfm.h global.h growth.h meshenum.h newton.h nutrient.h params.h pradi.h solver.h spmat.h tumeval.h utils.h rdsolve.h slesolve.h singleton.h

# LSolver
cord_SOURCES += lsolver/bicgsq.h lsolver/bicgstab.h lsolver/cblas.h lsolver/cghs.h lsolver/gmres.h lsolver/matrix.h

EXTRA_DIST = cord.ini gsl test scripts 
EXTRA_DIST += iniparser/AUTHORS iniparser/LICENSE iniparser/README iniparser/Makefile.in iniparser/INSTALL iniparser/Makefile.am iniparser/Makefile.original iniparser/src
EXTRA_DIST += lsolver/README lsolver/www/index.html lsolver/www/line_eyes.gif lsolver/www/owl2.gif

dist-hook:
	rm -rf `find $(distdir) -name .svn`


