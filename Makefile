PREFIX  = $(HOME)

# Dependence
# ==========

TOTAL=$(HOME)
NEUS=$(HOME)
MAD=$(HOME)
UNIC=$(HOME)


# Define CXX, CXXFLAGS, SOFLAGS & LIBS
# ====================================

OPT2       = -O2

ROOTCONFIG = root-config

ARCH      := $(shell $(ROOTCONFIG) --arch)
ALTCXX    := $(shell $(ROOTCONFIG) --cxx)
ROOTCFLAGS:= $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS  := $(shell $(ROOTCONFIG) --libs)

ARCHOK     = no

ifeq ($(ARCH),linux)
  CXX      = g++
  CXXFLAGS = $(OPT2) -Wall -fPIC
  SOFLAGS  = -shared
  ARCHOK   = yes
endif

ifeq ($(ARCH),linuxx8664gcc)
  CXX      = g++
  CXXFLAGS = $(OPT2) -Wall -fPIC
  SOFLAGS  = -shared
  ARCHOK   = yes
endif

ifeq ($(ARCH),win32gcc)
  CXX      = g++
  CXXFLAGS = $(OPT2) -Wall
  SOFLAGS  = -shared
  ARCHOK   = yes
endif

# If arch != the above options, a error message is given
ifeq ($(ARCHOK),no)
  $(error $(ARCH) invalid architecture)
endif

# In case that ALTCXX = 0, g++ will be used as CXX
ifneq ($(ALTCXX),)
  CXX = $(ALTCXX)
endif

# Finally, define CXXFLAGS & LIBS
CXXFLAGS+= $(ROOTCFLAGS)
CXXFLAGS+= -I$(UNIC)/include -I$(NEUS)/include -I$(MAD)/include -g
LIBS     = $(ROOTLIBS) -L$(NEUS)/lib -lNEUS -lTOTAL -L$(MAD) -lMAD


# Define things related to rootcling
# =================================

LIBNAME := CNNS
ROOTIFIED_SOURCE := $(LIBNAME)Dict.cc
ROOTIFIED_HEADER := $(ROOTIFIED_SOURCE:.cc=.h)
ROOTIFIED_OBJECT := $(ROOTIFIED_SOURCE:.cc=.o)


# Define SOURCES, HEADERS & OBJECTS 
# ==========================================

SRCS = $(wildcard *.C)
EXES = $(SRCS:.C=.exe)

SOURCES = $(filter-out $(ROOTIFIED_SOURCE), $(wildcard *.cc))
HEADERS = $(filter-out $(ROOTIFIED_HEADER) $(LINKDEF), $(wildcard *.h))
OBJECTS = $(SOURCES:.cc=.o)
DEPFILE = $(SOURCES:.cc=.d)

LINKDEF = LinkDef.h


# Define LIBRARY, ROOTMAP & variables to create them
# ==================================================

LIBRARY = lib$(LIBNAME).so
ROOTMAP = $(LIBRARY:.so=.rootmap)

SYMBOLS = `nm -CPu $(LIBRARY) |\
	  awk -F: '/^T/{printf("(%s)\n",$$1)}' |\
	  sort -u | tr '\n' '|'`
ALLLIBS = $(wildcard $(shell $(ROOTCONFIG) --libdir)/*.so)
DEPENDS = $(shell symbols=$(SYMBOLS); \
	  for so in $(ALLLIBS); \
	  do nm -CD --defined-only $$so |\
	  grep -E 'T ('$$symbols')::' > /dev/null &&\
	  echo $$so;\
	  done | sort -u | tr '\n' ' ')

# Action starts
# =============

# the first target is the default target, it depends on $(ROOTMAP)
# before "make all", make will include all other makefiles specified
# by the include command
all: $(EXES)
	@echo
	@echo "* Done!"
	@echo 

# include *.d files, which are makefiles defining dependencies between files
ifeq ($(filter info clean,$(MAKECMDGOALS)),)
  -include $(DEPFILE)
endif

# rules to create *.d files
%.d:%.cc
	@echo creating $@
	@set -e; rm -f $@; \
	  $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$ 

# lib$(LIBNAME).so depends on all *.o files.
#  The flag "-shared" is used to create shared libs
#  $@ represents the target, that is, lib$(LIBNAME).so
#  $^ represents all the prerequisites, i.e., all *.o files
$(LIBRARY): $(ROOTIFIED_OBJECT) $(OBJECTS)
	@echo
	@echo "* Creating shared library:"
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) -o $@ $^

# An xxx.o file depends on xxx.cc. It is created with the command:
# 	g++ -c xxx.cc -o xxx.o
# Since this is obvious, "make" automatically does it. 
# There is no need to explicitly write down the rules to do it.

# To use classes & functions directly in ROOT. One has to do the following
$(ROOTIFIED_SOURCE): $(HEADERS) $(LINKDEF)
	@echo 
	@echo "* Rootifying files:" 
	@rm -f $(ROOTIFIED_SOURCE) $(ROOTIFIED_HEADER) 
	rootcling -f $@ -cxxflags="$(CXXFLAGS)" -s lib$(LIBNAME) -rml lib$(LIBNAME) -rmf $(ROOTMAP) $^
	@echo 
	@echo "* Creating object files:" 

info: 
	@echo
	@echo "target:   $(LIBRARY) $(ROOTMAP)"
	@echo "sources:  $(SOURCES)"
	@echo "headers:  $(HEADERS)"
	@echo "objects:  $(OBJECTS)"
	@echo 
	@echo "compiler: $(CXX)"
	@echo "flags:    $(CXXFLAGS)"
	@echo "libs:     $(LIBS)"
	@echo

clean:
	$(RM) *.exe *.o *.d *.d.* *~ *Dict* $(ROOTMAP) $(LIBRARY)

tags:
	ctags --c-kinds=+p $(HEADERS) $(SOURCES)

install: $(LIBRARY)
	@echo
	@echo "* Installing library to PREFIX=$(PREFIX)"
	@echo -n "checking if $(PREFIX) exists..."
	@if [ -d $(PREFIX) ]; then \
	  echo "yes."; \
	else \
	  echo "no."; \
	  echo "mkdir $(PREFIX)..."; \
	  mkdir $(PREFIX); \
	fi
	@echo -n "copying lib$(LIBNAME).so to $(PREFIX)/lib..."
	@if [ -d $(PREFIX)/lib ]; then \
	  cp lib$(LIBNAME).so $(PREFIX)/lib; \
	  cp lib$(LIBNAME).rootmap $(PREFIX)/lib; \
	else \
	  mkdir $(PREFIX)/lib; \
	  cp lib$(LIBNAME).so $(PREFIX)/lib; \
	  cp lib$(LIBNAME).rootmap $(PREFIX)/lib; \
	fi; 
	@echo "done."; 
	@echo -n "copying *.h to $(PREFIX)/include/$(LIBNAME)..."
	@if [ -d $(PREFIX)/include ]; then \
	  if [ -d $(PREFIX)/include/$(LIBNAME) ]; then \
	    cp $(HEADERS) $(PREFIX)/include/$(LIBNAME); \
	  else \
	    mkdir $(PREFIX)/include/$(LIBNAME); \
	    cp $(HEADERS) $(PREFIX)/include/$(LIBNAME); \
	  fi; \
	else \
	  mkdir $(PREFIX)/include; \
	  mkdir -p $(PREFIX)/include/$(LIBNAME); \
	  cp *.h $(PREFIX)/include/$(LIBNAME); \
	fi
	@echo "done."
	@echo
	@echo "* Creating executables:"

uninstall:
	$(RM) -r $(PREFIX)/include/$(LIBNAME)
	$(RM) -r $(PREFIX)/lib/lib$(LIBNAME).so


$(EXES):%.exe:%.C install
	$(CXX) $< $(CXXFLAGS) -L. -l$(LIBNAME) $(LIBS) -o $@

.PHONY: all info tags clean
