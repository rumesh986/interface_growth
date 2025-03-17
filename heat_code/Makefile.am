# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= temps2 temps3 temps4 temps_poisson temps_single_unsteady double_unsteady moving_double_unsteady

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
temps2_SOURCES = temps2.cc
temps3_SOURCES = temps3.cc
temps4_SOURCES = temps4.cc
temps_poisson_SOURCES = temps_poisson.cc
temps_single_unsteady_SOURCES = temps_single_unsteady.cc
double_unsteady_SOURCES = double_unsteady.cc
moving_double_unsteady_SOURCES = moving_double_unsteady.cc

# Required libraries: Only the "generic" and "poisson" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. The generic sources also require the "external" libraries
# that are shipped with oomph-lib. The Fortran libraries, $(FLIBS), get 
# included just in case we decide to use a solver that involves Fortran 
# sources. 
# NOTE: The order in which libraries are specified is important! 
#       The most specific ones (that might involve references to other
#       libraries) need to be listed first, the more basic ones 
#       later. In this example we have (from right to left, i.e. from
#       general to specific):
#       -- The fortran libraries: They  are compiler specific and
#          obviously can't depend on any code that we (or others)
#          have written. $(FLIBS) is a variable that automake will translate
#          to the actual fortran libraries.
#       -- The external (third party) libraries: They cannot depend on 
#          any of our code. The variable $(EXTERNAL_LIBS) is defined
#          in the machine-generated file configure.ac in oomph-lib's 
#          home directory. 
#       -- Oomph-lib's generic library contains oomph-lib's fundamental 
#          objects which do not depend on any specific system of PDEs 
#          or element types.
#       -- Finally, oomph-lib's poisson library contains oomph-lib's 
#          Poisson elements which refer to objects from the generic
#          library. 
temps2_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
temps3_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
temps4_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
temps_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
temps_single_unsteady_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
double_unsteady_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
moving_double_unsteady_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  
