# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS=erf2d1 erf2d2 erf2d3 erf2d4 erf2d5 erf2d6 erf2d7 erf2d3p1 erf2d3p2 erf2d8 erf2d9 erf2d10 erf2d11 erf2d12 erf2d13 erf2d14 erf2d15 erf2d16 erf2d17 erf2d18

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
erf2d1_SOURCES = erf2d1.cc erf2d1.h solver.h
erf2d2_SOURCES = erf2d2.cc erf2d2.h solver.h two_layer_mesh.h
erf2d3_SOURCES = erf2d3.cc erf2d3.h solver.h two_layer_mesh.h
erf2d4_SOURCES = erf2d4.cc erf2d4.h solver.h two_layer_mesh.h
erf2d5_SOURCES = erf2d5.cc erf2d5.h solver.h two_layer_mesh.h
erf2d6_SOURCES = erf2d6.cc two_layer_mesh.h
erf2d7_SOURCES = erf2d7.cc two_layer_mesh.h
erf2d3p1_SOURCES = erf2d3p1.cc two_layer_mesh.h three_layer_mesh.h
erf2d3p2_SOURCES = erf2d3p2.cc two_layer_mesh.h three_layer_mesh.h
erf2d8_SOURCES = erf2d8.cc two_layer_mesh.h
erf2d9_SOURCES = erf2d9.cc two_phase_free_boundary_mesh.h
erf2d10_SOURCES = erf2d10.cc two_phase_free_boundary_mesh.h
erf2d11_SOURCES = erf2d11.cc two_phase_free_boundary_mesh.h
erf2d12_SOURCES = erf2d12.cc two_phase_free_boundary_mesh.h
erf2d13_SOURCES = erf2d13.cc two_phase_free_boundary_mesh.h
erf2d14_SOURCES = erf2d14.cc two_phase_free_boundary_mesh.h
erf2d15_SOURCES = erf2d15.cc two_phase_free_boundary_mesh.h
erf2d16_SOURCES = erf2d16.cc two_phase_free_boundary_mesh.h
erf2d17_SOURCES = erf2d17.cc two_phase_free_boundary_mesh.h
erf2d18_SOURCES = erf2d18.cc two_phase_free_boundary_mesh.h

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
erf2d1_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d2_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d3_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d4_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d5_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d6_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d7_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d3p1_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d3p2_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d8_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d9_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d10_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d11_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d12_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d13_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d14_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d15_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d16_LDADD = -L@libdir@ -lsolid -lconstitutive -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d17_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
erf2d18_LDADD = -L@libdir@ -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@  
