########################################################################
# (c) 2009-2011 Ecole Polytechnique Federale de Lausanne               #
# All rights reserved.                                                 #
#                                                                      #
# EPFL grants a non-exclusive and non-transferable license for non     #
# commercial use of the Software for education and research purposes   #
# only. Any other use of the Software is expressly excluded.           #
#                                                                      #
# Redistribution of the Software in source and binary forms, with or   #
# without modification, is not permitted.                              #
#                                                                      #
# Written by Engin Turetken and Jerome Berclaz.                        #
#                                                                      #
# http://cvlab.epfl.ch/research/body/surv                              #
# Contact <pom@epfl.ch> for comments & bug reports.                    #
########################################################################

LDFLAGS=
CXXFLAGS=-O3 -msse2 -mfpmath=sse

all: ksp 

ksp: ksp_graph.o ksp_computer.o main.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) *.o *~ ksp
