from libcpp.vector cimport vector 
from libcpp cimport bool
ctypedef vector[int] int_vec
ctypedef vector[float] float_vec


cdef extern from "ksp_computer.h": 
    cdef cppclass KShorthestPathComputer:
        int ComputeKShorthestNodeDisjointPaths(
                                                KShorthestPathGraph& PathGraph,
                                                int nMaxNoOfPaths, float fMaxPathLength,
                                                unsigned char* pLabelledNodeObjectMap,
                                                bool bTransformEdgeLengths = true);

cdef extern from "ksp_graph.h": 
    cdef cppclass KShorthestPathGraph:
        KShorthestPathGraph( 
                      float_vec pfData,
                      int nDataWidth,
                      int nDataHeight,
                      int nDataDepth,
                      int nNodeNeighborhoodSize,
                      int_vec pnSrcAndDstNeighIndices)
        int_vec getPath(int first_frame,int last_frame)
        int_vec getTest(int first_frame,int last_frame)

		
cdef class pyKShorthestPathGraph: 
    cdef KShorthestPathGraph* thisptr # hold a C++ instance
    def __cinit__( self,float_vec pfData,
                      int nDataWidth,
                      int nDataHeight,
                      int nDataDepth,
                      int nNodeNeighborhoodSize,
                      int_vec pnSrcAndDstNeighIndices):
        self.thisptr = new KShorthestPathGraph( pfData,
                      nDataWidth,
                      nDataHeight,
                      nDataDepth,
                      nNodeNeighborhoodSize,
                      pnSrcAndDstNeighIndices)
    def __dealloc__(self):
        del self.thisptr
    def getPath(self,int first_frame,int last_frame):
        return self.thisptr.getPath(first_frame,last_frame)
    def getTest(self,int first_frame,int last_frame):
        return self.thisptr.getTest(first_frame,last_frame)


cdef class pyKShorthestPathComputer: 
    cdef KShorthestPathComputer* thisptr # hold a C++ instance
    def __cinit__( self):
        self.thisptr = new KShorthestPathComputer()

    def __dealloc__(self):
        del self.thisptr