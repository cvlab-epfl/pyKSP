/************************************************************************/
/* (c) 2009-2011 Ecole Polytechnique Federale de Lausanne               */
/* All rights reserved.                                                 */
/*                                                                      */
/* EPFL grants a non-exclusive and non-transferable license for non     */
/* commercial use of the Software for education and research purposes   */
/* only. Any other use of the Software is expressly excluded.           */
/*                                                                      */
/* Redistribution of the Software in source and binary forms, with or   */
/* without modification, is not permitted.                              */
/*                                                                      */
/* Written by Engin Turetken.                                           */
/*                                                                      */
/* http://cvlab.epfl.ch/research/body/surv                              */
/* Contact <pom@epfl.ch> for comments & bug reports.                    */
/************************************************************************/

#ifndef K_SHORTHEST_PATH_COMPUTER_H
#define K_SHORTHEST_PATH_COMPUTER_H

#include "ksp_graph.h"

class KShorthestPathComputer
{
	
public:

  static int test( int a);
  // Returns the number of paths found before one of the stopping criteria is met
  static int ComputeKShorthestNodeDisjointPaths(
                                                KShorthestPathGraph& PathGraph,
                                                int nMaxNoOfPaths, float fMaxPathLength,	// Stopping criteria
                                                unsigned char* pLabelledNodeObjectMap,	// Path label map on the nodes of the original graph
                                                bool bTransformEdgeLengths = true);
 private:
	
  static void TransformShortestPath( 
                                    KShorthestPathGraph& PathGraph,
                                    int nSourceNode, int nTerminalNode,
                                    std::vector< int >& rPreds, 
                                    int* pnReversedNodeMap,  
                                    int* pnKShorthestPathPreds,
                                    int* pnTerminalNodePreds,
                                    int nCurrentPathIndex);
									
  static void TransformShortestPath_G( 
                                    KShorthestPathGraph::BaseGraphType & G,
                                    int nSourceNode, int nTerminalNode,
                                    std::vector< int >& rPreds, 
                                    int* pnReversedNodeMap,  
                                    int* pnKShorthestPathPreds,
                                    int* pnTerminalNodePreds,
                                    int nCurrentPathIndex);
	
};

#endif
