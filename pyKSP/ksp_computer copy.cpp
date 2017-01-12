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

#include "ksp_computer.h"

#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/dag_shortest_paths.hpp>


int KShorthestPathComputer::test( int a){
	return a+2;
}

int KShorthestPathComputer::ComputeKShorthestNodeDisjointPaths(
                                                                KShorthestPathGraph& PathGraph,
                                                               int nMaxNoOfPaths, float fMaxPathLength,
                                                               unsigned char* pLabelledNodeObjectMap, 
                                                               bool bTransformEdgeLengths)
{
  // Declarations
  //KShorthestPathGraph PathGraph = PathGraph_or;
  KShorthestPathGraph::BaseGraphType & G = PathGraph.GetBaseGraph();
  typedef boost::graph_traits < KShorthestPathGraph::BaseGraphType >::vertex_descriptor VertexDesc;
  typedef boost::graph_traits < KShorthestPathGraph::BaseGraphType >::edge_descriptor EdgeDesc;
  typedef std::pair<int, int> Edge;
  int nNoOfOrigNodes;
  int nNoOfCurrNodes;
  float* pfCostsInOrigGraph;
  int* pnReversedNodeMap;
  int* pnKShorthestPathPreds;
  int* pnTerminalNodePreds;
  int nPathCounter;
  int nNoOfPaths;
  int nBufferNode;
  int nSourceNode;
  int nTerminalNode; 
  boost::property_map <	KShorthestPathGraph::BaseGraphType, 
    boost::edge_weight_t >::type EdgeWeightMap;
  boost::graph_traits < KShorthestPathGraph::BaseGraphType >::edge_iterator Ei, Eend;
  std::vector< int >::iterator PredsIter;
  std::vector< int >::iterator PredsEndIter;
	
  //Initializations
  pfCostsInOrigGraph = new float[nMaxNoOfPaths];
  nNoOfOrigNodes = boost::num_vertices(G);
  pnReversedNodeMap = (int*)calloc(2 * nNoOfOrigNodes, sizeof(int));
  pnKShorthestPathPreds = new int[nNoOfOrigNodes];
  pnTerminalNodePreds = new int[nMaxNoOfPaths];
  nSourceNode = PathGraph.GetSrcNodeIndx();
  nTerminalNode = PathGraph.GetDstNodeIndx();
	
  // Finding the shorthest path (first) on the inputted graph possibly 
  // having negative edges, but no negatice cycles. The graph is fully
  // explored.
  std::vector< int > Preds(2 * nNoOfOrigNodes);
  std::vector< float > Dists(2 * nNoOfOrigNodes);
	
  boost::dag_shortest_paths(G, nSourceNode, 
                            boost::predecessor_map(
                                                   &Preds[0]).distance_map(&Dists[0]));
	
  for( nPathCounter = 0; nPathCounter < nMaxNoOfPaths; nPathCounter++ )
    {
      // Adding the cost of the shorthest path to the list.
      if( nPathCounter != 0)
        {
          if( bTransformEdgeLengths )
            {
              pfCostsInOrigGraph[nPathCounter] = 
                pfCostsInOrigGraph[nPathCounter-1] + Dists[nTerminalNode];
            }
          else
            {
              pfCostsInOrigGraph[nPathCounter] = Dists[nTerminalNode];
            }
        }
      else
        {
          pfCostsInOrigGraph[nPathCounter] = Dists[nTerminalNode];
        }
		
		
      // Checking if the max path length stopping criteria is met.
      if( pfCostsInOrigGraph[nPathCounter] >= fMaxPathLength )
        {
          break;
        }
		
      // Transform edge lengths (i.e., compute the canonical graph) so that 
      // no negative cost edges appear in the transformed graph.
      if( bTransformEdgeLengths )
        {			
          EdgeWeightMap = get(boost::edge_weight, G);
			
          for( boost::tie(Ei, Eend) = boost::edges(G); Ei != Eend; ++Ei) 
            {					
              const EdgeDesc& rBufferEdgeDesc = (*Ei);

              boost::put( EdgeWeightMap, rBufferEdgeDesc, 
                          boost::get( EdgeWeightMap, rBufferEdgeDesc) + Dists[rBufferEdgeDesc.m_source] - 
                          Dists[rBufferEdgeDesc.m_target]);
            }
        }
		
      // Transforming the last shorthest path (i.e., splitting the path nodes 
      // and reversing the path edges).
      TransformShortestPath( PathGraph, nSourceNode, nTerminalNode, Preds, 
                             pnReversedNodeMap, pnKShorthestPathPreds, 
                             pnTerminalNodePreds, nPathCounter);				
		
      // Finding the next shorthest path on the transformed graph with no 
      // negative edges.
      if( nPathCounter != (nMaxNoOfPaths-1) )
        {
          if( bTransformEdgeLengths )
            {
              boost::dijkstra_shortest_paths(
                                             G, nSourceNode, 
                                             boost::predecessor_map(
                                                                    &Preds[0]).distance_map(&Dists[0]));
            }
          else
            {
              // Initializing the pred and dist arrays for bellman ford.
              // This is a bit ridicilious since this initialization step is
              // not required by DAG and Dijkstra versions of the shorthest path
              // algorithms.
              nNoOfCurrNodes =  boost::num_vertices(G);
              for(int i = 0; i < nNoOfCurrNodes; i++)
                {
                  Preds[i] = i;
                }
              Dists.assign(nNoOfCurrNodes, std::numeric_limits<float>::max());
              // Specifying the source vertex
              Dists[nSourceNode] = 0;
				
              boost::bellman_ford_shortest_paths(
                                                 G, nSourceNode, boost::predecessor_map(
                                                                                        &Preds[0]).distance_map(&Dists[0]));
            }
        }
    }
  nNoOfPaths = nPathCounter;
	
  // Constructing the node based object labels
  for( nPathCounter = 1; nPathCounter <= nNoOfPaths; nPathCounter++)
    {
      nBufferNode = pnTerminalNodePreds[nPathCounter - 1];
      pLabelledNodeObjectMap[nBufferNode] = nPathCounter;
				
      while( nBufferNode != nSourceNode )
        {
          nBufferNode = pnKShorthestPathPreds[nBufferNode];
          pLabelledNodeObjectMap[nBufferNode] = nPathCounter;			
        }
    }
	
  //Deallocations
  delete[] pfCostsInOrigGraph;
  free(pnReversedNodeMap);
  delete[] pnKShorthestPathPreds;
  delete[] pnTerminalNodePreds;
	
  return nNoOfPaths;
}



void KShorthestPathComputer::TransformShortestPath( 
                                                   KShorthestPathGraph& PathGraph,
                                                   int nSourceNode, int nTerminalNode,
                                                   std::vector< int >& rPreds, 
                                                   int* pnReversedNodeMap,  
                                                   int* pnKShorthestPathPreds,
                                                   int* pnTerminalNodePreds,
                                                   int nCurrentPathIndex)
{
  // Declerations and initializations
  KShorthestPathGraph::BaseGraphType & G = PathGraph.GetBaseGraph();
  boost::graph_traits< KShorthestPathGraph::BaseGraphType >::out_edge_iterator Ei, Eend;
  boost::edge_weight_t EdgeWeightPropType;
  typedef boost::graph_traits< KShorthestPathGraph::BaseGraphType >::edge_descriptor EdgeDesc;
  bool bDummyEdgeCheck;
  EdgeDesc BufferEdge;
  int nNoOfCurrentNodes;
  int nCurrentIndx2 = 0;
  int nSuccIndx = nTerminalNode;
  int nCurrentIndx1 = rPreds[nSuccIndx];
  int nPredIndx = rPreds[nCurrentIndx1];
	
	
  if( pnReversedNodeMap[nCurrentIndx1] <= 0 )
    {		
      pnTerminalNodePreds[nCurrentPathIndex] = nCurrentIndx1;
    }
  else
    {
      pnTerminalNodePreds[nCurrentPathIndex] = pnKShorthestPathPreds[nPredIndx];
    }
	
  while( nCurrentIndx1 != nSourceNode )
    {	
		
      if( pnReversedNodeMap[nCurrentIndx1] <= 0 )
        {
          // If this node is not included in any one of the previously 
          // computed shorthest paths (Pm), then split this node and
          // reverse its successor edge in the path.
          pnKShorthestPathPreds[nSuccIndx] = nCurrentIndx1;
			
          // If this node had been previously included in any one of the previously 
          // computed shorthest paths, and hence had been splitted but then merged 
          // again (the second vertex is not deallocated and left isolated 
          // in the graph), then do not create a new vertex but use the one 
          // that had been previously allocated.
          if( pnReversedNodeMap[nCurrentIndx1] < 0 )
            {
              nNoOfCurrentNodes = - pnReversedNodeMap[nCurrentIndx1];
            }
          else
            {
              nNoOfCurrentNodes = boost::num_vertices(G) + 1;
              boost::add_vertex(G);
            }
			
          nCurrentIndx2 = nNoOfCurrentNodes;
          pnReversedNodeMap[nCurrentIndx1] = nNoOfCurrentNodes;
          pnReversedNodeMap[nCurrentIndx2] = nNoOfCurrentNodes;
			
          // Splitting the node
			
          // Transforming all outgoing edges of the first node
          for( boost::tie(Ei, Eend) = boost::out_edges(nCurrentIndx1, G); Ei != Eend; ++Ei) 
            {
              boost::add_edge(nCurrentIndx2, (*Ei).m_target, get( EdgeWeightPropType, G, (*Ei)), G);
            }
			
          // Deleting all outgoing edges of the first node
          boost::clear_out_edges(nCurrentIndx1, G);
			
          // Adding an edge between the two nodes
          boost::add_edge(nCurrentIndx2, nCurrentIndx1, 0, G);
			
          // Reversing (both direction and sign) the successor edge 
          // outgoing from this node
          boost::tie(BufferEdge, bDummyEdgeCheck) = boost::edge(nCurrentIndx2, nSuccIndx, G);
          boost::add_edge(nSuccIndx, nCurrentIndx2, 
                          -boost::get( EdgeWeightPropType, G, BufferEdge), G);
          boost::remove_edge(nCurrentIndx2, nSuccIndx, G);
			
          nSuccIndx = nCurrentIndx1;
          nCurrentIndx1 = nPredIndx;
          nPredIndx = rPreds[nCurrentIndx1];
        }
      else
        {	
          // If entered any one of the previous shorthest paths, continue till
          // exiting it.
			
          pnKShorthestPathPreds[nSuccIndx] = pnKShorthestPathPreds[nPredIndx];
			
          // Do not merge the first node and only reverse the successor edge.
          boost::tie(BufferEdge, bDummyEdgeCheck) = boost::edge(nCurrentIndx1, nSuccIndx, G);
          boost::add_edge(nSuccIndx, nCurrentIndx1, 
                          -boost::get( EdgeWeightPropType, G, BufferEdge), G);
          boost::remove_edge(nCurrentIndx1, nSuccIndx, G);
			
          // Also reverse the pred edge.
          boost::tie(BufferEdge, bDummyEdgeCheck) = boost::edge(nPredIndx, nCurrentIndx1, G);
          boost::add_edge(nCurrentIndx1, nPredIndx, 
                          -boost::get( EdgeWeightPropType, G, BufferEdge), G);
          boost::remove_edge(nPredIndx, nCurrentIndx1, G);
			
          nSuccIndx = nCurrentIndx1;
          nCurrentIndx1 = nPredIndx;
          nCurrentIndx2 = rPreds[nCurrentIndx1];
          nPredIndx = rPreds[nCurrentIndx2];
			
          while( pnReversedNodeMap[nCurrentIndx1] == pnReversedNodeMap[nCurrentIndx2] )
            {
              // Deleting the edge between the parts
              boost::remove_edge(nCurrentIndx2, nCurrentIndx1, G);
				
              // Merging the node by copying all outgoing edges from node2
              // and directing (edge reversing step) the incoming edge of 
              // node2 to node1.
              for( boost::tie(Ei, Eend) = boost::out_edges(nCurrentIndx2, G); Ei != Eend; ++Ei) 
                {
                  boost::add_edge(nCurrentIndx1, (*Ei).m_target, boost::get( EdgeWeightPropType, G, (*Ei)), G);
                }
				
              boost::tie(BufferEdge, bDummyEdgeCheck) = boost::edge(nPredIndx, nCurrentIndx2, G);
              boost::add_edge(nCurrentIndx1, nPredIndx, 
                              -boost::get( EdgeWeightPropType, G, BufferEdge), G);
              boost::remove_edge(nPredIndx, nCurrentIndx2, G);
				
              // No incoming edge of the second node remains. The outgoing 
              // edges will be removed so that it is isolated from the rest 
              // of the graph.
              boost::clear_out_edges(nCurrentIndx2, G);
				
              pnReversedNodeMap[nCurrentIndx1] = - pnReversedNodeMap[nCurrentIndx1];
              pnReversedNodeMap[nCurrentIndx2] = 0;
				
              nSuccIndx = nCurrentIndx1;
              nCurrentIndx1 = nPredIndx;
              nCurrentIndx2 = rPreds[nCurrentIndx1];
              nPredIndx = rPreds[nCurrentIndx2];
            }
			
          nSuccIndx = nCurrentIndx1;
          nCurrentIndx1 = nCurrentIndx2;
          // Pred is the same
			
        }
    }
	
  pnKShorthestPathPreds[nSuccIndx] = nSourceNode;
	
  // Reverse the last edge that emenates from the source
  boost::tie(BufferEdge, bDummyEdgeCheck) = 
    boost::edge(nCurrentIndx1, nSuccIndx, G);
  boost::add_edge(nSuccIndx, nCurrentIndx1, 
                  -boost::get( EdgeWeightPropType, G, BufferEdge), G);
  boost::remove_edge(nCurrentIndx1, nSuccIndx, G);
}
