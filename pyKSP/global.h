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
/* Written by Engin Turetken and Jerome Berclaz.                        */
/*                                                                      */
/* http://cvlab.epfl.ch/research/body/surv                              */
/* Contact <pom@epfl.ch> for comments & bug reports.                    */
/************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

// constants to avoid numerical instabilities
#define MIN_OCCUR_PROB 1e-6
#define MAX_OCCUR_PROB 0.999999

// depth of the graph, that is number of cells that can be traveled in
// one time frame
#define DEFAULT_DEPTH     1

// maximum number of trajectories
#define DEFAULT_MAX_TRAJ  255

// this constant forces the algorithm to ignore very short trajectories
#define MAX_PATH_LENGTH   -45

#endif
