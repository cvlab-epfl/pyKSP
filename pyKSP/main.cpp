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

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "global.h"
#include "ksp_graph.h"
#include "ksp_computer.h"


static struct option long_options[] = {
  {"help", 0, 0, 'h'},
  {"output", 1, 0, 'o'},
  {0, 0, 0, 0}
};

static const char *explanation[] = {
  "display help",
  "output file name"
};

void print_help(char *program) {
  printf("Usage: %s [OPTIONS] <data_file>\n", program);
  printf("  where OPTIONS are:\n");
  int index = 0;
  while (long_options[index].name != 0) {
    printf("    ");
    if (long_options[index].val != 0)
      printf("-%c", char(long_options[index].val));
    printf("\t--%s  \t%s\n", long_options[index].name, explanation[index]);
    index++;
  }
}


int main(int argc, char **argv)
{
  printf("--------------------------------------------------------------------------\n");
  printf("Multiple Object Tracker using K-Shortest Paths. (c) 2009-2011 CVLab - EPFL\n");
  printf("All rights reserved.\n");
  printf("Written by Engin Turetken and Jerome Berclaz.\n");
  printf("--------------------------------------------------------------------------\n\n");
  char c;
  const char *config_file_name = NULL;
  const char *output_file_name = "ksp-out.dat";
  size_t buffer_size = 4096;
  char *buffer = new char[buffer_size];

  int first_frame         = 0;
  int last_frame          = 0;
  int grid_width          = 0;
  int grid_height         = 0;
  int grid_size           = 0;
  int nbr_frames          = 0;
  int neighborhood_size   = 0;
  unsigned int nbr_nodes  = 0;
  int depth               = DEFAULT_DEPTH;
  int max_traj            = DEFAULT_MAX_TRAJ;
  char *input_format      = NULL;
  float *input_data       = NULL;
  std::vector<int> access_points;

  // read command line parameters
  while (1) {
    int option_index = 0;

    c = getopt_long (argc, argv, "ho:", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      switch (option_index) {
      default:
        printf("Unrecognized option.\n");
        exit(0);
      }

      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
    case 'h':
      print_help(argv[0]);
      break;
    case 'o':
      output_file_name = optarg;
      break;
    default:
      printf("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if (optind < argc) {
    config_file_name = argv[optind++];
    if (optind < argc) {
      printf("Unknown options: ");
      while (optind < argc)
        printf("%s ", argv[optind++]);
      printf("\nThey will be ignored\n");
    }
  }

  if (!config_file_name) {
    printf("Error: no input file specified.\n");
    print_help(argv[0]);
    return 0;
  }


  // open config file
  FILE *config_file = fopen(config_file_name, "r");
  if (!config_file) {
    printf("Error: unable to open configuration file '%s'.\n", config_file_name);
    return 1;
  }
  
  int res;
  do {
    res = getline(&buffer, &buffer_size, config_file);
    if (res > 0) {
      char *command = strtok(buffer, " ");
      if (command[0] != '#' && command[0] != '\n') { // ignore comments and empty lines
      
        if (strcmp(command, "GRID") == 0) {
          grid_width = atoi(strtok(NULL, " "));
          grid_height = atoi(strtok(NULL, " "));
        }
        else if (strcmp(command, "FRAMES") == 0) {
          first_frame = atoi(strtok(NULL, " "));
          last_frame  = atoi(strtok(NULL, " "));
        }
        else if (strcmp(command, "ACCESS_POINTS") == 0) {
          char *ap = strtok(NULL, " ");
          char *cell = strtok(ap, ",");
          while (cell) {
            access_points.push_back(atoi(cell));
            cell = strtok(NULL, ",");
          }
        }
        else if (strcmp(command, "DEPTH") == 0) {
          depth = atoi(strtok(NULL, " "));
        }
        else if (strcmp(command, "MAX_TRAJ") == 0) {
          max_traj = atoi(strtok(NULL, " "));
        }
        else if (strcmp(command, "INPUT_FORMAT") == 0) {
          char *inf = strtok(NULL, " \n");
          input_format = new char[strlen(inf) + 1];
          strcpy(input_format, inf);
        }
        else {
          printf("Unknown command '%s' in configuration file '%s'.\n", command,
                 config_file_name);
          return 1;
        }
      }
    }
  } while (res >= 0);
  fclose(config_file);  

  nbr_frames = last_frame - first_frame + 1;
  if (nbr_frames <= 0) {
    printf("Error: No frame is going to be processed.\n");
    return 1;
  }
  neighborhood_size = depth * 2 + 1;
  grid_size = grid_width * grid_height;
  nbr_nodes = grid_size * nbr_frames;


  // print information message
  printf("Input file: '%s'\n", config_file_name);
  printf("===\n");
  printf("Grid: %dx%d\n", grid_width, grid_height);
  printf("%d frames: %d - %d\n", nbr_frames, first_frame, last_frame);
  printf("Depth: %d\n", depth);
  printf("Maximum number of trajectories: %d\n", max_traj);
  printf("%d access point(s)\n", int(access_points.size()));
  printf("Input files: '%s'\n", input_format);
  printf("Output file: '%s'\n===\n\n", output_file_name);

  
  // reading input files
  printf("Reading input files...");
  fflush(stdout);
  input_data = new float[nbr_nodes];
  float *frame_data = input_data;
  int frame;
  for (int f=first_frame; f<=last_frame; f++) {
    sprintf(buffer, input_format, f);
    FILE *frame_file = fopen(buffer, "r");
    if (!frame_file) {
      printf("Error: unable to open input file '%s'\n", buffer);
      return 1;
    }
    for (int i=0; i<grid_size; i++) {
      if (fscanf(frame_file, "%d %f\n", &frame, frame_data + i) != 2) {
        printf("Error while reading input file '%s', position %d\n", buffer, i);
        return 1;
      }
    }
    fclose(frame_file);
    frame_data += grid_size;
  }
  printf("\t\t\t[ok]\n");

  
  // taking the log of probabilities
  printf("Taking the log of probabilities...");
  fflush(stdout);
  float min_prob_log = -log( MIN_OCCUR_PROB / (1 - MIN_OCCUR_PROB) );
  float max_prob_log = -log( MAX_OCCUR_PROB / (1 - MAX_OCCUR_PROB) );
  float proba;
  for ( unsigned int i = 0; i < nbr_nodes; i++ ) {
    proba = input_data[i];
    
    if ( proba < MIN_OCCUR_PROB ) input_data[i] = min_prob_log;
    else if ( proba > MAX_OCCUR_PROB ) input_data[i] = max_prob_log;
    else input_data[i] = -log( proba / (1 - proba) );
  }
  printf("\t[ok]\n");


  // constructing the graph
  printf("Constructing the graph...");
  fflush(stdout);
  KShorthestPathGraph ksp_graph(input_data,
                                grid_width,
                                grid_height,
                                nbr_frames,
                                neighborhood_size,
                                access_points);
  printf("\t\t[ok]\n");


/////TEST
//KShorthestPathGraph ksp_graph1(ksp_graph.loc_pfData,ksp_graph.loc_nDataWidth,ksp_graph.loc_nDataHeight,ksp_graph.loc_nDataDepth,neighborhood_size,ksp_graph.loc_pnSrcAndDstNeighIndices); 

/////  

  // running the K-Shorthest node-disjoint paths algorithm                     
  printf("Running optimization...");
  fflush(stdout);
//KShorthestPathGraph ksp_graph1(input_data,
//                                grid_width,
//                                grid_height,
//                                nbr_frames,
//                                neighborhood_size,
//                                access_points);
								
								

//std::cout<<&ksp_graph<<" "<<&ksp_graph1<< std::endl;
////TEST with object

std::vector<int>  v = ksp_graph.getPath(0,20);

/////



								
unsigned char *labeled_objects = (unsigned char*)calloc(ksp_graph.GetNoOfVertices(), 
                                                          sizeof(unsigned char));
int nbr_paths = KShorthestPathComputer::ComputeKShorthestNodeDisjointPaths(ksp_graph, 
                                                                             max_traj, 
                                                                             MAX_PATH_LENGTH, 
                                                                             labeled_objects);
//																			 
																			 
 printf("\t\t\t[ok]\n");
  printf("\nFound K = %d shortest paths.\n", nbr_paths);
  
unsigned char *labeled_objects1 = (unsigned char*)calloc(ksp_graph.GetNoOfVertices(), 
                                                          sizeof(unsigned char));
  
 int nbr_paths1 = KShorthestPathComputer::ComputeKShorthestNodeDisjointPaths(ksp_graph, 
                                                                             max_traj, 
                                                                             MAX_PATH_LENGTH, 
                                                                             labeled_objects1);
  printf("\t\t\t[ok]\n");
  printf("\nFound K = %d shortest paths.\n", nbr_paths);

 
 
  // saving results
  printf("\nSaving results...");
  fflush(stdout);
  
  FILE *output_file = fopen(output_file_name, "w");	
  if (!output_file) {
    printf("Error: unable to open file '%s' for writing.\n", output_file_name);
    return 1;
  }
  fprintf(output_file, "%d\n", nbr_paths);

  unsigned int *positions = new unsigned int[nbr_paths];
  unsigned char *write_ptr = labeled_objects;
  for (int f=first_frame; f<=last_frame; f++) {
    for (int i=0; i<nbr_paths; i++) positions[i] = -1;
    for (int i=0; i<grid_size; i++) {
      if (write_ptr[i]) positions[write_ptr[i] - 1] = i;
    }
    fprintf(output_file, "%d\t", f);
    for (int i=0; i<nbr_paths; i++) fprintf(output_file, "%d\t", positions[i]);
    fprintf(output_file, "\n");

    write_ptr += grid_size;
  }
  fclose(output_file);
  delete [] positions;
  printf("\t\t\t[ok]\n");


  // releasing memory
  free(labeled_objects);
  delete [] input_data;
  delete [] buffer;
  if (input_format) delete [] input_format;
  
  return 0;
}

// Local Variables:
// mode: c++
// compile-command: "make -C ."
// End:
