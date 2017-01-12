This is a python wrapper written by Pierre Baque for 
“Multiple Object Tracker Using K-Shortest Paths”

if you want to modify the source code and recompile the wrapper, use the command:

python setup.py build_ext --inplace



Multiple Object Tracker Using K-Shortest Paths
==============================================
http://cvlab.epfl.ch/software/ksp


INTRODUCTION
------------

  This code implements a multiple object tracker based on the
  k-shortest paths algorithm. Its input consists in a set of
  probabilistic occupancy maps, that is, for every time frame, a set
  of occupancy probabilities, one for each of the potential target
  locations. Such input data is usually produced by an object
  detector.

  The KSP object tracking algorithm is able to track an unknown and
  varying number of objects. Unlike other related methods, it operates
  on the full set of potential target locations, and not just on the
  detections themselves. This characteristic allows it to handle
  missing detections and false positives well.

INSTALLATION AND TEST
---------------------

  The program should compile on any modern GNU/Linux system. Its only
  external library requirement is the Boost library, which needs to be
  installed along with its headers. To compile the program, simply
  use the 'make' command.
        
  To run a simple test on occupancy maps provided in the archive just
  type

    ./ksp test.ksp

  It will generate a result file 'ksp-out.dat' in the current
  directory. An alternate output file name can be specified with the
  option -o.

CONFIGURATION FILE
------------------

  The program relies on a configuration file that contains all the
  parameters needed for applying the tracker. Below is a description
  of the syntax of this file. Note that the lines starting with a
  character '#' are considered as comments and ignored.

  GRID <grid_width> <grid_height>

    The program assumes that the potential target locations are
    distributed as a regular grid. This command specifies the
    dimensions of this grid.

  FRAMES <first_frame> <last_frame>

    This command specifies the time window of the sequence to process.

  ACCESS_POINTS <pos1,pos2,pos3,...>

    This command specifies which - if any - cells of the grid can act
    as entrance and exit point, so that the number of objects can
    vary. In an indoor people tracking application, those would
    typically be locations close to doors, from which people can enter
    or exit. The set is given as a comma separated list of cells. A
    cell position is numbered by its coordinates in the grid:
                position = x + y * grid_width
    where the coordinates are defined in
                [0, grid_width[ x [0, grid_height[ .
    If there are no access point in the scene, this command should be
    commented out.

  DEPTH <depth>

    The depth represents the maximum allowed movement in the grid
    between two consecutive time frames. It is measured in grid
    cells. A value of 1 means that you can reach the next grid cell at
    most, while 2 means that you can travel 2 grid cells in any
    direction.

  MAX_TRAJ <maximum_number_of_trajectories>

    This value represents an upper bound on the number of trajectories
    retrieved by the algorithm. For implementation reasons, it cannot
    be larger than 255.

  INPUT_FORMAT <input_file_format>

    This string specifies the file format of the input data. Each time
    frame is stored in a different file, and every file contains, at
    each line, the number of the grid cell and its probability of
    occupancy. Thus the input file format should contain a special
    character %d that will be replaced by the frame number, in the
    standard printf way.

  An example configuration file named 'test.ksp' is provided for
  convenience. It performs tracking on the probabilistic occupancy
  maps stored in the folder 'proba'. Note that this data has been
  generated with the POM people detector available at
  http://cvlab.epfl.ch/software/pom .

REFERENCES
----------

  For more information about the KSP algorithm, please check the
  following article:

  Jerome Berclaz, Francois Fleuret, Engin Turetken and Pascal Fua,
  "Multiple Object Tracking using K-Shortest Paths Optimization", IEEE
  Transactions on Pattern Analysis and Machine Intelligence (TPAMI),
  2011, to be published.

LICENSE
-------

  The source code is available for academic purposes only and is
  distributed under a proprietary non-commercial license. Please read
  the included file 'license.txt' for the full terms of the
  license. If you are interested in using this algorithm in a
  commercial product, you can contact us to purchase a commercial
  license.

CONTACT
-------

  Please mail pom@epfl.ch for bug reports, comments and questions.
