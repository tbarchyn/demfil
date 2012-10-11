// Generic filter program for performing 'focal statistics' in parallel with OpenMP
// 02 Jan 2012

// -------------------------------------------------------------------------------
// GLOBAL VARIABLES
const int max_nrow = 8000;              // these need to be adjusted if compiling on
const int max_ncol = 8000;              // systems with limited memory capacity
double in [max_nrow][max_ncol];         // input and output arrays, static allocation
double out [max_nrow][max_ncol];

const int max_filsize = 1001;           // this must be odd!
const int cen_i = max_filsize / 2;      // set the center coordinates
const int cen_j = max_filsize / 2;
bool fil[max_filsize][max_filsize];     // boolean filter mask, true = included, false = excluded

int nrows, ncols;                       // global constants for number of rows and number of columns

ostringstream infile;                   // input file
ostringstream outfile;                  // output file
ostringstream funcode;                  // function code
double rad;                             // radius of filter circle
double nodataflag;                      // no data flag value
char xllcorner[100];                    // set character array for projection parameters
char yllcorner[100];
char cellsize[100];
double nontoxic_frac = 1.0;             // fraction of nontoxic values required
                                        // defaults to 1.0
