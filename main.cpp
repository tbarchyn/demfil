// Generic filter program for performing 'focal statistics' in parallel with OpenMP
// 02 Jan 2012

/*
WARNING!! This program has no warranty, you are using it at your own risk!

Notes:
This program uses OpenMP to split the processing into chunks that are run
in separate threads. The granularity of the processing can be adjusted by modifying
the 'chunksize'. Reduced granularity will reduce the processing wait time at
the end of the loop while the loop waits for the last chunk to be finished.

Compiler flags (recomended) for windows:
-Wall -pedantic -O3 -fopenmp

Compiler flags (recomended) for linux:
-Wall -pedantic -O3 -fopenmp -no-stack-protector

Linked libraries:
On Codeblocks you have to make sure the library 'libgomp-1.dll' is in the local copy
of MinGW to compile with OpenMP. If it is missing, download the GUI installer for MinGW
from the internet and replace the entire MinGW folder with the updated version that should
include the missing library. On windows, a number of different libraries need to be in the
same directory as this exe for it to work, although this could vary with linker settings:
libgcc_s_dw2-1.dll, libgomp-1.dll, libstdc++-6.dll, pthreadGC2.dll

Notes: methods
ArcGIS performs focal statistics with similar methods to these. ArcGIS 10 has significantly
improved the performance of the ArcGIS engine for focal statistics; however, there could
be issues with accuracy, there originally was a bug with the function when ArcGIS 10 launched -
I'd check ArcGIS closely.

Arguments:
1) input file name
2) radius of test circle in cells
3) function code:
    m = mean
    s = sum
    f = minimum
    c = maximum
4) output file
5) required nontoxic proportion: the proportion of the filter circle required
    to be non-missing to output a value. This is optional, it is always set to 1.0,
    meaning all the circle is required to have values to report the output. The value must be between
    0.0 and 1.0. e.g., if nontoxic proportion is 1.0 and there is one missing value in the
    sliding circle, the output value will be missing (-9999.0).

All arguments are space separated in both linux and windows, for example, to run the program in
windows (assuming you compiled it with binary name: filter.exe) with the input file 'input.asc',
with a mean filter with radius 30 cells and output file of 'output.asc', you would type:

filter.exe input.asc 30 m output.asc

Or in Linux, sometimes you need to add './' to the beginning:

./filter.exe input.asc 30 m output.asc

*/

#define CHUNKSIZE 100      // define parallel chunksize for dynamic scheduling in OpenMP

#include <string.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <math.h>
#include <omp.h>            // note: for windows OpenMP requires special libraries, not
                            // found in stripped down versions of MinGW

using namespace std;

// Include header files
#include "tfil_globals.hpp"         // global variable declarations
#include "ascii_readwrite.hpp"      // functions for reading and writing ArcGIS ascii files
#include "tfil_func.hpp"            // main filter function

void print_man()
{
    // Print some help on the arguments if there are issues with the input file
    cout << "This program requires 4 arguments:\n"
        << "1) input file name (no spaces!), ArcGIS ASCII raster format\n"
        << "2) radius of filter circle in cells\n"
        << "3) function code, a single letter that is one of the following:\n"
        << "  m = mean\n  s = sum\n  f = minimum (floor)\n  c = maximum (ceiling)\n"
        << "4) output file name (no spaces!), ArcGIS ASCII raster format\n"
        << "5) optional last argument is proportion of filter window required to report a value\n\n"
        << "Example:\nI want to filter the file 'test.asc', with a mean filter with circle\n"
        << "with radius 30 cells, and output file name 'oput.asc', I also don't\n"
        << "care if up to half of the filter circle is missing data.\n"
        << "The name of the compiled version of this program is 'filter.exe'.\n"
        << "I would then type the following into the command line and press enter:\n\n"
        << "filter.exe test.asc 30 m oput.asc 0.5\n\n" << endl;
}

// -------------------------------------------------------------------------------
// MAIN
int main(int nArgs, char *pszArgs[])
{
    /* Arguments
    input file = string (no spaces!!!!)
    radius of filter circle = double
    function code = single character
    output file = string (no spaces!!!!!)
    toxicity code = float, the proportion of the circle that must be present
                    to report a value in the output raster
    */
    // Argument check
    if (nArgs < 5)
    {
        cout << "ERROR: not enough arguments!" << endl;
        print_man();
        exit(5);
    }

    // Read in the arguments
    infile << pszArgs[1];
    rad = atof (pszArgs[2]);
    funcode << pszArgs[3];
    outfile << pszArgs[4];
    if (nArgs > 5)      // ensure, if we are going to record the nontoxic fraction
    {                   // that the user actually input the value
        nontoxic_frac = atof (pszArgs[5]);
    }
    else
    {
        nontoxic_frac = 1.0;
    }
    // Print arguments to the console (defensive)
    // Determine present time
    time_t nowTime;
    struct tm * timeString;
    time (&nowTime);
    timeString = localtime (&nowTime);

    cout << "Welcome to the Tom's Filter program, Time: " << asctime(timeString)<< endl;
    cout << "This program was developed by the Hugenholtz Research Team" << endl;
    cout << "at the University of Lethbridge, Lethbridge, AB, Canada" << endl;
    cout << "Version compiled at: " << __TIMESTAMP__ << endl;
    cout << "This program has no warranty! It may not work as expected!" << endl;
    cout << "Arguments:\n  Input file: " << infile.str().c_str() << endl;
    cout << "  Radius of filter circle: " << rad << endl;
    cout << "  Function code: " << funcode.str().c_str() << endl;
    cout << "  Output file: " << outfile.str().c_str() << endl;
    cout << "  Required nontoxic fraction: " << nontoxic_frac << endl;

    init_tfil();                // initialize
    read_ArcAscii_double();     // read in the data from the file
    run_tfil();                 // run
    oput_ArcAscii_float();      // output the data in ArcAscii format

    return 0;
}




