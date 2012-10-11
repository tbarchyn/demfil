// Generic read/write functions for ArcGIS Ascii files
// 02 Jan 2012

// -------------------------------------------------------------------------------
// READ ARCGIS ASCII FUNCTION
void read_ArcAscii_double ()
{
    /*
    The ArcGIS Ascii file has a header that consists of 6 rows of header info
    The file format is not standardized, so this tool may go haywire, but I
    believe it will work safely with arcGIS Ascii rasters created by this
    program and created by ArcGIS 10.

    Global variable requirements:
    This program will also write to global variables:
    in = input array, which is initialized with dimensions much bigger than necessary
    nrows = number of rows
    ncols = number of columns
    ascii_header = a string that contains GIS info, this gets written out again in the output file
    infile = ostringstream with the location of the file
    xllcorner, yllcorner, cellsize, nodataflag = projection parameters
    max_ncol, max_nrow = maximum size of in array
    */
    cout << "-------------------------------------------------------------" << endl;
    cout << "Beginning ArcGIS Ascii file read . . ." << endl;
    FILE *pFile;
    int scan = 0;       // dummy scan variable
    pFile = fopen ( infile.str().c_str() , "r");
    if (pFile == NULL)
    {
        cout << "ERROR: cannot find input file!" << endl;
        exit (10);
    }

    //=========================================================================================
    // Search 1: look for the number of columns
    int fail_cntr = 0;
    char read1[4];
    do
    {
        scan = fscanf (pFile, "%s", read1);
        fail_cntr++;
        if (fail_cntr == 100) { cout << "FILE READ FAILURE!, need 'ncols'" << endl; exit(2); }
    }
    while ( strcmp (read1, "ncols") != 0 &&
           strcmp (read1, "NCOLS") != 0);

    // The next integer should be the number of columns
    scan = fscanf (pFile, "%d", &ncols);

    //=========================================================================================
    // Search 2: look for the number of rows
    rewind(pFile);      // rewind to the beginning again
    fail_cntr = 0;
    char read2[4];
    do
    {
        scan = fscanf (pFile, "%s", read2);
        fail_cntr++;
        if (fail_cntr == 100) { cout << "FILE READ FAILURE!, need 'nrows'" << endl; exit(2); }
    }
    while ( strcmp (read2, "nrows") != 0 &&
           strcmp (read2, "NROWS") != 0);
    // The next integer should be the number of rows
    scan = fscanf (pFile, "%d", &nrows);

    //=========================================================================================
    // Search 3: look for the xllcorner
    rewind(pFile);      // rewind to the beginning again
    fail_cntr = 0;
    char read3[9];
    do
    {
        scan = fscanf (pFile, "%s", read3);
        fail_cntr++;
        if (fail_cntr == 100) { cout << "FILE READ FAILURE!, need 'xllcorner'" << endl; exit(2); }
    }
    while ( strcmp (read3, "xllcorner") != 0 &&
           strcmp (read3, "XLLCORNER") != 0);
    // The next string should be the xll corner
    scan = fscanf (pFile, "%s", xllcorner);

    //=========================================================================================
    // Search 4: look for the yllcorner
    rewind(pFile);      // rewind to the beginning again
    fail_cntr = 0;
    char read4[9];
    do
    {
        scan = fscanf (pFile, "%s", read4);
        fail_cntr++;
        if (fail_cntr == 100) { cout << "FILE READ FAILURE!, need 'yllcorner'" << endl; exit(2); }
    }
    while ( strcmp (read4, "yllcorner") != 0 &&
           strcmp (read4, "YLLCORNER") != 0);
    // The next string should be the yll corner
    scan = fscanf (pFile, "%s", yllcorner);

    //=========================================================================================
    // Search 5: look for the cellsize
    rewind(pFile);      // rewind to the beginning again
    fail_cntr = 0;
    char read5[9];
    do
    {
        scan = fscanf (pFile, "%s", read5);
        fail_cntr++;
        if (fail_cntr == 100) { cout << "FILE READ FAILURE!, need 'yllcorner'" << endl; exit(2); }
    }
    while ( strcmp (read5, "cellsize") != 0 &&
           strcmp (read5, "CELLSIZE") != 0);
    // The next float should be the yll corner
    scan = fscanf (pFile, "%s", cellsize);

    //=========================================================================================
    // Search 6: look for the nodata flag value (which should be -9999.0)
    rewind(pFile);      // rewind to the beginning again
    fail_cntr = 0;
    char read6[9];
    do
    {
        scan = fscanf (pFile, "%s", read6);
        fail_cntr++;
        if (fail_cntr == 100) { cout << "FILE READ FAILURE!, need 'nodata flag'" << endl; exit(2); }
    }
    while ( strcmp (read6, "nodata_value") != 0 &&
           strcmp (read6, "NODATA_value") != 0 &&
           strcmp (read6, "NODATA_VALUE") != 0);
    // The next float should be the yll corner
    scan = fscanf (pFile, "%lf", &nodataflag);

    //=========================================================================================
    // Check the size of the array to make sure its not too big
    if (nrows > max_nrow || ncols > max_ncol)
    {
        cout << "ERROR!!!: too many rows or columns: contact Tom and/or recompile with larger memory allocation" << endl;
        exit (7);
    }

    //=========================================================================================
    // OK, now we should be situated correctly to start reading in the body of the file
    // Read in the body of the file starting here, the next piece of text should be the start of the file
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            scan = fscanf (pFile, "%lf", &in[i][j]);
            if (scan != 1)
            {
                cout << "ERROR #2: problem with input file" << endl;
                exit(2);
            }
        }
    }
    fclose (pFile);

    // Check and potentially correct a nodata flag that is something other than -9999.0
    if (nodataflag != -9999.0)
    {
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                if (in[i][j] == nodataflag)
                {
                    in[i][j] = -9999.0;
                }
            }
        }
        cout << "WARNING: your Arc ASCII file has a nodata value of " << nodataflag << endl;
        cout << "Please note: I've changed it to -9999.0" << endl;
        nodataflag = -9999.0;
    }

    // Print the operation to the console
    time_t nowTime;
    struct tm * timeString;
    time (&nowTime);
    timeString = localtime (&nowTime);

    cout << "FILE read into memory successfully, Time: " << asctime(timeString) << endl;
    cout << "Number of rows: " << nrows << endl;
    cout << "Number of columns: " << ncols << endl;
    cout << "XLL corner: " << xllcorner << endl;
    cout << "YLL corner: " << yllcorner << endl;
    cout << "Cellsize: " << cellsize << endl;
    cout << "NODATA_value: " << nodataflag << endl;
    cout << "First number read: " << in[0][0] << endl;
    cout << "-------------------------------------------------------------" << endl;
}

// -------------------------------------------------------------------------------
// OUTPUT FUNCTION
void oput_ArcAscii_float()
{
    /*
    Arguments: the main argument is a stringstream of the output file

    The ArcGIS Ascii file has a header that consists of 6 rows of header info
    The file format is not standardized, so this tool may go haywire, but I
    believe it will work safely with arcGIS Ascii rasters created by this
    program and created by ArcGIS 10.

    Global requirements:
    This program will write out files from the 'out' array.
    The number of rows and columns are required in objects 'nrows' and 'ncols'
    The GIS data are required, in the following global variables
    xllcorner = x lower left corner
    yllcorner = y lower left cornter
    cellsize = cellsize
    nodataflag = should be -9999.0
    */
    cout << "-------------------------------------------------------------" << endl;
    cout << "Beginning ArcGIS Ascii file output . . ." << endl;
    // output the file
    FILE * pFile;
    pFile = fopen ( outfile.str().c_str() , "w");

    // Write the header
    ostringstream header;
    header << "ncols " << ncols << "\n";
    header << "nrows " << nrows << "\n";
    header << "xllcorner " << xllcorner << "\n";
    header << "yllcorner " << yllcorner << "\n";
    header << "cellsize " << cellsize << "\n";
    header << "NODATA_value " << nodataflag << "\n";

    // write header to the file stream
    fprintf (pFile, "%s", header.str().c_str());

    // Now, write the rest of the file out
    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < (ncols - 2); j++)
        {
            fprintf (pFile, "%f ", out[i][j]);
        }
        fprintf (pFile, "%f", out[i][ncols - 2]);
        fprintf (pFile, " %f", out[i][ncols - 1]);
        fprintf (pFile, "%s", "\n");    // endline character
    }
    fclose (pFile);

    // Print the operation to the console
    time_t nowTime;
    struct tm * timeString;
    time (&nowTime);
    timeString = localtime (&nowTime);

    cout << "Header values:\n" << header.str().c_str();
    cout << "FILE output to hard disk successfully, Time: " << asctime(timeString) << endl;
    cout << "-------------------------------------------------------------" << endl;
}



