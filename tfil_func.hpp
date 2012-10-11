// Generic filter program for performing 'focal statistics' in parallel with OpenMP
// 02 Jan 2012

// -------------------------------------------------------------------------------
// INITIALIZE FUNCTION: prepares the globals for setting the data
void init_tfil()
{
    // First, set all arrays to null parameters, just to be safe
    for (int i = 0; i < max_nrow; i++)
    {
        for (int j = 0; j < max_ncol; j++)
        {
            in[i][j] = -9999.0;
            out[i][j] = -9999.0;
        }
    }
    for (int i = 0; i < max_filsize; i++)
    {
        for (int j = 0; j < max_filsize; j++)
        {
            fil[i][j] = false;
        }
    }
    // Read in the arcgis Ascii file
}

// -------------------------------------------------------------------------------
// RUN FUNCTION
void run_tfil()
{
    //=========================================================================================
    // Create the filter boolean array
    // This array is 'true' if within the filter 'window'
    // First check the filter radius, it cannot be greater than the size of the array
    if (rad > cen_i)
    {
        cout << "INVALID Filter radius, recompile with larger static filter allocation!" << endl; exit (3);
    }

    double dist = 0.0;          // distance from focal cell: measured center to center
    int mask_sum = 0;           // the mask sum is the number of occurences of 'true'
    int min_i = max_filsize;    // begin by setting this to a very high value
    for (int i = 0; i < max_filsize; i++)
    {
        for (int j = 0; j < max_filsize; j++)
        {
            // Calculate difference between focal cell and location, and tag if less than radius
            // Note that the following comparison is 'less than or equal to': this matters for edge cells
            dist = double (sqrt ( ( (i - cen_i) * (i - cen_i) ) + ( (j - cen_j) * (j - cen_j) ) ) );
            if (dist <= rad)
            {
                fil[i][j] = true;
                mask_sum++;             // add one to the total
                if (i < min_i)
                {
                    min_i = i;
                }
            }
            else
            {
                fil[i][j] = false;
            }
        }
    }

    // Calculate the number of required values from each filter window
    // use ceiling to be conservative with this function
    int req_valcount = (int) ceil(nontoxic_frac * mask_sum);

    int edge_guard = cen_i - min_i;         // determine the number of cells to 'guard' on the edges

    // Pre-calculate start and finish coords for input array, these are subsequently used in
    // for loops with '<' conditionals (see below), thus, the loop will end one short of the
    // ending coordinates, leaving a strip of nodatas on the edge of the grids
    const int i_st = edge_guard;
    const int i_end = nrows - edge_guard;
    const int j_st = edge_guard;
    const int j_end = ncols - edge_guard;

    // Check to ensure the start and finish variables are not out of bounds!!
    if (i_st < 0 || i_st >= nrows || i_end < 0 || i_end >= nrows)
    {
        cout << "INVALID Filter radius!" << endl; exit (3);
    }
    if (j_st < 0 || j_st >= ncols || j_end < 0 || j_end >= ncols)
    {
        cout << "INVALID Filter radius!" << endl; exit (3);
    }

    // Pre-calculate starting and ending points for the filter mask array
    // The '+1' at the end of the ending coordinate is required because it is
    // part of a 'for loop' with a '<' conditional evaluation.
    const int i_f_st = min_i;
    const int i_f_end = cen_i + edge_guard + 1;
    const int j_f_st = min_i;
    const int j_f_end = cen_j + edge_guard + 1;

    // NEW algorithm also sets arrays with the coordinates for the trailing edge of the sliding
    // window and the leading edge of the window. The idea is to speed up calculations by keeping
    // the running sum and simply adding and subtracting from the front and back of the filter
    // window. This reduces the number of references substantially, and should speed up the
    // program significantly.
    // set the lookup arrays
    int trailing_i [max_filsize];
    int trailing_j [max_filsize];
    int leading_i [max_filsize];
    int leading_j [max_filsize];

    // now, loop through the filter and set the lookup arrays
    int i_tr = 0;   // thumb coordinate for the lookup, starting with zero
    for (int i = i_f_st; i < i_f_end; i++)
    {
        for (int j = j_f_st; j < j_f_end; j++)
        {
            // Starting a new row, we'll run across the filter and record the edge coords
            // Set the coordinates, this will be assessed when the function has a focal
            // cell that is one cell to the right, so, we need to subtract one from it.
            // These values are the offset from the focal cell, so some will be negative
            if (fil[i][j] && !fil[i][j-1])
            {
                trailing_i [i_tr] = i - cen_i;
                trailing_j [i_tr] = j - cen_j - 1;
            }
            // check for the leading coordinate
            if (fil[i][j] && !fil[i][j+1])
            {
                leading_i [i_tr] = i - cen_i;
                leading_j [i_tr] = j - cen_j;
            }
        }
        i_tr++;
    }
    int len_lkups = i_tr;       // save the length of the lookup array

    cout << "Beginning calculations with function code: " << funcode.str().c_str() << endl;

    //=========================================================================================
    // Calculation module: MEAN
    if (strcmp (funcode.str().c_str(), "m") == 0 || strcmp (funcode.str().c_str(), "M") == 0)
    {
        cout << "EXECUTING: mean . . ." << endl;
        // Use OpenMP to split the rows up into small chunks that are farmed
        // out to available processers dynamically as they are available.
        #pragma omp parallel for schedule (dynamic, CHUNKSIZE)
        for (int i = i_st; i < i_end; i++)
        {
            // Prepare some private variables, note that 'j' is also private to each processer
            double sub_val = 0.0;       // temp variables to for the sliding window part
            double add_val = 0.0;
            double runsum = 0.0;        // running sum
            int i_in, j_in;             // thumb coordinates
            int nontoxic_cntr = 0;         // nontoxic counter to track good values

            // START NEW ROW CALC SEQUENCE HERE
            // If this is a new row, we have to thumb over the whole filter mask
            // and properly calculate the mean and runsum
            runsum = 0.0;               // running sum for mean calculation
            int j = j_st;               // set 'j' to starting column
            i_in = i - edge_guard;      // coordinates that thumb over the input grid
            j_in = j - edge_guard;      // starting from minimum values, progressively looping up
            nontoxic_cntr = 0;             // nontoxic counter starts at zero

            // Main filter loop, thumbs over the mask, while updating input grid coords simultaneously
            for (int i_fil = i_f_st; i_fil < i_f_end; i_fil++)
            {
                j_in = j - edge_guard;      // reset j_in back to beginning column
                for (int j_fil = j_f_st; j_fil < j_f_end; j_fil++)
                {
                    // Check to see if the filter mask is 'true'
                    if ( fil[i_fil][j_fil] )
                    {
                        if (in[i_in][j_in] != -9999.0)  // check toxicity
                        {
                            nontoxic_cntr++;   // advance the nontoxic counter
                            //Add to the running sum, if inside the filter mask
                            runsum = runsum + in[i_in][j_in];
                        }
                    }
                    j_in++;     // increment input thumb coordinate
                }
                i_in++;     // increment input thumb coordinate
            }
            if (nontoxic_cntr < req_valcount)        // check to see if we have enough values
            {
                out[i][j] = -9999.0;        // fill in with 'nodata'
            }
            else
            {
                out[i][j] = runsum / nontoxic_cntr;      // calculate mean and store before moving on
            }
            // END NEW ROW CALC SEQUENCE HERE

            // ROW LOOP: continue to the right
            for (int j = (j_st + 1); j < j_end; j++)
            {
                // START SHIFTING CALC SEQUENCE HERE: the rest of the calcs will be this type
                // Loop down the lookups and add and subtract values
                for (int i_tr = 0; i_tr < len_lkups; i_tr++)
                {
                    // Perform the subtraction from the running sum
                    sub_val = in [ (i + trailing_i[i_tr]) ][ (j + trailing_j[i_tr]) ];
                    if (sub_val != -9999.0)
                    {
                        nontoxic_cntr--;            // decrement the nontoxic counter
                        runsum = runsum - sub_val; // subtract val from running sum
                    }
                    // Perform the addition to the running sum
                    add_val = in [ (i + leading_i[i_tr]) ][ (j + leading_j[i_tr]) ];
                    if (add_val != -9999.0)
                    {
                        nontoxic_cntr++;            // increment the nontoxic counter
                        runsum = runsum + add_val;  // add value to the running sum
                    }
                }
                // Now, record the value in the output array, if we are non-toxic
                if (nontoxic_cntr >= req_valcount)        // check toxic counter
                {
                    out[i][j] = runsum / nontoxic_cntr;        // fill in with 'nodata'
                }
                // END SHIFTING CALC SEQUENCE HERE: on to the next column
            }
        }
    }


    //=========================================================================================
    // Calculation module: SUM
    else if (strcmp (funcode.str().c_str(), "s") == 0 || strcmp (funcode.str().c_str(), "S") == 0)
    {
        cout << "EXECUTING: sum . . ." << endl;
        // Use OpenMP to split the rows up into small chunks that are farmed
        // out to available processers dynamically as they are available.
        #pragma omp parallel for schedule (dynamic, CHUNKSIZE)
        for (int i = i_st; i < i_end; i++)
        {
            // Prepare some private variables, note that 'j' is also private to each processer
            double sub_val = 0.0;       // temp variables to for the sliding window part
            double add_val = 0.0;
            double runsum = 0.0;        // running sum
            int i_in, j_in;             // thumb coordinates
            int nontoxic_cntr = 0;         // nontoxic counter to track good values

            // START NEW ROW CALC SEQUENCE HERE
            // If this is a new row, we have to thumb over the whole filter mask
            // and properly calculate the mean and runsum
            runsum = 0.0;               // running sum for mean calculation
            int j = j_st;               // set 'j' to starting column
            i_in = i - edge_guard;      // coordinates that thumb over the input grid
            j_in = j - edge_guard;      // starting from minimum values, progressively looping up
            nontoxic_cntr = 0;             // nontoxic counter starts at zero

            // Main filter loop, thumbs over the mask, while updating input grid coords simultaneously
            for (int i_fil = i_f_st; i_fil < i_f_end; i_fil++)
            {
                j_in = j - edge_guard;      // reset j_in back to beginning column
                for (int j_fil = j_f_st; j_fil < j_f_end; j_fil++)
                {
                    // Check to see if the filter mask is 'true'
                    if ( fil[i_fil][j_fil] )
                    {
                        if (in[i_in][j_in] != -9999.0)  // check toxicity
                        {
                            nontoxic_cntr++;   // advance the nontoxic counter
                            //Add to the running sum, if inside the filter mask
                            runsum = runsum + in[i_in][j_in];
                        }
                    }
                    j_in++;     // increment input thumb coordinate
                }
                i_in++;     // increment input thumb coordinate
            }
            if (nontoxic_cntr < req_valcount)        // check to see if we have enough values
            {
                out[i][j] = -9999.0;        // fill in with 'nodata'
            }
            else
            {
                out[i][j] = runsum;      // calculate sum and store before moving on
            }
            // END NEW ROW CALC SEQUENCE HERE

            // ROW LOOP: continue to the right
            for (int j = (j_st + 1); j < j_end; j++)
            {
                // START SHIFTING CALC SEQUENCE HERE: the rest of the calcs will be this type
                // Loop down the lookups and add and subtract values
                for (int i_tr = 0; i_tr < len_lkups; i_tr++)
                {
                    // Perform the subtraction from the running sum
                    sub_val = in [ (i + trailing_i[i_tr]) ][ (j + trailing_j[i_tr]) ];
                    if (sub_val != -9999.0)
                    {
                        nontoxic_cntr--;            // decrement the nontoxic counter
                        runsum = runsum - sub_val; // subtract val from running sum
                    }
                    // Perform the addition to the running sum
                    add_val = in [ (i + leading_i[i_tr]) ][ (j + leading_j[i_tr]) ];
                    if (add_val != -9999.0)
                    {
                        nontoxic_cntr++;            // increment the nontoxic counter
                        runsum = runsum + add_val;  // add value to the running sum
                    }
                }
                // Now, record the value in the output array, if we are non-toxic
                if (nontoxic_cntr >= req_valcount)        // check toxic counter
                {
                    out[i][j] = runsum;      // calculate sum and store before moving on
                }
                // END SHIFTING CALC SEQUENCE HERE: on to the next column
            }
        }
    }


    //=========================================================================================
    // Calculation module: MINIMUM
    else if (strcmp (funcode.str().c_str(), "f") == 0 || strcmp (funcode.str().c_str(), "F") == 0)
    {
        cout << "EXECUTING: minimum . . ." << endl;
        // Use OpenMP to split the rows up into small chunks that are farmed
        // out to available processers dynamically as they are available.
        #pragma omp parallel for schedule (dynamic, CHUNKSIZE)
        for (int i = i_st; i < i_end; i++)
        {
            // Prepare some private variables, note that 'j' is also private to each processer
            double sub_val = 0.0;       // temp variables to for the sliding window part
            double add_val = 0.0;
            double minval = 0.0;        // minimum value
            int i_in, j_in;             // thumb coordinates
            int nontoxic_cntr = 0;      // nontoxic counter to track good values
            int min_i = -1;              // track coordinates of the minimum value in the filter
            int min_j = -1;

            // START NEW ROW CALC SEQUENCE HERE
            // If this is a new row, we have to thumb over the whole filter mask
            // and properly calculate the mean and runsum
            minval = 99999999.9;        // set min value to a high value
            int j = j_st;               // set 'j' to starting column
            i_in = i - edge_guard;      // coordinates that thumb over the input grid
            j_in = j - edge_guard;      // starting from minimum values, progressively looping up
            nontoxic_cntr = 0;             // nontoxic counter starts at zero

            // Main filter loop, thumbs over the mask, while updating input grid coords simultaneously
            for (int i_fil = i_f_st; i_fil < i_f_end; i_fil++)
            {
                j_in = j - edge_guard;      // reset j_in back to beginning column
                for (int j_fil = j_f_st; j_fil < j_f_end; j_fil++)
                {
                    // Check to see if the filter mask is 'true'
                    if ( fil[i_fil][j_fil] )
                    {
                        if (in[i_in][j_in] != -9999.0)  // check toxicity
                        {
                            nontoxic_cntr++;   // advance the nontoxic counter
                            if (in[i_in][j_in] < minval)
                            {
                                minval = in[i_in][j_in];    // record new minimum value
                                min_i = i_in;               // record the coordinates
                                min_j = j_in;
                            }
                        }
                    }
                    j_in++;     // increment input thumb coordinate
                }
                i_in++;     // increment input thumb coordinate
            }
            if (nontoxic_cntr >= req_valcount)        // check to see if we have enough values
            {
                out[i][j] = minval;      // record the output value
            }
            // END NEW ROW CALC SEQUENCE HERE

            int i_lkup, j_lkup;     // set variables to record lookup coordinates
            // ROW LOOP: continue to the right
            for (int j = (j_st + 1); j < j_end; j++)
            {
                // START SHIFTING CALC SEQUENCE HERE: the rest of the calcs will be this type
                bool redo_normal = false;   // set flag to false
                // Loop down the lookups and assess the values as they come up
                for (int i_tr = 0; i_tr < len_lkups; i_tr++)
                {
                    // Jot down the lookup coordinates
                    i_lkup = i + trailing_i[i_tr];
                    j_lkup = j + trailing_j[i_tr];
                    sub_val = in [ i_lkup ][ j_lkup ];
                    if (sub_val != -9999.0)
                    {
                        nontoxic_cntr--;            // decrement the nontoxic counter
                        if (i_lkup == min_i && j_lkup == min_j)
                        {
                            redo_normal = true;     // if one of the values on the trailing
                                                    // edge is the minimum value, we have to
                                                    // redo the algorithm normally to find min
                        }
                    }
                    // Check the leading values
                    i_lkup = i + leading_i[i_tr];
                    j_lkup = j + leading_j[i_tr];
                    add_val = in [ i_lkup ][ j_lkup ];
                    if (add_val != -9999.0)
                    {
                        nontoxic_cntr++;            // increment the nontoxic counter
                        // record a new low value, but don't bother if we're already going
                        // to redo the focal cell
                        if (!redo_normal && add_val < minval)
                        {
                            minval = in[i_lkup][j_lkup];    // record new minimum value
                            min_i = i_lkup;
                            min_j = j_lkup;
                        }
                    }
                }
                // REDO normally: if the minimum value was found on the trailing
                // edge of the sliding window, we need to find a new minimum value normally
                if (redo_normal)
                {
                    minval = 99999999.9;        // set min value to a high value
                    i_in = i - edge_guard;      // coordinates that thumb over the input grid
                    j_in = j - edge_guard;      // starting from minimum values, progressively looping up
                    nontoxic_cntr = 0;          // nontoxic counter starts at zero

                    // Main filter loop, thumbs over the mask, while updating input grid coords simultaneously
                    for (int i_fil = i_f_st; i_fil < i_f_end; i_fil++)
                    {
                        j_in = j - edge_guard;      // reset j_in back to beginning column
                        for (int j_fil = j_f_st; j_fil < j_f_end; j_fil++)
                        {
                            // Check to see if the filter mask is 'true'
                            if ( fil[i_fil][j_fil] )
                            {
                                if (in[i_in][j_in] != -9999.0)  // check toxicity
                                {
                                    nontoxic_cntr++;   // advance the nontoxic counter
                                    if (in[i_in][j_in] < minval)
                                    {
                                        minval = in[i_in][j_in];    // record new minimum value
                                        min_i = i_in;
                                        min_j = j_in;
                                    }
                                }
                            }
                            j_in++;     // increment input thumb coordinate
                        }
                        i_in++;     // increment input thumb coordinate
                    }
                }
                // Now, record the value in the output array, if we are non-toxic
                if (nontoxic_cntr >= req_valcount)        // check to see if we have enough values
                {
                    out[i][j] = minval;      // record the output value
                }
                // END SHIFTING CALC SEQUENCE HERE: on to the next column
            }
        }
    }

    //=========================================================================================
    // Calculation module: MAXIMUM
    else if (strcmp (funcode.str().c_str(), "c") == 0 || strcmp (funcode.str().c_str(), "C") == 0)
    {
        cout << "EXECUTING: maximum . . ." << endl;
        // Use OpenMP to split the rows up into small chunks that are farmed
        // out to available processers dynamically as they are available.
        #pragma omp parallel for schedule (dynamic, CHUNKSIZE)
        for (int i = i_st; i < i_end; i++)
        {
            // Prepare some private variables, note that 'j' is also private to each processer
            double sub_val = 0.0;       // temp variables to for the sliding window part
            double add_val = 0.0;
            double maxval = 0.0;        // minimum value
            int i_in, j_in;             // thumb coordinates
            int nontoxic_cntr = 0;      // nontoxic counter to track good values
            int max_i = -1;              // track coordinates of the maximum values in the filter
            int max_j = -1;

            // START NEW ROW CALC SEQUENCE HERE
            // If this is a new row, we have to thumb over the whole filter mask
            // and properly calculate the mean and runsum
            maxval = -99999999.9;        // set max value to a low value
            int j = j_st;               // set 'j' to starting column
            i_in = i - edge_guard;      // coordinates that thumb over the input grid
            j_in = j - edge_guard;      // starting from minimum values, progressively looping up
            nontoxic_cntr = 0;             // nontoxic counter starts at zero

            // Main filter loop, thumbs over the mask, while updating input grid coords simultaneously
            for (int i_fil = i_f_st; i_fil < i_f_end; i_fil++)
            {
                j_in = j - edge_guard;      // reset j_in back to beginning column
                for (int j_fil = j_f_st; j_fil < j_f_end; j_fil++)
                {
                    // Check to see if the filter mask is 'true'
                    if ( fil[i_fil][j_fil] )
                    {
                        if (in[i_in][j_in] != -9999.0)  // check toxicity
                        {
                            nontoxic_cntr++;   // advance the nontoxic counter
                            if (in[i_in][j_in] > maxval)
                            {
                                maxval = in[i_in][j_in];    // record new max value
                                max_i = i_in;               // record the coordinates
                                max_j = j_in;
                            }
                        }
                    }
                    j_in++;     // increment input thumb coordinate
                }
                i_in++;     // increment input thumb coordinate
            }
            if (nontoxic_cntr >= req_valcount)        // check to see if we have enough values
            {
                out[i][j] = maxval;      // record the output value
            }
            // END NEW ROW CALC SEQUENCE HERE

            int i_lkup, j_lkup;     // set variables to record lookup coordinates
            // ROW LOOP: continue to the right
            for (int j = (j_st + 1); j < j_end; j++)
            {
                // START SHIFTING CALC SEQUENCE HERE: the rest of the calcs will be this type
                bool redo_normal = false;   // set flag to false
                // Loop down the lookups and assess the values as they come up
                for (int i_tr = 0; i_tr < len_lkups; i_tr++)
                {
                    // Jot down the lookup coordinates
                    i_lkup = i + trailing_i[i_tr];
                    j_lkup = j + trailing_j[i_tr];
                    sub_val = in [ i_lkup ][ j_lkup ];
                    if (sub_val != -9999.0)
                    {
                        nontoxic_cntr--;            // decrement the nontoxic counter
                        if (i_lkup == max_i && j_lkup == max_j)
                        {
                            redo_normal = true;     // if one of the values on the trailing
                                                    // edge is the maximum value, we have to
                                                    // redo the algorithm normally to find max
                        }
                    }
                    // Check the leading values
                    i_lkup = i + leading_i[i_tr];
                    j_lkup = j + leading_j[i_tr];
                    add_val = in [ i_lkup ][ j_lkup ];
                    if (add_val != -9999.0)
                    {
                        nontoxic_cntr++;            // increment the nontoxic counter
                        // record a new low value, but don't bother if we're already going
                        // to redo the focal cell
                        if (!redo_normal && add_val > maxval)
                        {
                            maxval = in[i_lkup][j_lkup];    // record new maximum value
                            max_i = i_lkup;
                            max_j = j_lkup;
                        }
                    }
                }
                // REDO normally: if the minimum value was found on the trailing
                // edge of the sliding window, we need to find a new minimum value normally
                if (redo_normal)
                {
                    maxval = -99999999.9;        // set min value to a high value
                    i_in = i - edge_guard;      // coordinates that thumb over the input grid
                    j_in = j - edge_guard;      // starting from minimum values, progressively looping up
                    nontoxic_cntr = 0;          // nontoxic counter starts at zero

                    // Main filter loop, thumbs over the mask, while updating input grid coords simultaneously
                    for (int i_fil = i_f_st; i_fil < i_f_end; i_fil++)
                    {
                        j_in = j - edge_guard;      // reset j_in back to beginning column
                        for (int j_fil = j_f_st; j_fil < j_f_end; j_fil++)
                        {
                            // Check to see if the filter mask is 'true'
                            if ( fil[i_fil][j_fil] )
                            {
                                if (in[i_in][j_in] != -9999.0)  // check toxicity
                                {
                                    nontoxic_cntr++;   // advance the nontoxic counter
                                    if (in[i_in][j_in] > maxval)
                                    {
                                        maxval = in[i_in][j_in];    // record new maximum value
                                        max_i = i_in;
                                        max_j = j_in;
                                    }
                                }
                            }
                            j_in++;     // increment input thumb coordinate
                        }
                        i_in++;     // increment input thumb coordinate
                    }
                }
                // Now, record the value in the output array, if we are non-toxic
                if (nontoxic_cntr >= req_valcount)        // check to see if we have enough values
                {
                    out[i][j] = maxval;      // record the output value
                }
                // END SHIFTING CALC SEQUENCE HERE: on to the next column
            }
        }
    }

    else
    {
        cout << "ERROR: I couldn't recognize your function code??" << endl;
    }
    cout << "Ending calculations with function code: " << funcode.str().c_str() << endl;
}





