/*
F5_module.cpp:
Class for calling FAST5 statistics modules.

*/

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "HDF5_module.h"
#include "ComFunction.h"
#include "H5Cpp.h"
using namespace H5;

const H5std_string FILE_NAME("/home/perdomoj/github/LongReadSum/output/SimpleCompound.h5");
const H5std_string DATASET_NAME("dset");
const int          NX   = 4; // dataset dimensions
const int          NY   = 6;
const int          RANK = 2;

int runExample(void)
{
    std::cout << "Running HDF5 Module..." << std::endl;
    // Try block to detect exceptions raised by any of the calls inside it
    try {
        // Turn off the auto-printing when failure occurs so that we can
        // handle the errors appropriately
        Exception::dontPrint();

        // Create a new file using the default property lists.
        H5File file(FILE_NAME, H5F_ACC_TRUNC);

        // Create the data space for the dataset.
        hsize_t dims[2]; // dataset dimensions
        dims[0] = NX;
        dims[1] = NY;
        DataSpace dataspace(RANK, dims);

        // Create the dataset.
        DataSet dataset = file.createDataSet(DATASET_NAME, PredType::STD_I32BE, dataspace);

    } // end of try block

    // catch failure caused by the H5File operations
    catch (FileIException error) {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch (DataSetIException error) {
        error.printErrorStack();
        return -1;
    }

    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException error) {
        error.printErrorStack();
        return -1;
    }

    return 0; // successfully terminated
}
