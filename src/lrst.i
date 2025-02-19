/*
lrst.i: SWIG module defining the Python wrapper for our C++ modules
[Ref.: https://www.swig.org/Doc1.3/Python.html]
 */

%module lrst  // Module name

// Headers
%{
#include <string>
#include <vector>
#include "input_parameters.h"
#include "output_data.h"
#include "module_caller.h"
%}

%apply long long { int64_t };  // Maps int64_t to long long in Python
%apply unsigned long long { uint64_t };  // Maps uint64_t to unsigned long long in Python
%apply long { long int };  // Maps long int to long in Python

// Map uint64_t arrays to Python lists
%typemap(out) uint64_t[ANY] {
    PyObject *list = PyList_New($1_dim0);
    for (int i = 0; i < $1_dim0; i++) {
        PyList_SetItem(list, i, PyLong_FromUnsignedLongLong($1[i]));
    }
    $result = list;
}

// Map std::string arrays to Python lists
%typemap(out) std::string[ANY] {
    PyObject *list = PyList_New($1_dim0);
    for (int i = 0; i < $1_dim0; i++) {
        PyList_SetItem(list, i, PyUnicode_FromString($1[i].c_str()));
    }
    $result = list;
}

// Map std::map<std::string, std::map<int32_t, std::tuple<char, char, double,
// int, bool>>> to Python dictionary
%typemap(out) std::map<std::string, std::map<int32_t, std::tuple<char, char, double, int, bool>>> {
    PyObject *dict = PyDict_New();
    for (auto const &it : $1) {
        PyObject *inner_dict = PyDict_New();
        for (auto const &inner_it : it.second) {
            PyObject *tuple = PyTuple_Pack(5, 
                                           PyUnicode_FromStringAndSize(&std::get<0>(inner_it.second), 1), 
                                           PyUnicode_FromStringAndSize(&std::get<1>(inner_it.second), 1), 
                                           PyFloat_FromDouble(std::get<2>(inner_it.second)),
                                           PyLong_FromLong(std::get<3>(inner_it.second)),
                                           PyBool_FromLong(std::get<4>(inner_it.second)));
            PyDict_SetItem(inner_dict, PyLong_FromLong(inner_it.first), tuple);
        }
        PyDict_SetItem(dict, PyUnicode_FromString(it.first.c_str()), inner_dict);
    }
    $result = dict;
}

// Map std::map<char, int> to Python dictionary
%typemap(out) std::map<char, int> {
    PyObject *dict = PyDict_New();
    for (auto const &it : $1) {
        PyDict_SetItem(dict, PyUnicode_FromStringAndSize(&it.first, 1), PyLong_FromLong(it.second));
    }
    $result = dict;
}

%include <std_string.i>
%include <stdint.i>
%include <std_vector.i>

%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(Int2DVector) std::vector<std::vector<int>>;
%template(StringVector) std::vector<std::string>;
%template(CharVector) std::vector<char>;

// These are the header functions wrapped by our lrst module (Like an 'import')
%include "input_parameters.h"  // Contains InputPara for passing parameters to C++
%include "output_data.h"  // Contains data structures for storing statistics for each file type
%include "module_caller.h"  // Functions for calling the C++ statistics computation modules
