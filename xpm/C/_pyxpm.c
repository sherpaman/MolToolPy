#include <Python.h>
#include <numpy/arrayobject.h>
//#include <numpy/NumPy_macros.h>
#include <stdio.h>
#include "X11/xpm.h"
//#include "pyxpm.h"



static PyObject *pyxpm_read(PyObject *self, PyObject *args){
    
    char* fn;
    int i;
    
    XpmImage XPM;
    XpmInfo  XPM_info;
    
    
    npy_intp dims[2];
    
    PyIntObject *py_int;
    
    if (!PyArg_ParseTuple(args,"s", &fn)) 
        return NULL;
    
    //printf("Start Reading file : %s\n", fn);
    
    XpmReadFileToXpmImage(fn, &XPM, &XPM_info);
    
    printf("Done!\n");
    printf("%d x %d \n",XPM.height,XPM.width);
    
    // Create a Python List to return back
    
    PyListObject *ret = PyList_New((Py_ssize_t) 6);
    
    // Start Populating the return list with info from the XPM
    
    PyList_SetItem(ret, (Py_ssize_t) 0, PyInt_FromLong(XPM.height));
    PyList_SetItem(ret, (Py_ssize_t) 1, PyInt_FromLong(XPM.width));
    PyList_SetItem(ret, (Py_ssize_t) 2, PyInt_FromLong(XPM.cpp));
    PyList_SetItem(ret, (Py_ssize_t) 3, PyInt_FromLong(XPM.ncolors));
    
    dims[0] = XPM.height;
    dims[1] = XPM.width; 
    
    //printf(dims[0],dims[1]);

    /* 
     * List that will contain the Color Table of XPM:
     *  LIST = [ [ string0, color0] ,[string1, color1] ... [stringN, colorN] ]
     */
    PyListObject *colorTab = PyList_New((Py_ssize_t) XPM.ncolors);
    for (i=0;i<XPM.ncolors;i++){
        PyListObject *color = PyList_New((Py_ssize_t) 2);
        PyList_SetItem(color, (Py_ssize_t) 0, PyString_FromString(XPM.colorTable[i].string));
        PyList_SetItem(color, (Py_ssize_t) 1, PyString_FromString(XPM.colorTable[i].c_color));
        PyList_SetItem(colorTab, (Py_ssize_t) i, color); 
    }
    // Insert also the Color Table on the return list
    PyList_SetItem(ret,(Py_ssize_t) 4, colorTab);

    // Python Array that will contain all the data
    PyArrayObject *data = PyArray_SimpleNewFromData(2, dims, NPY_INT, XPM.data);
    
    // Insert the array in the return list
    PyList_SetItem(ret,(Py_ssize_t) 5, data);
    
    return ret;
}

static char module_docstring[] =
    "This module provides an interface for readinf XPM files using C.";
static char read_docstring[] =
    "Reads a XPM file.";

static PyMethodDef pyxpm_methods[] = {
    {"read", pyxpm_read, METH_VARARGS, read_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_pyxpm(void)
{
    PyObject *m = Py_InitModule3("_pyxpm", pyxpm_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load `numpy` functionality. */
    import_array();
};
