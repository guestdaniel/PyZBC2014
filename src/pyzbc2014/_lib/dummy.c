// src/pyzbc2014/_lib/dummy.c
#include <Python.h>
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "dummy_ext_for_platlib", NULL, -1, NULL, NULL, NULL, NULL, NULL
};
PyMODINIT_FUNC PyInit_dummy_ext_for_platlib(void) {
    return PyModule_Create(&moduledef);
}