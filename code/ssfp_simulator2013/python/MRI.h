#ifndef __PYX_HAVE__MRI
#define __PYX_HAVE__MRI


#ifndef __PYX_HAVE_API__MRI

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(PyObject) *MxySignal(double, double, double, double, int, double, double, double, double, double);

#endif /* !__PYX_HAVE_API__MRI */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initMRI(void);
#else
PyMODINIT_FUNC PyInit_MRI(void);
#endif

#endif /* !__PYX_HAVE__MRI */
