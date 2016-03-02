#ifndef Types_H
#define Types_H
#include <vector>

//#define HIGH_PRECISION

#ifdef HIGH_PRECISION
    typedef double ScalarType;
#else
    typedef float ScalarType;
#endif


#endif /* Types_H */
