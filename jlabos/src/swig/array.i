namespace plb {
    template<typename T> class MultiNTensorField3D;
}

%typemap(in)     (PRECOMP_T * IN_ARRAY1, int DIM1) {
        /* Functions from jni.h */
            $1 = (PRECOMP_T *) JCALL2(ARRAY_ELEMENTS, jenv, $input, 0);
                $2 = (int)    JCALL1(GetArrayLength,       jenv, $input);
}
%typemap(jni)    (PRECOMP_T * IN_ARRAY1, int DIM1) "jPRECOMP_TArray"
%typemap(jtype)  (PRECOMP_T * IN_ARRAY1, int DIM1) "PRECOMP_T[]"
%typemap(jstype) (PRECOMP_T * IN_ARRAY1, int DIM1) "PRECOMP_T[]"
%typemap(javain) (PRECOMP_T * IN_ARRAY1, int DIM1) "$javainput"

/*  %typemap(jtype) (PRECOMP_T* INPLACE_ARRAY1, int DIM1) */

/* Specify signature of method to handle */
/*%apply (PRECOMP_T * IN_ARRAY1, int DIM1)   { (PRECOMP_T * in_value, int size) };*/


/*%apply (PRECOMP_T* INPLACE_ARRAY1, int DIM1) {(PRECOMP_T* result, int size)};*/
/*%apply (PRECOMP_T* IN_ARRAY1, int DIM1) {(PRECOMP_T* result, int size)};*/


/*%apply (PRECOMP_T* IN_ARRAY1, int DIM1) {(PRECOMP_T* alpha, int size)};*/
/*%apply (PRECOMP_T* IN_ARRAY1, int DIM1) {(PRECOMP_T* alpha, int size)};*/


/*%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* externalVector, int numDimIs2)};*/
/*%apply (FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* externalVector, int numDimIs2)};*/

/*%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* velocity, int numDimIs2)};*/
/*%apply (FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* velocity, int numDimIs2)};*/

/*%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* populations, int numDimIsQ)};*/
/*%apply (FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* alpha, int size)};*/
/*

%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* externalVector, int numDimIs2)};
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* velocity, int numDimIs2)};
%apply(FLOAT_T* IN_ARRAY1, int DIM1) {(FLOAT_T* populations, int numDimIsQ)};
*/
