
struct idxStruct
{
   int  from;
   int  to;
};


SEXP fastMatch(SEXP x, SEXP y, SEXP xidx, SEXP yidx, SEXP xolength, SEXP tol);