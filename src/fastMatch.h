
struct idxStruct
{
   int  from;
   int  to;
};

int lowerBound(double x,double *val,int first, int length);

int upperBound(double x,double *val,int first, int length);

SEXP fastMatch(SEXP x, SEXP y, SEXP xidx, SEXP yidx, SEXP xolength, SEXP tol);