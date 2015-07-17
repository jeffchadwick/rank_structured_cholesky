int pcg_pysolver( const char* basename );
void pcg_pysolve(int nrow, const int* jc, const int* ir, const double* pr,
                 const double* rhs, double* solution);
