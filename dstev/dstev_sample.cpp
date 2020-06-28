// dstev


double zz[d*d];
double work[2*d-2];
char jobz='V';
int info;
dstev_( jobz, d, &a[0],&b[0],  zz, d, work, info);



double* zz = new double [restrict_dim*restrict_dim];
double* work = new double[(2*restrict_dim-2)];
char jobz='V';
int info;
dstev_( jobz, restrict_dim, &a[0],&b[0],  zz, restrict_dim, work, info);
