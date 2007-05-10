#define DBL_MEMCPY(dest,src,n) memcpy((dest),(src),(n)*sizeof(double))
#define DBL_ZERO_MEMSET(dest,n) memset((dest),0,(n)*sizeof(double))
#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != 0) *(sp) = (s);} while(0)
