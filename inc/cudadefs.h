#ifdef DEBUG
#define CUDA_CALL(x) do { if((x) != cudaSuccess) {\
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CUDA_CALL(x) (x)
#endif /*DEBUG*/
#ifdef DEBUG
#define CURAND_CALL(x) do { if((x) != CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CURAND_CALL(x) (x)
#endif /*DEBUG*/
#ifdef DEBUG
#define CUFFT_CALL(x) do { int err=(x); if(err != CUFFT_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
	printf("Error code: %d\n",err);\
        exit(EXIT_FAILURE);}} while(0)
#else
#define CUFFT_CALL(x) (x)
#endif /*DEBUG*/
