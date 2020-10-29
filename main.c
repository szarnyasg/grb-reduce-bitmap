#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING 1

#include <GraphBLAS.h>
#include <LAGraph.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#undef NDEBUG
#include <assert.h>
#include <string.h>

#define LAGRAPH_FREE_ALL            \
{                                   \
    /* TODO */                      \
}


// override to exit on problem to use it in OpenMP region
#undef LAGRAPH_ERROR
#define LAGRAPH_ERROR(message,info)                                         \
{                                                                           \
    fprintf (stderr, "LAGraph error: %s\n[%d]\n%s\nFile: %s Line: %d\n",    \
        message, info, GrB_error ( ), __FILE__, __LINE__) ;                 \
    LAGRAPH_FREE_ALL ;                                                      \
    exit (info) ;                                                         \
}

int main() {
    LAGraph_init();
    GxB_Global_Option_set(GxB_BURBLE, true);

    // print version
    char *date, *compile_date, *compile_time;
    int version[3];
    GxB_Global_Option_get(GxB_LIBRARY_VERSION, version);
    GxB_Global_Option_get(GxB_LIBRARY_DATE, &date);
    GxB_Global_Option_get(GxB_LIBRARY_COMPILE_DATE, &compile_date);
    GxB_Global_Option_get(GxB_LIBRARY_COMPILE_TIME, &compile_time);
    printf("Library version %d.%d.%d\n", version[0], version[1], version[2]);
    printf("Library date: %s\n", date);
    printf("Compiled at %s on %s\n", compile_time, compile_date);

    LAGraph_set_nthreads(1);

    uint64_t indices[] = {1, 2, 4, 5, 7, 11, 12, 13, 15, 18, 19, 20, 27, 32, 33, 35, 37, 41, 46, 50, 52, 53, 55, 57, 58, 61, 62, 63, 65, 66, 69, 70, 72, 73, 74, 75, 78, 79, 81, 84, 86, 87, 90, 91, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 107, 108, 109, 110, 115, 116, 117, 118, 120, 123, 129, 131, 132, 133, 134, 136, 140, 145, 146, 149, 152, 153, 154, 156, 158, 159, 160, 161, 163, 164, 165, 166, 168, 169, 172, 176, 177, 181, 184, 186, 187, 189, 191, 193, 194, 195, 197, 200, 201, 202, 203, 204, 205, 208, 209, 210, 211, 216, 217, 218, 219, 224, 225, 229, 230, 232, 235, 236, 238, 239, 242, 243};
    uint64_t nvals = sizeof(indices)/sizeof(uint64_t);
    uint64_t n = 1000;
    for (uint64_t i = 0; i < nvals; i++) {
        indices[i]--;
    }

    GrB_Matrix emptyMatrix = NULL;
    GrB_Matrix_new(&emptyMatrix, GrB_UINT64, n, n);

    GrB_Vector v = NULL;
    GrB_Vector_new(&v, GrB_UINT64, n);
    GrB_Vector_build(v, indices, indices, nvals, GrB_PLUS_UINT64);

    GrB_Vector w = NULL;
    GrB_Vector_new(&w, GrB_UINT64, n);
    GrB_Vector_assign_UINT64(w, v, GrB_NULL, 0, GrB_ALL, 0, GrB_NULL);

    GrB_Matrix_reduce_Monoid(w, v, GrB_NULL, GrB_PLUS_MONOID_UINT64, emptyMatrix, GrB_NULL);

    // Cleanup
    LAGr_free(&v);
    LAGr_free(&w);
    LAGr_free(&emptyMatrix);

    LAGraph_finalize();

    return 0;
}
