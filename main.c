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
//    GxB_Global_Option_set(GxB_BURBLE, true);

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

    GrB_Vector v = NULL;
    GrB_Vector_new(&v, GrB_BOOL, 1073);
    GrB_Vector_setElement(v, true, 3);

    GrB_Matrix M = NULL;
    FILE* f = NULL;
    f = fopen("../matrix.mm", "r");
    LAGraph_mmread(&M, f);

    GxB_Vector_fprint(v, "v0", GxB_SUMMARY, stdout);
    GrB_vxm(v, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, v, M, GrB_NULL);
    GxB_Vector_fprint(v, "v1", GxB_SUMMARY, stdout);
    GrB_vxm(v, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, v, M, GrB_NULL);
    GxB_Vector_fprint(v, "v2", GxB_SUMMARY, stdout);

    // Cleanup
    LAGraph_finalize();

    return 0;
}
