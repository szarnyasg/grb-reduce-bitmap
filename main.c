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


void advance_wavefront(GrB_Matrix HasCreator, GrB_Matrix HasCreatorTransposed,  GrB_Matrix ReplyOf, GrB_Matrix ReplyOfTransposed,  GrB_Matrix Knows, GrB_Vector frontier, GrB_Vector next, GrB_Vector seen, GrB_Index numPersons, int64_t comment_lower_limit, GrB_Index NumComments) {
    if (comment_lower_limit == -1) {
        GrB_vxm(next, seen, NULL, GxB_ANY_PAIR_BOOL, frontier, Knows, GrB_DESC_RC);
    } else {
        GxB_Scalar limit;
        GxB_Scalar_new(&limit, GrB_INT64);
        GxB_Scalar_setElement_INT64(limit, comment_lower_limit);

        // build selection matrix based on the frontier's content
        GrB_Matrix Sel;
        GrB_Matrix_new(&Sel, GrB_BOOL, numPersons, numPersons);

        GrB_Index nvals1;
        GrB_Vector_nvals(&nvals1, frontier);
        GrB_Index *I = (GrB_Index*) (LAGraph_malloc(nvals1, sizeof(GrB_Index)));
        bool *X = (bool*) (LAGraph_malloc(nvals1, sizeof(bool)));

        GrB_Index nvals2 = nvals1;
        GrB_Vector_extractTuples_BOOL(I, X, &nvals2, frontier);
        assert(nvals1 == nvals2);
        GrB_Matrix_build_BOOL(Sel, I, I, X, nvals1, GrB_LOR);

        GrB_Matrix M2;
        GrB_Matrix_new(&M2, GrB_BOOL, numPersons, NumComments);
        GrB_mxm(M2, NULL, NULL, GxB_ANY_PAIR_BOOL, Sel, HasCreatorTransposed, NULL);

        // direction 1
        GrB_Matrix M3a;
        GrB_Matrix_new(&M3a, GrB_UINT64, numPersons, NumComments);
        GrB_mxm(M3a, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, M2, ReplyOf, NULL);

        GrB_Matrix Interactions1;
        GrB_Matrix_new(&Interactions1, GrB_UINT64, numPersons, numPersons);
        GrB_mxm(Interactions1, Knows, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, M3a, HasCreator, NULL);


#ifndef NDEBUG
//            GxB_Matrix_fprint(Interactions1, "Interactions1 before select", GxB_SUMMARY, stdout);
#endif
//            GxB_Matrix_select(Interactions1, NULL, NULL, GxB_GT_THUNK, Interactions1, limit, NULL);
#ifndef NDEBUG
//            GxB_Matrix_fprint(Interactions1, "Interactions1 after select", GxB_SUMMARY, stdout);
#endif

        // direction 2
        GrB_Matrix M3b;
        GrB_Matrix_new(&M3b, GrB_UINT64, numPersons, NumComments);
        GrB_mxm(M3b, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, M2, ReplyOfTransposed, NULL);


        GrB_Matrix Interactions2;
        GrB_Matrix_new(&Interactions2, GrB_UINT64, numPersons, numPersons);
        GrB_mxm(Interactions2, Interactions1, GrB_NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, M3b, HasCreator, NULL);


#ifndef NDEBUG
//        GxB_Matrix_fprint(Interactions2, "Interactions2 before select", GxB_SUMMARY, stdout);
#endif
#ifndef NDEBUG
//            GxB_Matrix_fprint(Interactions2, "Interactions2 after select", GxB_SUMMARY, stdout);
#endif

        // Interactions1 = Interactions1 * Interactions2
        GrB_Matrix_eWiseMult_BinaryOp(Interactions1, NULL, NULL, GrB_MIN_UINT64, Interactions1, Interactions2, NULL);
        GxB_Matrix_select(Interactions1, NULL, NULL, GxB_GT_THUNK, Interactions1, limit, NULL);
#ifndef NDEBUG
//            GxB_Matrix_fprint(Interactions1, "Interactions final", GxB_SUMMARY, stdout);
#endif
        GrB_Matrix_reduce_BinaryOp(next, NULL, NULL, GxB_PAIR_BOOL, Interactions1, GrB_DESC_T0);

    }
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

    GrB_Matrix Knows, ReplyOf, HasCreator;
    LAGraph_binread(&Knows, "../knows.grb");
    LAGraph_binread(&ReplyOf, "../replyOf.grb");
    LAGraph_binread(&HasCreator, "../hasCreator.grb");

    GrB_Index NumComments;
    GrB_Matrix_nrows(&NumComments, ReplyOf);

    GrB_Index numPersons;
    GrB_Matrix_nrows(&numPersons, Knows);

    // hard-code input params for o1k data set

    GrB_Index p1 = 319;
    GrB_Index p2 = 809;
    int64_t comment_lower_limit = 1;

    // declare

    GrB_Vector frontier1 = NULL;
    GrB_Vector frontier2 = NULL;
    GrB_Vector next1 = NULL;
    GrB_Vector next2 = NULL;
    GrB_Vector intersection1 = NULL;
    GrB_Vector intersection2 = NULL;
    GrB_Vector seen1 = NULL;
    GrB_Vector seen2 = NULL;

    // create

    GrB_Vector_new(&frontier1, GrB_BOOL, numPersons);
    GrB_Vector_new(&frontier2, GrB_BOOL, numPersons);
    GrB_Vector_new(&next1, GrB_BOOL, numPersons);
    GrB_Vector_new(&next2, GrB_BOOL, numPersons);
    GrB_Vector_new(&intersection1, GrB_BOOL, numPersons);
    GrB_Vector_new(&intersection2, GrB_BOOL, numPersons);
    
    // init

    GrB_Vector_setElement_BOOL(frontier1, true, p1);
    GrB_Vector_setElement_BOOL(frontier2, true, p2);
    GrB_Vector_dup(&seen1, frontier1);
    GrB_Vector_dup(&seen2, frontier2);

    int distance = 0;

    if (p1 == p2) {
        distance = 0;
    } else {

        GrB_Matrix HasCreatorTransposed, ReplyOfTransposed;
        if (comment_lower_limit == -1) {
            HasCreatorTransposed = NULL;
            ReplyOfTransposed = NULL;
        } else {
            GrB_Matrix_new(&HasCreatorTransposed, GrB_UINT64, numPersons, NumComments);
            GrB_transpose(HasCreatorTransposed, NULL, NULL, HasCreator, NULL);
            GrB_Matrix_new(&ReplyOfTransposed, GrB_UINT64, NumComments, NumComments);
            GrB_transpose(ReplyOfTransposed, NULL, NULL, ReplyOf, NULL);
        }

        for (GrB_Index level = 1; level < numPersons / 2 + 1; level++) {
            advance_wavefront(HasCreator, HasCreatorTransposed, ReplyOf, ReplyOfTransposed, Knows, frontier1, next1, seen1, numPersons, comment_lower_limit, NumComments);

            GrB_Index next1nvals;
            GrB_Vector_nvals(&next1nvals, next1);
            if (next1nvals == 0) {
                distance = -1;
                break;
            }

            GrB_Vector_eWiseMult_BinaryOp(intersection1, NULL, NULL, GrB_LAND, next1, frontier2, NULL);

            GrB_Index intersection1_nvals;
            GrB_Vector_nvals(&intersection1_nvals, intersection1);
            if (intersection1_nvals > 0) {
                distance = level * 2 - 1;
                break;
            }

//                GrB_vxm(next2, seen2, NULL, GxB_ANY_PAIR_BOOL, frontier2, A, GrB_DESC_RC);
            advance_wavefront(HasCreator, HasCreatorTransposed, ReplyOf, ReplyOfTransposed, Knows, frontier2, next2, seen2, numPersons, comment_lower_limit, NumComments);

            GrB_Vector_eWiseMult_BinaryOp(intersection2, NULL, NULL, GrB_LAND, next1, next2, NULL);

            GrB_Index intersection2_nvals;
            GrB_Vector_nvals(&intersection2_nvals, intersection2);
            if (intersection2_nvals > 0) {
                distance = level * 2;
                break;
            }

            GrB_Index next2nvals;
            GrB_Vector_nvals(&next2nvals, next2);
            if (next2nvals == 0) {
                distance = -1;
                break;
            }

            GrB_eWiseAdd_Vector_BinaryOp(seen1, NULL, NULL, GrB_LOR, seen1, next1, NULL);
            GrB_eWiseAdd_Vector_BinaryOp(seen2, NULL, NULL, GrB_LOR, seen2, next2, NULL);

            GrB_Vector_dup(&frontier1, next1);
            GrB_Vector_dup(&frontier2, next2);
        }
    }

    printf("Distance: %d\n", distance);

    // Cleanup
    LAGraph_finalize();

    return 0;
}
