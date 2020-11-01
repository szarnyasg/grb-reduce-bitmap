#define LAGRAPH_EXPERIMENTAL_ASK_BEFORE_BENCHMARKING 1

#include <GraphBLAS.h>
#include <LAGraph.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define LAGRAPH_FREE_ALL            \
{                                   \
    /* TODO */                      \
}

// "!!" means that the operation is expected to be expensive

void advance_wavefront(GrB_Matrix HasCreator, GrB_Matrix ReplyOf, GrB_Matrix Knows, GrB_Vector frontier, GrB_Vector next, GrB_Vector seen, GrB_Index numPersons, GrB_Index numComments, int64_t comment_lower_limit) {
    if (comment_lower_limit == -1) {
        GrB_vxm(next, seen, NULL, GxB_ANY_PAIR_BOOL, frontier, Knows, GrB_DESC_RC);
    } else {
        GxB_Scalar limit;
        GxB_Scalar_new(&limit, GrB_INT64);
        GxB_Scalar_setElement_INT64(limit, comment_lower_limit);

        // build selection matrix based on the frontier's content

        // turn vector "frontier" into a selection matrix "FrontierSel"
        GrB_Index numFrontierPersons;
        GrB_Vector_nvals(&numFrontierPersons, frontier);
        GrB_Index *I = (GrB_Index*) (LAGraph_malloc(numFrontierPersons, sizeof(GrB_Index)));
        bool *X = (bool*) (LAGraph_malloc(numFrontierPersons, sizeof(bool)));
        GrB_Index nvalsTmp = numFrontierPersons;
        GrB_Vector_extractTuples(I, X, &nvalsTmp, frontier);
        assert(numFrontierPersons == nvalsTmp);

        printf("-------------- start of performance bottleneck --------------\n");
        GrB_Matrix CommentsOfFrontierPeople;
        GrB_Matrix_new(&CommentsOfFrontierPeople, GrB_BOOL, numPersons, numComments);

        // <option0>: use selection matrix and mxm to extract row
        printf(">>>>>>>>>>>>>>>>>>>> option 0\n");
        GrB_Matrix FrontierSel0;
        GrB_Matrix_new(&FrontierSel0, GrB_BOOL, numPersons, numPersons);
        GrB_Matrix_build(FrontierSel0, I, I, X, numFrontierPersons, GrB_LOR);
        printf(" !!");
        GrB_Matrix CommentsOfFrontierPeopleT;
        GrB_Matrix_new(&CommentsOfFrontierPeopleT, GrB_BOOL, numComments, numPersons);
        GrB_mxm(CommentsOfFrontierPeopleT, NULL, NULL, GxB_ANY_PAIR_BOOL, HasCreator, FrontierSel0, NULL);
        GrB_transpose(CommentsOfFrontierPeople, NULL, NULL, CommentsOfFrontierPeopleT, NULL);
        printf("\n");
        // </option0>

        // <option1>: use selection matrix and mxm to extract row
        printf(">>>>>>>>>>>>>>>>>>>> option 1\n");
        GrB_Matrix FrontierSel1;
        GrB_Matrix_new(&FrontierSel1, GrB_BOOL, numPersons, numPersons);
        GrB_Matrix_build(FrontierSel1, I, I, X, numFrontierPersons, GrB_LOR);
        printf(" !!");
        GrB_mxm(CommentsOfFrontierPeople, NULL, NULL, GxB_ANY_PAIR_BOOL, FrontierSel1, HasCreator, GrB_DESC_T1);
        // </option1>

        // <option2>: use extract on the transposed variant of HasCreator then assign the extracted matrix
        GrB_Matrix CommentsOfFrontierPeopleExtracted2;
        printf(">>>>>>>>>>>>>>>>>>>> option 2\n");
        GrB_Matrix_new(&CommentsOfFrontierPeopleExtracted2, GrB_BOOL, numFrontierPersons, numComments);
        printf(" !!");
        GrB_Matrix_extract(CommentsOfFrontierPeopleExtracted2, NULL, NULL, HasCreator, I, numFrontierPersons, GrB_ALL, numComments, GrB_DESC_T0);
        GrB_Matrix_assign(CommentsOfFrontierPeople, NULL, NULL, CommentsOfFrontierPeopleExtracted2, I, numFrontierPersons, GrB_ALL, numComments, NULL);
        // </option2>

        // <option3>: use extract on HasCreator then assign the transposed variant of the extracted matrix
        printf(">>>>>>>>>>>>>>>>>>>> option 3\n");
        GrB_Matrix CommentsOfFrontierPeopleExtracted3;
        GrB_Matrix_new(&CommentsOfFrontierPeopleExtracted3, GrB_BOOL, numComments, numFrontierPersons);
        printf(" !!");
        GrB_Matrix_extract(CommentsOfFrontierPeopleExtracted3, NULL, NULL, HasCreator, GrB_ALL, numComments, I, numFrontierPersons, NULL);
        GrB_Matrix_assign(CommentsOfFrontierPeople, NULL, NULL, CommentsOfFrontierPeopleExtracted3, I, numFrontierPersons, GrB_ALL, numComments, GrB_DESC_T0);
        // </option3>

        // <option4>: precompute HasCreator^T, then do extract/assign
        printf(">>>>>>>>>>>>>>>>>>>> option 4\n");
        GrB_Matrix HasCreatorT1;
        GrB_Matrix_new(&HasCreatorT1, GrB_BOOL, numPersons, numComments);
        printf(" !!");
        GrB_transpose(HasCreatorT1, NULL, NULL, HasCreator, NULL);
        GrB_Matrix CommentsOfFrontierPeopleExtracted4;
        GrB_Matrix_new(&CommentsOfFrontierPeopleExtracted4, GrB_BOOL, numFrontierPersons, numComments);
        GrB_Matrix_extract(CommentsOfFrontierPeopleExtracted4, NULL, NULL, HasCreatorT1, I, numFrontierPersons, GrB_ALL, numComments, NULL);
        GrB_Matrix_assign(CommentsOfFrontierPeople, NULL, NULL, CommentsOfFrontierPeopleExtracted4, I, numFrontierPersons, GrB_ALL, numComments, NULL);
        // </option4>

        // <option5>: precompute HasCreator^T, then do mxm with the selection matrix
        printf(">>>>>>>>>>>>>>>>>>>> option 5\n");
        GrB_Matrix HasCreatorT2;
        GrB_Matrix_new(&HasCreatorT2, GrB_BOOL, numPersons, numComments);
        printf(" !!");
        GrB_transpose(HasCreatorT2, NULL, NULL, HasCreator, NULL);
        GrB_Matrix FrontierSel2;
        GrB_Matrix_new(&FrontierSel2, GrB_BOOL, numPersons, numPersons);
        GrB_Matrix_build(FrontierSel2, I, I, X, numFrontierPersons, GrB_LOR);
        GrB_mxm(CommentsOfFrontierPeople, NULL, NULL, GxB_ANY_PAIR_BOOL, FrontierSel2, HasCreatorT2, NULL);
        // </option5>

        printf("--------------- end of performance bottleneck ---------------\n");

        // direction 1
        GrB_Matrix RepliesToCommentsOfFrontierPeople;
        GrB_Matrix_new(&RepliesToCommentsOfFrontierPeople, GrB_UINT64, numPersons, numComments);
        GrB_mxm(RepliesToCommentsOfFrontierPeople, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, CommentsOfFrontierPeople, ReplyOf, NULL);

        GrB_Matrix InteractionsToComments;
        GrB_Matrix_new(&InteractionsToComments, GrB_UINT64, numPersons, numPersons);
        GrB_mxm(InteractionsToComments, Knows, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, RepliesToCommentsOfFrontierPeople, HasCreator, NULL);

        // direction 2
        GrB_Matrix RepliesFromCommentsOfFrontierPeople;
        GrB_Matrix_new(&RepliesFromCommentsOfFrontierPeople, GrB_UINT64, numPersons, numComments);
        printf(" !!");
        GrB_mxm(RepliesFromCommentsOfFrontierPeople, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, CommentsOfFrontierPeople, ReplyOf, GrB_DESC_T1);

        GrB_Matrix InteractionsFromComments;
        GrB_Matrix_new(&InteractionsFromComments, GrB_UINT64, numPersons, numPersons);
        GrB_mxm(InteractionsFromComments, InteractionsToComments, GrB_NULL, GrB_PLUS_TIMES_SEMIRING_UINT64, RepliesFromCommentsOfFrontierPeople, HasCreator, NULL);

        // InteractionsToComments = InteractionsToComments * InteractionsFromComments
        GrB_eWiseMult(InteractionsToComments, InteractionsToComments, NULL, GrB_MIN_UINT64, InteractionsToComments, InteractionsFromComments, NULL);
        GxB_select(InteractionsToComments, NULL, NULL, GxB_GT_THUNK, InteractionsToComments, limit, NULL);
        GrB_reduce(next, NULL, NULL, GxB_PAIR_BOOL, InteractionsToComments, GrB_DESC_T0);
    }
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

    GrB_Matrix Knows, ReplyOf, HasCreator;
    LAGraph_binread(&Knows, "../knows.grb");
    LAGraph_binread(&ReplyOf, "../replyOf.grb");
    LAGraph_binread(&HasCreator, "../hasCreator.grb");

    GrB_Index numComments;
    GrB_Matrix_nrows(&numComments, ReplyOf);

    GrB_Index numPersons;
    GrB_Matrix_nrows(&numPersons, Knows);

    // hard-coded input params
    // try to select persons who have a path between them with the prescribed comment_lower_limit

    GrB_Index p1 = 31621;
    GrB_Index p2 = 79481;
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

    GrB_Vector_setElement(frontier1, true, p1);
    GrB_Vector_setElement(frontier2, true, p2);
    GrB_Vector_dup(&seen1, frontier1);
    GrB_Vector_dup(&seen2, frontier2);

    int distance = 0;

    // measurem processing time using LAGraph_tic/-toc
    double tic [2] ;
    LAGraph_tic (tic) ;

    if (p1 == p2) {
        distance = 0;
    } else {
        for (GrB_Index level = 1; level < numPersons / 2 + 1; level++) {
            // advance first wavefront
            advance_wavefront(HasCreator, ReplyOf, Knows, frontier1, next1, seen1, numPersons, numComments, comment_lower_limit);

            GrB_Index next1nvals;
            GrB_Vector_nvals(&next1nvals, next1);
            if (next1nvals == 0) {
                distance = -1;
                break;
            }

            GrB_eWiseMult(intersection1, NULL, NULL, GrB_LAND, next1, frontier2, NULL);

            GrB_Index intersection1_nvals;
            GrB_Vector_nvals(&intersection1_nvals, intersection1);
            if (intersection1_nvals > 0) {
                distance = level * 2 - 1;
                break;
            }

            // advance second wavefront
            advance_wavefront(HasCreator, ReplyOf, Knows, frontier2, next2, seen2, numPersons, numComments, comment_lower_limit);

            GrB_eWiseMult(intersection2, NULL, NULL, GrB_LAND, next1, next2, NULL);

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

            GrB_eWiseAdd(seen1, NULL, NULL, GrB_LOR, seen1, next1, NULL);
            GrB_eWiseAdd(seen2, NULL, NULL, GrB_LOR, seen2, next2, NULL);

            GrB_Vector_dup(&frontier1, next1);
            GrB_Vector_dup(&frontier2, next2);
        }
    }
    double elapsed = LAGraph_toc(tic);

    printf("Distance: %d\n", distance);
    printf("Processing time %12.3f sec\n", elapsed);

    // Cleanup
    LAGraph_finalize();

    return 0;
}
