#include "ec.h"
#include "isog.h"
#include <assert.h>

void
ec_eval_five(ec_curve_t *image,
              const ec_isog_odd_t *phi,
              ec_point_t *points,
              unsigned short length)
{

    ec_point_t R, A24, B24, A3;
    ec_point_t stack[POWER_OF_5];
    unsigned int stack_index[POWER_OF_5];
    int index = 0, nstack = 0, ii = 0;
    digit_t m;

    copy_point(&R, &phi->ker_minus);

    AC_to_A24(&A24, &phi->curve);
    fp2_sub(&A3.z, &A24.x, &A24.z);
    fp2_copy(&A3.x, &A24.x);

    // Traverse Tree
    for (int row = 1; row < phi->degree[1]; row ++) {
        while (index < phi->degree[1] - row) {
            copy_point(&stack[nstack], &R);
            stack_index[nstack++] = index;
            m = STRATEGY5[ii++];
            for (int j = 0; j < m; j++) {
                xMUL_FIVE(&R, &R, &A3, &A24);
            }
            index += m;
        }
        kps(1, R, A24);
        xisog(&B24, 1, A24);
        copy_point(&A24, &B24);
        fp2_sub(&A3.z, &A24.x, &A24.z);
        fp2_copy(&A3.x, &A24.x);

        for(int i = 0; i < nstack; i++) {
            xeval(&stack[i], 1, stack[i], A24);
        }

        for(int i = 0; i < length; i++) {
            xeval(&points[i], 1, points[i], A24);
        }

        copy_point(&R, &stack[nstack - 1]);
        index = stack_index[nstack - 1];
        nstack--;
    }
    kps(1, R, A24);
    xisog(&B24, 1, A24);
    copy_point(&A24, &B24);

    for(int i = 0; i < length; i++) {
        xeval(&points[i], 1, points[i], A24);
    }

    A24_to_AC(image, &A24);
    image->is_A24_computed_and_normalized = 0;
}