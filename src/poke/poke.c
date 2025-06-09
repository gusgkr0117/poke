#include <hd.h>

int main() {
    int res = 1;

    randombytes_init((unsigned char *)"some", (unsigned char *)"string", 128);

    printf("Running hd module unit tests\n");

    res = hd_chain_test();

    if (res == 0) {
        printf("hd module unit tests passed\n");
    } else {
        printf("hd module unit tests failed\n");
    }

    return res;
}