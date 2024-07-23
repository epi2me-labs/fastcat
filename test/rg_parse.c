#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../src/common.h"

typedef struct {
    char* runid;
    char* basecall_model;
    char* mod_model;
    char* barcode;
    char* suffix;
} TestCase;


int compare(char* str1, char* str2) {
    if (str1 == NULL && str2 == NULL) {
        return 0;
    }
    if (str1 == NULL || str2 == NULL) {
        return 1;
    }
    return strncmp(str1, str2, max(strlen(str1), strlen(str2)));
}


int main() {
    char *runid = "ef1af1ab8967cb20ca30dbeca93fd66592bf4619";
    char *basecall_model = "basecall_model_name@v1.2.3";
    char *mod_model_name = "basecall_model_name@v1.2.3_5mCG_5hmCG@v1";
    char *barcode = "barcode01";
    char *suffix = "-1A2B3C4D";

    TestCase cases[] = {
        {runid, basecall_model, mod_model_name, barcode, suffix},
        {runid, basecall_model, mod_model_name, barcode, NULL},
        {runid, basecall_model, mod_model_name, NULL, suffix},
        {runid, basecall_model, mod_model_name, NULL, NULL},
        {runid, basecall_model, NULL, barcode, suffix},
        {runid, basecall_model, NULL, barcode, NULL},
        {runid, basecall_model, NULL, NULL, suffix},
        {runid, basecall_model, NULL, NULL, NULL},
    };

    int fails = 0;
    for (int i = 0; i < 8; i++) {

        char* read_group = calloc(400, sizeof(char));
        read_group = strcpy(read_group, cases[i].runid);
        if (cases[i].basecall_model != NULL) {
            read_group = strcat(read_group, "_");
            read_group = strcat(read_group, cases[i].basecall_model);
        }
        if (cases[i].mod_model != NULL) {
            read_group = strcat(read_group, "_");
            read_group = strcat(read_group, cases[i].mod_model);
        }
        if (cases[i].barcode != NULL) {
            read_group = strcat(read_group, "_");
            read_group = strcat(read_group, cases[i].barcode);
        }
        if (cases[i].suffix != NULL) {
            read_group = strcat(read_group, cases[i].suffix);
        }
        printf("Test case %d: %s\n", i, read_group);

        readgroup* info = create_rg_info(read_group);

        int fail = 0;
        fail += compare(info->runid, cases[i].runid) != 0;
        fail += compare(info->basecaller, cases[i].basecall_model) != 0;
        fail += compare(info->modcaller, cases[i].mod_model) != 0;
        fail += compare(info->barcode, cases[i].barcode) != 0;

        if (fail) {
            fails++;
            printf("  Failed\n");
            printf("    Expected: %s %s %s %s\n", cases[i].runid, cases[i].basecall_model, cases[i].mod_model, cases[i].barcode);
            printf("         Got: %s %s %s %s\n", info->runid, info->basecaller, info->modcaller, info->barcode);
        }

        free(info);
        free(read_group);
        printf("\n");
    }

    if (fails == 0) {
        printf("All tests passed\n");
    } else {
        printf("%d tests failed\n", fails);
    }
    return fails != 0;
}
