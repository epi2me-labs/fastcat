// bamstats program

#include <err.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>


enum command_mode {
    MODE_HELP = 0,
    MODE_BUILD,
    MODE_FETCH,
    MODE_DUMP,
    MODE_INVALID };
static const enum command_mode ncommand = MODE_INVALID;

enum command_mode get_mode(const char *modestr) {
    if (0 == strcmp(modestr, "help")) {return MODE_HELP;}
    if (0 == strcmp(modestr, "build")) {return MODE_BUILD;}
    if (0 == strcmp(modestr, "fetch")) {return MODE_FETCH;}
    if (0 == strcmp(modestr, "dump")) {return MODE_DUMP;}
    return MODE_INVALID;
}

void fprint_commands(void);

const char *mode_string(const enum command_mode mode) {
    switch (mode) {
    case MODE_HELP:
        return "help";
    case MODE_BUILD:
        return "build";
    case MODE_FETCH:
        return "fetch";
    case MODE_DUMP:
        return "dump";
    case MODE_INVALID:
        fprint_commands();
        errx(EXIT_FAILURE, "Invalid subcommand\n");
    default:
        errx(EXIT_FAILURE, "bamindex failure -- report bug\n");
    }

    return NULL;
}

const char *mode_description(const enum command_mode mode) {
    switch (mode) {
    case MODE_HELP:
        return "Print general help or help about a subcommand.";
    case MODE_BUILD:
        return "Build a BAM index.";
    case MODE_FETCH:
        return "Fetch records from a BAM using an index.";
    case MODE_DUMP:
        return "Dump an index fetch to text.";
    case MODE_INVALID:
        fprint_commands();
        errx(EXIT_FAILURE, "Invalid subcommand\n");
    default:
        errx(EXIT_FAILURE, "bamindex failure -- report bug\n");
    }

    return NULL;
}

void fprint_commands(void) {
    for (enum command_mode i = 0; i < ncommand; i++) {
        fprintf(
            stderr, "* bamindex %-14s%s\n", mode_string(i), mode_description(i));
    }
}
int main_build(int argc, char *argv[]);
int main_fetch(int argc, char *argv[]);
int main_dump(int argc, char *argv[]);

int main(int argc, char *argv[]) {

    if (argc == 1) {
        // Called as program name on it's own
        fprint_commands();
        return EXIT_SUCCESS;
    }

    int ret = EXIT_FAILURE;
    switch (get_mode(argv[1])) {
    case MODE_HELP:
        fprint_commands();
        break;
    case MODE_BUILD:
        ret = main_build(argc - 1, argv + 1);
        break;
    case MODE_FETCH:
        ret = main_fetch(argc - 1, argv + 1);
        break;
    case MODE_DUMP:
        ret = main_dump(argc - 1, argv + 1);
        break;
    default:
        ret = EXIT_FAILURE;
        warnx("Unrecognised subcommand %s\n", argv[1]);
    }

    return ret;
}

