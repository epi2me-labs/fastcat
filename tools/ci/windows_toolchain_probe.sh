#!/usr/bin/env bash

set -euo pipefail

export PATH=/ucrt64/bin:/usr/bin:/bin:"$PATH"
cc=/ucrt64/bin/gcc
probe_dir=build/windows-toolchain-probe

mkdir -p "$probe_dir"

cat > "$probe_dir/printf.c" <<'EOF'
#include <stddef.h>
#include <stdio.h>

int main(void)
{
    printf("%zu\n", (size_t)1);
}
EOF

cat > "$probe_dir/fc_argp.c" <<'EOF'
#include "core/fc_argp.h"
#include <stddef.h>

int main(void)
{
    argp_error(NULL, "count: %zu", (size_t)1);
}
EOF

cat > "$probe_dir/regex.c" <<'EOF'
#include <regex.h>

int main(void)
{
    regex_t regex;
    return regcomp(&regex, "fastcat", REG_EXTENDED);
}
EOF

echo "MSYSTEM=${MSYSTEM:-unset}"
echo "CC=$cc"
echo "cc=$(command -v cc)"
cc -dumpmachine
echo "gcc=$(command -v gcc)"
gcc -dumpmachine
echo "ucrt=$($cc -print-file-name=libucrt.a)"
pacman -Q mingw-w64-ucrt-x86_64-libsystre

"$cc" -Wall -Wextra -Werror "$probe_dir/printf.c" -o "$probe_dir/printf.exe"
"$probe_dir/printf.exe"
echo "ordinary percent-zu: passed"

"$cc" -Isrc -Ithird_party/argp-standalone -Wall -Wextra -Werror \
    -c "$probe_dir/fc_argp.c" -o "$probe_dir/fc_argp.o"
echo "fastcat argp percent-zu: passed"

"$cc" -Wall -Wextra -Werror -static "$probe_dir/regex.c" \
    -lregex -ltre -lintl -liconv -o "$probe_dir/regex.exe"
echo "HTSlib POSIX regex support: passed"
