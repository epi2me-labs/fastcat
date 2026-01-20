#ifndef FASTCAT_FC_DUST_H
#define FASTCAT_FC_DUST_H

#include <stddef.h>
#include <stdint.h>

double fc_dust_fraction(const uint8_t *seq, size_t len, int t, int w);

#endif
