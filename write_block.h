#ifndef _WRITE_BLOCK_H
#define _WRITE_BLOCK_H

#include "struct.h"
#include <memory>

/* prototypes */
static bool condcount(int a);
void print_bc(FILE *fw, const std::unique_ptr<Block> &b,int num);
#endif

