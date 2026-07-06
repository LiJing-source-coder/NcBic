#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "struct.h"
#include <vector>
#include <memory>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_dec_float_100;

/* Return current checkpoint phase (0=none, 1=seeds done, 2=cluster partial) */
int checkpoint_probe_phase(void);

/* Load edge_list and rec_num from checkpoint (phase >= 1) */
void checkpoint_load_seeds(int &rec_num, Edge **&edge_list);

/* Load cluster loop state from checkpoint (phase >= 2) */
void checkpoint_load_cluster_state(int &current_i,
                                   std::vector<bool> &allincluster,
                                   std::vector<std::unique_ptr<Block>> &bb);

/* Save full checkpoint (overwrites existing). phase=1 or 2. */
void checkpoint_save(int phase,
                     Edge **edge_list,
                     int rec_num,
                     int current_i,
                     const std::vector<bool> &allincluster,
                     const std::vector<std::unique_ptr<Block>> &bb);

/* Remove checkpoint file after successful completion */
void checkpoint_remove(void);

/* True if -C was set to >0 */
bool checkpoint_enabled(void);

/* Dump blocks from checkpoint to FILE* and return number of blocks written */
int checkpoint_dump_blocks(FILE* fw);

#endif
