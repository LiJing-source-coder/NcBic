#include "checkpoint.h"
#include "write_block.h"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <string>

/* On-disk checkpoint layout (binary, little-endian assumed on x86_64) */
#pragma pack(push, 1)
struct CkptHeader {
    char     magic[8];          /* "NCBICKPT" */
    uint32_t version;           /* 1 */
    uint32_t phase;             /* 1 or 2 */
    char     fn[LABEL_LEN];
    int      rows;
    int      cols;
    int      col_width;
    int      row_maxwidth;
    int      row_minwidth;
    int      row_maxinit;
    double   r_threshold;
    double   seed_repvalue;
    double   tolerance;
    double   pvalue;
    double   filter;
    int      sch_block;
    int      rpt_block;
    uint64_t rec_num;
    uint64_t current_i;
    uint64_t bb_size;
};
#pragma pack(pop)

static const char CKPT_MAGIC[8] = {'N','C','B','I','C','K','P','T'};

static const char* checkpoint_filename(void)
{
    static char* cached = NULL;
    if (!cached)
        cached = addSuffix(po->FN, ".ckpt");
    return cached;
}

static void fill_header(CkptHeader* h, int phase,
                        int rec_num, int current_i,
                        std::size_t bb_size)
{
    memset(h, 0, sizeof(*h));
    memcpy(h->magic, CKPT_MAGIC, 8);
    h->version     = 1;
    h->phase       = phase;
    strncpy(h->fn, po->FN, LABEL_LEN-1);
    h->rows        = ::rows;
    h->cols        = ::cols;
    h->col_width   = po->COL_WIDTH;
    h->row_maxwidth= po->ROW_MAXWIDTH;
    h->row_minwidth= po->ROW_MINWIDTH;
    h->row_maxinit = po->ROW_MAXINIT;
    h->r_threshold = po->R_THRESHOLD;
    h->seed_repvalue= po->SEED_RePValue;
    h->tolerance   = po->TOLERANCE;
    h->pvalue      = po->pvalue;
    h->filter      = po->FILTER;
    h->sch_block   = po->SCH_BLOCK;
    h->rpt_block   = po->RPT_BLOCK;
    h->rec_num     = static_cast<uint64_t>(rec_num);
    h->current_i   = static_cast<uint64_t>(current_i);
    h->bb_size     = static_cast<uint64_t>(bb_size);
}

static bool same_param(const CkptHeader* h)
{
    if (strcmp(h->fn, po->FN) != 0)
        return false;
    if (h->col_width    != po->COL_WIDTH)    return false;
    if (h->row_maxwidth != po->ROW_MAXWIDTH) return false;
    if (h->row_minwidth != po->ROW_MINWIDTH) return false;
    if (h->row_maxinit  != po->ROW_MAXINIT)  return false;
    if (fabs(h->r_threshold - po->R_THRESHOLD) > 1e-9) return false;
    if (fabs(h->seed_repvalue - po->SEED_RePValue) > 1e-9) return false;
    if (fabs(h->tolerance - po->TOLERANCE) > 1e-9) return false;
    if (fabs(h->pvalue - po->pvalue) > 1e-9) return false;
    if (fabs(h->filter - po->FILTER) > 1e-9) return false;
    if (h->sch_block != po->SCH_BLOCK) return false;
    if (h->rpt_block != po->RPT_BLOCK) return false;
    return true;
}

/* ------------------------------------------------------------------ */

int checkpoint_probe_phase(void)
{
    const char* ckpt_fn = checkpoint_filename();
    FILE* fp = fopen(ckpt_fn, "rb");
    if (!fp) return 0;

    CkptHeader h;
    if (fread(&h, sizeof(h), 1, fp) != 1) {
        fclose(fp);
        return 0;
    }
    fclose(fp);

    if (memcmp(h.magic, CKPT_MAGIC, 8) != 0 || h.version != 1)
        return 0;

    if (!same_param(&h)) {
        errAbort("Checkpoint parameters mismatch; remove %s to start fresh",
                 ckpt_fn);
    }
    return static_cast<int>(h.phase);
}

void checkpoint_load_seeds(int &rec_num, Edge **&edge_list)
{
    const char* ckpt_fn = checkpoint_filename();
    FILE* fp = mustOpen(ckpt_fn, "rb");

    CkptHeader h;
    if (fread(&h, sizeof(h), 1, fp) != 1)
        errAbort("Failed to read checkpoint header");

    if (h.phase < 1)
        errAbort("Checkpoint phase < 1, cannot resume seeds");

    if (h.rows != ::rows || h.cols != ::cols)
        errAbort("Checkpoint rows/cols mismatch (%d,%d) vs current (%d,%d)",
                 h.rows, h.cols, ::rows, ::cols);

    rec_num = static_cast<int>(h.rec_num);
    edge_list = (Edge **)malloc(sizeof(Edge *) * rec_num);
    if (!edge_list) errAbort("Out of memory loading checkpoint seeds");

    for (int i = 0; i < rec_num; i++) {
        edge_list[i] = (Edge *)malloc(sizeof(Edge));
        if (fread(edge_list[i], sizeof(Edge), 1, fp) != 1)
            errAbort("Failed to read checkpoint edge %d", i);
    }
    fclose(fp);
    progress("Checkpoint loaded: %d seeds restored", rec_num);
}

void checkpoint_load_cluster_state(int &current_i,
                                   std::vector<bool> &allincluster,
                                   std::vector<std::unique_ptr<Block>> &bb)
{
    const char* ckpt_fn = checkpoint_filename();
    FILE* fp = mustOpen(ckpt_fn, "rb");

    CkptHeader h;
    if (fread(&h, sizeof(h), 1, fp) != 1)
        errAbort("Failed to read checkpoint header");

    if (h.phase < 2)
        errAbort("Checkpoint phase < 2, cannot resume clustering");

    if (h.rows != ::rows || h.cols != ::cols)
        errAbort("Checkpoint rows/cols mismatch (%d,%d) vs current (%d,%d)",
                 h.rows, h.cols, ::rows, ::cols);

    /* skip edge_list section */
    long edge_section = static_cast<long>(h.rec_num * sizeof(Edge));
    if (fseek(fp, edge_section, SEEK_CUR) != 0)
        errAbort("Failed to seek past checkpoint seeds");

    current_i = static_cast<int>(h.current_i);

    allincluster.resize(Mrows);
    for (int j = 0; j < Mrows; j++) {
        char c;
        if (fread(&c, 1, 1, fp) != 1)
            errAbort("Failed to read checkpoint allincluster");
        allincluster[j] = (c != 0);
    }

    bb.clear();
    for (std::size_t bidx = 0; bidx < h.bb_size; bidx++) {
        std::unique_ptr<Block> b(new Block());
        for (int c = 0; c < cols; c++) {
            int v;
            if (fread(&v, sizeof(int), 1, fp) != 1)
                errAbort("Failed to read checkpoint block BTP");
            b->BTP[c] = v;
        }
        for (int r = 0; r < Mrows; r++) {
            char c;
            if (fread(&c, 1, 1, fp) != 1)
                errAbort("Failed to read checkpoint block genes");
            b->genes[r] = (c != 0);
        }
        uint32_t plen = 0;
        if (fread(&plen, sizeof(plen), 1, fp) != 1)
            errAbort("Failed to read checkpoint p-value length");
        std::string pstr(plen, '\0');
        if (plen > 0) {
            if (fread(&pstr[0], 1, plen, fp) != plen)
                errAbort("Failed to read checkpoint p-value string");
        }
        b->p_value = cpp_dec_float_100(pstr);
        if (fread(&b->size, sizeof(int), 1, fp) != 1)
            errAbort("Failed to read checkpoint block size");
        bb.push_back(std::move(b));
    }

    fclose(fp);
    progress("Checkpoint loaded: resume from seed %d/%lu, %lu blocks recovered",
             current_i, static_cast<unsigned long>(h.rec_num),
             static_cast<unsigned long>(bb.size()));
}

void checkpoint_save(int phase,
                     Edge **edge_list,
                     int rec_num,
                     int current_i,
                     const std::vector<bool> &allincluster,
                     const std::vector<std::unique_ptr<Block>> &bb)
{
    const char* ckpt_fn = checkpoint_filename();
    char* tmp_fn = addSuffix(ckpt_fn, ".tmp");
    FILE* fp = fopen(tmp_fn, "wb");
    if (!fp) {
        err("Cannot open %s for checkpoint", tmp_fn);
        free(tmp_fn);
        return;
    }

    CkptHeader h;
    fill_header(&h, phase, rec_num, current_i, bb.size());

    if (fwrite(&h, sizeof(h), 1, fp) != 1)
        errAbort("Failed to write checkpoint header");

    /* always write edge_list so phase-2 files are self-contained */
    if (edge_list && rec_num > 0) {
        for (int i = 0; i < rec_num; i++) {
            if (fwrite(edge_list[i], sizeof(Edge), 1, fp) != 1)
                errAbort("Failed to write checkpoint edge %d", i);
        }
    }

    if (phase >= 2) {
        for (int j = 0; j < Mrows; j++) {
            char c = allincluster[j] ? 1 : 0;
            if (fwrite(&c, 1, 1, fp) != 1)
                errAbort("Failed to write checkpoint allincluster");
        }
        for (const auto& bptr : bb) {
            const Block* b = bptr.get();
            for (int c = 0; c < cols; c++) {
                int v = b->BTP[c];
                if (fwrite(&v, sizeof(int), 1, fp) != 1)
                    errAbort("Failed to write checkpoint BTP");
            }
            for (int r = 0; r < Mrows; r++) {
                char c = b->genes[r] ? 1 : 0;
                if (fwrite(&c, 1, 1, fp) != 1)
                    errAbort("Failed to write checkpoint genes");
            }
            std::string pstr = b->p_value.convert_to<std::string>();
            uint32_t plen = static_cast<uint32_t>(pstr.length());
            if (fwrite(&plen, sizeof(plen), 1, fp) != 1)
                errAbort("Failed to write checkpoint p-value length");
            if (plen > 0) {
                if (fwrite(pstr.c_str(), 1, plen, fp) != plen)
                    errAbort("Failed to write checkpoint p-value string");
            }
            if (fwrite(&b->size, sizeof(int), 1, fp) != 1)
                errAbort("Failed to write checkpoint block size");
        }
    }

    fflush(fp);
    fclose(fp);

    if (rename(tmp_fn, ckpt_fn) != 0)
        errAbort("Failed to rename checkpoint temp file to %s", ckpt_fn);

    free(tmp_fn);
    progress("Checkpoint saved to %s (phase %d)", ckpt_fn, phase);
}

void checkpoint_remove(void)
{
    const char* ckpt_fn = checkpoint_filename();
    remove(ckpt_fn);
}

bool checkpoint_enabled(void)
{
    return po->CHECKPOINT_INTERVAL > 0;
}

int checkpoint_dump_blocks(FILE* fw)
{
    const char* ckpt_fn = checkpoint_filename();
    FILE* fp = fopen(ckpt_fn, "rb");
    if (!fp) {
        err("Cannot open checkpoint file %s", ckpt_fn);
        return 0;
    }

    CkptHeader h;
    if (fread(&h, sizeof(h), 1, fp) != 1) {
        err("Failed to read checkpoint header from %s", ckpt_fn);
        fclose(fp);
        return 0;
    }

    if (memcmp(h.magic, CKPT_MAGIC, 8) != 0 || h.version != 1) {
        err("Invalid checkpoint file %s", ckpt_fn);
        fclose(fp);
        return 0;
    }

    if (h.phase < 2) {
        err("Checkpoint phase < 2, no blocks to dump in %s", ckpt_fn);
        fclose(fp);
        return 0;
    }

    if (h.rows != ::rows || h.cols != ::cols) {
        err("Checkpoint rows/cols mismatch (%d,%d) vs current (%d,%d)",
            h.rows, h.cols, ::rows, ::cols);
        fclose(fp);
        return 0;
    }

    /* skip edge_list section */
    long edge_section = static_cast<long>(h.rec_num * sizeof(Edge));
    if (fseek(fp, edge_section, SEEK_CUR) != 0) {
        err("Failed to seek past checkpoint seeds");
        fclose(fp);
        return 0;
    }

    /* skip allincluster */
    if (fseek(fp, Mrows, SEEK_CUR) != 0) {
        err("Failed to seek past checkpoint allincluster");
        fclose(fp);
        return 0;
    }

    int blocks_written = 0;
    for (std::size_t bidx = 0; bidx < h.bb_size; bidx++) {
        std::unique_ptr<Block> b(new Block());
        for (int c = 0; c < cols; c++) {
            int v;
            if (fread(&v, sizeof(int), 1, fp) != 1) {
                err("Failed to read checkpoint block BTP");
                fclose(fp);
                return blocks_written;
            }
            b->BTP[c] = v;
        }
        for (int r = 0; r < Mrows; r++) {
            char c;
            if (fread(&c, 1, 1, fp) != 1) {
                err("Failed to read checkpoint block genes");
                fclose(fp);
                return blocks_written;
            }
            b->genes[r] = (c != 0);
        }
        uint32_t plen = 0;
        if (fread(&plen, sizeof(plen), 1, fp) != 1) {
            err("Failed to read checkpoint p-value length");
            fclose(fp);
            return blocks_written;
        }
        std::string pstr(plen, '\0');
        if (plen > 0) {
            if (fread(&pstr[0], 1, plen, fp) != plen) {
                err("Failed to read checkpoint p-value string");
                fclose(fp);
                return blocks_written;
            }
        }
        b->p_value = cpp_dec_float_100(pstr);
        if (fread(&b->size, sizeof(int), 1, fp) != 1) {
            err("Failed to read checkpoint block size");
            fclose(fp);
            return blocks_written;
        }
        print_bc(fw, b, static_cast<int>(bidx));
        blocks_written++;
    }

    fclose(fp);
    progress("Dumped %d blocks from %s", blocks_written, ckpt_fn);
    return blocks_written;
}
