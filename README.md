# NcBic

NcBic detects local non-constant bicluster patterns for complex disease analysis, with demonstrated applications in gene expression and MRI volumetric data such as ADNI Alzheimer's studies.

## Installation

Before installation, please make sure that **g++ (GCC)** supporting **C++11**, **OpenMP**, and **Boost** libraries have already been installed on your computer.

**Extract the package**

```bash
unzip NcBic.zip
```

**Move to the NcBic folder**

```bash
cd NcBic
```

**Make an executable program**

```bash
make
```

Then you can use it with the following options.

## Usage

```bash
./TransBic -i filename -o 1
./TransBic -i filename [argument list]
```

## Input

`-i` : input file must be the following tab-delimited format

```
o        cond1    cond2    cond3
gene1     2.4      3.5     -2.4
gene2    -2.1      0.0      1.2
```

| Option | Description | Default |
|--------|-------------|---------|
| `-x` | set the max number of rows of biclusters may have (`int`) | 0.5 × size of the input matrix |
| `-z` | set the min number of rows of biclusters may have (`int`) | 5 |
| `-g` | set the total p-value of removing proportion of seeds (`double`) | 0.6827 |
| `-t` | set the max number of rows to initial the common trend (`int`) | 60 |
| `-r` | set min value of the similarity ratio r for forming a new BTP-pattern (`double`) | 0.8 |

## Output

| Option | Description | Default |
|--------|-------------|---------|
| `-o` | number of blocks to report | 1 |
| `-k` | minimum column width of the block | 5% of columns, minimum 2 columns |
| `-f` | filtering overlapping blocks | 1 (do not remove any blocks) |
| `-c` | consistency level of the block (0.5–1.0] | 0.95 |
| `-p` | p-value threshold of the output block (0–0.05] | 0.001 |
| `-C` | checkpoint interval in minutes (`int`) | 0 (disabled) |
| `-D` | dump blocks from checkpoint file and exit | — |
