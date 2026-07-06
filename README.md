# NcBic
NcBic detects local non-constant bicluster patterns for complex disease analysis, with demonstrated applications in gene expression and MRI volumetric data such as ADNI Alzheimer's studies.

Installation
Before installation, please make sure that g++ (GCC) supporting C++11, OpenMP, and Boost libraries have already been installed on your computer.
Extract the package
	% unzip NcBic.zip
Move to the NcBic folder
	% cd NcBic
Make an executable programe
	% make
Then you can use it with following options.

+
===================================================================
[Usage]
$ ./TransBic -i filename -o 1 
$ ./TransBic -i filename [argument list] 
===================================================================
[Input]
-i : input file must be the following tab-delimited format
     data 
     -------------------------------------
     o        cond1    cond2    cond3
     gene1     2.4      3.5     -2.4
     gene2    -2.1      0.0      1.2
     -------------------------------------
-x : set the max number of rows of biclusters may have(data type int)
	 default: 0.5*size of the input matrix
-z : set the min number of rows of biclusters may have(data type int)
	 default: 5
-g : set the total p-value of removing proportion of seeds (data type double)
	 default: 0.6827
-t : set the max number of rows to initial  the common trend (data type int)
     default: 60
-r : set min value of the similarity ratio r for forming a new BTP-pattern (data type double)
     default: 0.8
===================================================================
[Output]
-o : number of blocks to report, default: 1
-k : minimum column width of the block,
     default: 5% of columns, minimum 2 columns
-f : filtering overlapping blocks,
     default: 1 (do not remove any blocks)
-c : consistency level of the block (0.5-1.0].
     default: 0.95
-p : p-value threshold of the output block (0-0.05], default: 0.001
-C : checkpoint interval in minutes (data type int),
     default: 0 (disabled)
-D : dump blocks from checkpoint file and exit
===================================================================
