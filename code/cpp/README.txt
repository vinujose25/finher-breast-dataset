============
Introduction
============

Fast C++ implementation of Wilcoxon test-based probeset detection calling algorithm for PM-probe (Perfect-Match-probe) only HG-U219 array using anti-genomic probes as a background. The anti-genomic probes were introduced as a cost-effective replacement for the traditional MM (Mis-Match) probes as a background.
Due to the lack of MM probes, the detection calling algorithm associated with MAS5 normalization will not work. Hence an updated algorithm is necessary for detection calling.

The R implementation of the Wilcoxon test performs poorly in looping algorithms.
Parallel implementation of the Wilcoxon test in R would be a solution, but it has not been attempted as the code is intended to run in a single-core processor. The current implementation supports only single-core computation.

C++ implementation of the Wilcoxn test available from the alglib library, version 3.8.2, was used (see https://www.alglib.net/).
The R scripts associated with this detection calling algorithm, "dabg2call.r" and "mas5gccall.r", are available in the folder, code/r/.

=========================================
Description of the file/folder/C++ script
=========================================

1. dabg2call.cpp
1.1. C++ implementation of detection calling algorithm based on “mannwhitneyutest” (Wilcoxon rank-sum test for two independent samples) functionfrom the alglib library. Note that the R script “dabg2call.R” depends on the shared object generated from this script (“dabg2call.cpp”).

2. dabg2call.h
2.1. Header file associated with dabg2call.cpp.

3. dabg2call.o
3.1. Object file generated from dabg2call.cpp.

4. dabg2call.so
4.1. Shared object file generated from dabg2call.cpp.

5. mas5gccall.cpp
5.1. C++ implementation of detection calling algorithm based on “wilcoxonsignedranktest” (Wilcoxon signed-rank test for two matched/paired samples) function from the alglib library. Note that the R script “mas5gccall.R” depends on the shared object generated from this script (“mas5gccall.cpp”).

6. mas5gccall.o
6.1. Object file generated from mas5gccall.cpp.

7. mas5gccall.so
7.1. Shared object file generated from mas5gccall.cpp.

8. alglib-3.8.2
8.1. Local copy of the alglib library, version 3.8.2, downloaded from https://www.alglib.net/

9. wilcox.test.default.r
9.1. The default R implementation of the Wilcoxon test is used as a reference.

10. getDABG2.RData
10.1. Data to test the detection-calling algorithm.
10.2. R-object descriptions:
10.2.1. antiGenomic - list of anti-genomic probesets in HG-U219 array.
10.2.2. bgCor.rma.u219.ff.sen - AffyBatch object of HG-U219 expression profiles.
10.2.3. hgu219.probeInfo.gcIndexList – List of, per probeset, probe GC content and probe index in intensity matrix. This object is required for the detection call algorithm to work.

