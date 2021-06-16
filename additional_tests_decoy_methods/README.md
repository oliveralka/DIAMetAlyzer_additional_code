Collection of scripts and tools to generate a target-decoy assay library.

Either direcly from the assay library (DecoyGeneratorMetaboTool_2) or by using the spectra information based on SIRIUS fragmentation information and passatutto in a KNIME worklfow or on the command line (KNIME).

The DecoyGeneratorMetaboTool_2.0 is able to generate decoys directly from the assay library. 

Here, different methods can be applied: 
1: rtperm - retention time permutation (same precursor and fragments - other RT) 
2: mz (same precursor and fragments - other mz or other window))
3: linmassperm - same precursor mz, but with shifted by mass difference fragments
4: swathwindow - precursor mz to other mass window
5: fragdb - use the same precursor, but fragments with lower mz as prec from all the fragment masses existing in the assay annotatet library.

