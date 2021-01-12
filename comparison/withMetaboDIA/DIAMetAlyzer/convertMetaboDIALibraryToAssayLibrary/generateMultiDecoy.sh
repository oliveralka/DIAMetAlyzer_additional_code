#!/bin/bash

# TSV

# rtperm
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20200930_tdlib_library_pos_rtperm_250s.tsv -m rtperm --rt_mindist 250 --retentiontime 1200

# fragdb 
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_fragdb_test.tsv -m fragdb

# shufflefrag
##python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_shufflefrag.tsv -m shufflefrag 

# sw_perm
##python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_swperm.tsv -m sw_perm -sw swathwin.txt

# linmzperm 
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_linmzperm.tsv -m linmzperm 



#PQP

# rtperm
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_rtperm_250s.pqp -m rtperm --rt_mindist 250 --retentiontime 1200

# fragdb 
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20201205_DDA_msconvert_mzXML_10ppm_min_s00001_pf00.tsv -out 20201205_DDA_tdlib_msconvert_mzXML_10ppm_min_s00001_pf00.pqp -m fragdb
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20201205_DDA_msconvert_mzXML_10ppm_min_s02_pf00.tsv -out 20201205_DDA_tdlib_msconvert_mzXML_10ppm_min_s02_pf00.pqp -m fragdb
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20201205_DDA_msconvert_mzXML_10ppm_min_s05_pf00.tsv -out 20201205_DDA_tdlib_msconvert_mzXML_10ppm_min_s05_pf00.pqp -m fragdb

# shufflefrag
##python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_shufflefrag.pqp -m shufflefrag 

# sw_perm
##python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_swperm.pqp -m sw_perm -sw swathwin.txt

# linmzperm 
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20200929_assay_library_pos.tsv -out 20190927_tdlib_woprec_linmzperm.pqp -m linmzperm 

