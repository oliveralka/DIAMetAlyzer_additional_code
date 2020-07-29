#!/bin/bash

source activate py37

# TSV

# rtperm
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_rtperm.tsv -m rtperm

# fragdb 
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_fragdb.tsv -m fragdb

# fragdb with precursor CH2 shift 
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_fragdb_ch2prec.tsv -m fragdb -ch2prec

# shufflefrag
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_shufflefrag.tsv -m shufflefrag 

# shufflefrag with precursor CH2 shift
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_shufflefrag_ch2prec.tsv -m shufflefrag -ch2prec

# sw_perm
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_swperm.tsv -m sw_perm -sw swathwin.txt

# linmzperm 
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_linmzperm.tsv -m linmzperm 

# linmzperm with precursor CH2 shift
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_linmzperm_ch2prec.tsv -m linmzperm -ch2prec

#PQP

# rtperm
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_rtperm.pqp -m rtperm

# fragdb 
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_fragdb.pqp -m fragdb

# fragdb with precursor CH2 shift 
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_fragdb_ch2prec.pqp -m fragdb -ch2prec

# shufflefrag
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_shufflefrag.pqp -m shufflefrag 

# shufflefrag with precursor CH2 shift
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_shufflefrag_ch2prec.pqp -m shufflefrag -ch2prec

# sw_perm
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_swperm.pqp -m sw_perm -sw swathwin.txt

# linmzperm 
python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_linmzperm.pqp -m linmzperm 

# linmzperm with precursor CH2 shift
#python ./DecoyGeneratorMetaboTool_2.0.py -as ./20190927_reference_lib_woprec.tsv -out 20190927_tdlib_woprec_linmzperm_ch2prec.pqp -m linmzperm -ch2prec
