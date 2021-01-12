#!/bin/bash

#SBATCH --job-name=diaumpire_metabodia_data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 24
#SBATCH --time=48:00:00
#SBATCH --mem=192G
#SBATCH --output=/nfs/wsi/abi/scratch/alka/DIAUmpire/20201207_diaumpire_mtbls_ourparam_mdia.stdout
#SBATCH --error=/nfs/wsi/abi/scratch/alka/DIAUmpire/20201207_diaumpire_mtbls_ourparam_mdia.stderr

# Careful number of threads are set in the params files of DIAUmpire!  

echo "DIAUmpire analyis"

mzXML_dir="/nfs/wsi/abi/scratch/alka/DIAUmpire/data_mtbls_417_ourparam"
params_file="/nfs/wsi/abi/scratch/alka/DIAUmpire/params/params_diaumpire_mtbls417_ourparam.se_params"
umpire_jar="/home/alka/software/DIAUmpire/v2.1.6/v2.1.6/DIA_Umpire_SE.jar"
java_executable="/home/alka/software/jre1.8.0_121/bin/java"

for filename in $mzXML_dir/*.mzXML; do
  echo $filename
  echo $params_file

  $java_executable -jar -Xmx8G $umpire_jar $filename $params_file

done

echo "All done. Have a nice one!" >&2

