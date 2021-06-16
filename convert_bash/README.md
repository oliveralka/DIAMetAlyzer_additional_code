This directory contains the scripts 
used for the conversion from .wiff to .mzML files.

Raw DIA files (.wiff) were converted to .mzML files without any additional processing (msconvert_only.sh)
using msconvert (Proteowizard).

Raw DDA files were centroided using the qtofpeakpicker and then converted to .mzML using msconvert (Proteowizard) (qtofpeakpicker_msconvert_bash.sh).

The test_bash.sh script shows how the script could be addapted to convert all files in a specific directoy. 
