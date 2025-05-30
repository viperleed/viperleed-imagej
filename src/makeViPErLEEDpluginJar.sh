#!/bin/sh
# Creates the ViPErLEED_ImageJ_plugins.jar file with the ImageJ LEED plugins.
# This file must be called in the plugins directory with the sources and class files.

tmpfile=/tmp/$$makeViPErLEEDplugin_extracted_byte_tmp
for i in *.class
do
    #major version number is the 7th byte; extract it
    dd status='none' skip=7 count=1 bs=1 if=$i of=$tmpfile
    #we want java6, major version number byte decimal 50, ascii '2'
    if [ `cat $tmpfile` != "2" ] ; then
        echo "WARNING: $i is NOT compiled for Java 6 (major version 50):"
    fi
done

for i in Open_Aida_LEED_Video.class Open_LEED_Movie.class LEED_Spot_Tracker.class LEED_Curve_Editor.class LEED_Data_Quality_Statistics.class LEED_R_Factor_Between_Datasets.class LEED_Stitch_Datasets.class LEED_Average_Datasets.class LEED_IV_Curve_Interpolation.class LEED_IV_Curve_Tools.class LEED_Commands.class
do
if [ ! -f $i ]; then
    echo "File $i not found"
    exit 1
fi
done

mkdir -p source/
cp -p *LEED*.java Leed*.java source/
zip --quiet ../ViPErLEED_.zip 00README.txt plugins.config Open_Aida_LEED_Video*.class Open_LEED_Movie*.class LEED_Spot_Tracker*.class LEED_Curve_Editor*.class LEED_R_Factor*.class LEED_Data_Quality_Statistics*.class LEED_R_Factor_Between_Datasets*.class LEED_IV_Curve_Tools*.class Leed*.class Averaging_LEED_Stack*.class Linear_Fit_LEED_Stack*.class LEED_Average_Datasets*.class LEED_IV_Curve_Interpolation*.class LEED_Stitch_Datasets*.class LEED_Commands*.class LICENSE*.txt source/*
rm source/*.java
mv ../ViPErLEED_.zip ../ViPErLEED_ImageJ_plugins.jar

rm -f $tmpfile



