# viperleed-imagej
ViPErLEED-ImageJ is a collection of ImageJ plugins for spot tracking and intensity evaluation as well as selection and editing of I(V) curves created by the Spot Tracker.
The output of the I(V) Curve Editor can be used as experimental data (`EXPBEAMS.csv`) for structure optimization.

## Installation
- Download ImageJ from
  https://imagej.nih.gov/ij/download.html  or
  https://mirror.imagej.net/download.html

Open ImageJ, update with `Help>Update ImageJ` to the latest version (currently, 1.53i required). Start ImageJ again and check the version.
If `Help>Update ImageJ` did not work (typically because of missing write permissions to the ImageJ directory), take the `ij.jar` file from
  http://wsr.imagej.net/download/daily-build/
and place it in the ImageJ directory (replace the old one). On MacIntosh OS X, you have to right-click the ImageJ icon and select `Show Package Contents`. The `ij.jar` file is in `Contents/Resources/Java`.

Place the `ViPErLEED_ImageJ_plugins.jar` file into ImageJ/plugins or an immediate subdirectory thereof and restart ImageJ.
The ViPErLEED plugins will appear in Plugins>ViPErLEED:
- Open Aida LEED video (opens `.vid` files of the legacy AIDA LEED data aquisition system)
- LEED Spot Tracker
- LEED I(V) Curve Editor

## Usage
The Help texts of the Spot Tracker and I(V) Curve Editor will provide some basic explanations on these plugins.
