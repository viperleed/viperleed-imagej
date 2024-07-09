# viperleed-imagej
ViPErLEED-ImageJ is a collection of ImageJ plugins for spot tracking and intensity evaluation as well as selection and editing of I(V) curves created by the Spot Tracker.
The output of the I(V) Curve Editor can be used as experimental data (`EXPBEAMS.csv`) for structure optimization.

## Installation
- Download ImageJ from [https://imagej.net/ij/download.html](https://imagej.net/ij/download.html)

- Open ImageJ, update with `Help>Update ImageJ` to the latest version (currently, 1.54k required; select "daily build" if it is not on the list yet). Start ImageJ again and check the version.
If `Help>Update ImageJ` did not work (typically because of missing write permissions to the ImageJ directory), take the `ij.jar` file from [http://wsr.imagej.net/download/daily-build/](http://wsr.imagej.net/download/daily-build/) and place it in the ImageJ directory (replace the old one). On MacIntosh OS X, you have to right-click the ImageJ icon and select `Show Package Contents`. The `ij.jar` file is in `Contents/Resources/Java`.

- Place the `ViPErLEED_ImageJ_plugins.jar` file (from this repository) into ImageJ/plugins or a subdirectory thereof and restart ImageJ. The ViPErLEED plugins will appear in Plugins>ViPErLEED.

- If you want to have the ViPErLEED plugins available with a single mouse click, add `run("LEED Commands");` to the "AutoRun" macro of your ImageJ StartupMacros file (in the ImageJ/macros directory) or to the list of startup commands by selecting Edit>Options>Startup from the ImageJ Menu.

## Usage
The Help texts of the plugins will provide some basic explanations on these plugins.
For more information, see the documentation available at [viperleed.org](https://www.viperleed.org/content/imagej_plugins.html) and the paper at [arXiv:2406.18413](https://arxiv.org/abs/2406.18413).
