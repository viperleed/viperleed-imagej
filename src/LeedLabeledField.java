import ij.IJ;
import ij.gui.GenericDialog;
import java.util.Vector;
import java.awt.Component;
import java.awt.TextField;
import java.awt.Color;
import java.io.File;
import java.io.FileReader;


/**
 *  This class is used for enabling and disabling more than one
 *  dialog components with a single call.
 *  Currently only supports a numeric field of a GenericDialog,
 *  where the associated label is enabled/disabled together with
 *  input field:
 *     genericDialog.addNumericField(...);
 *     LeedLabeledField leedInputSet = LeedLabeledField.numeric(genericDialog);
 *     ...
 *     boolean inputAllowed = ...
 *     leedInputSet.setEnabled(inputAllowed);
 */

/** This ImageJ plugin is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  This means that the collection of ViPErLEED ImageJ plugins is free software:
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation,
 *  either version 3 of the License, or (at your option) any later version.
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  The collection of ViPErLEED ImageJ plugins is distributed
 *  in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  The authors may decide later to put part of the auxiliary code in this work into
 *  the public domain, to allow incorporation into ImageJ if desired.
 *  (ImageJ is in the public domain.)
 *
 *  The help text is licensed under the Creative Commons Attribution 4.0 (CC BY 4.0) license.
 *
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research 7, 013006 (2025).
 *  @author Michael Schmid, IAP/TU Wien, 2024-2025
 */


public class LeedLabeledField {
    /** The array of Component that should be enabled/disabled together */
    private Component[] components;

    /** Stores the components for later enabline/disabling.
     *  @param components An array of java.awt.Components; must have at least one element
     *  and must not contain null. */
    private LeedLabeledField(Component[] components) {
        this.components = components;
    }

    /** Creates a LeedLabeledField for a numeric field.
     *  This must be called immediately after gd.addNumericField
     *  @param gd The GenericDialog whose addNumericField method was called. */
    public static LeedLabeledField numeric(GenericDialog gd) {
        Component label = gd.getLabel();
        Component field = (Component)gd.getNumericFields().lastElement();
        Component[] components = new Component[]{field, label};
        return new LeedLabeledField(components);
    }

    /** Enables or disables (grays out) the Components.
     *  @param b If true, the components are enabled; otherwise the components are disabled (grayed out) */
    public void setEnabled(boolean b) {
        for (Component c : components) {
            c.setEnabled(b);
        }
    }

    /** Returns whether the first of the components is enabled */
    public boolean isEnabled() {
        return components[0].isEnabled();
    }

    /** Returns the text if the components contain a TextField (note that a GenericDialog
     *  numeric field is a text field), or null */
    public String getText() {
        for (Component c : components)
            if (c instanceof TextField)
                return ((TextField)c).getText();
        return null;
    }
}
