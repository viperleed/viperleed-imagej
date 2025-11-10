import ij.IJ;
import ij.util.Tools;
import java.util.concurrent.atomic.AtomicIntegerArray;


/**
 *  A float array with atomic set/get and also a method to ensure a minimum value for
 *  an array element in a multi-threaded environment.
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
 *  @author Michael Schmid, IAP/TU Wien, 2020-2025
 */



public class LeedAtomicFloatArray {
    AtomicIntegerArray iArray;

    /** Creates a new LeedAtomicFloatArray of given length with given initial values for all elements */
    public LeedAtomicFloatArray(int length, float initialValue) {
        iArray = new AtomicIntegerArray(length);
	for (int i=0; i<length; i++)
	    iArray.set(i, Float.floatToIntBits(initialValue));
    }

    /** Returns the element at position i */
    public float get(int i) {
	return Float.intBitsToFloat(iArray.get(i));
    }

    /** Sets the element at position i to the given value */
    public void set(int i, float value) {
	iArray.set(i, Float.floatToIntBits(value));
    }

    /** Sets the element at position i to the given value unless the value of that element is already smaller */
    public void setMinimumValue(int i, float value) {
	int iCurrent;
	do {
	    iCurrent = iArray.get(i);
	    if (value >= Float.intBitsToFloat(iCurrent)) return;
	} while (!iArray.compareAndSet(i, /*expected=*/iCurrent, /*new=*/Float.floatToIntBits(value)));
    }

}
