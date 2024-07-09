/** This class implements an expandable float array similar
 *  to an ArrayList or Vector but more efficient because it uses primitive types.
 *  Calls to this class are not synchronized.
 */

/** This code is part of the ViPErLEED package for LEED I(V) analysis.
 *  Licensed under GNU General Public License v3.0 or later (GPL-3.0-or-later),
 *  https://www.gnu.org/licenses/gpl-3.0.html
 *  The authors may decide later to put part of the auxiliary code in this work into the public domain,
 *  to allow incorporation into ImageJ if desired (ImageJ is in the public domain).
 *  When using and/or modifying this program for scientific work, please cite
 *  the paper describing it:
 *  M. Schmid, F. Kraushofer, A. M. Imre, T. Kißlinger, L. Hammer, U. Diebold, and M. Riva,
 *  ViPErLEED package II: Spot tracking, extraction and processing of I(V) curves,
 *  Phys. Rev. Research, 2024. 
 *  @author Michael Schmid, IAP/TU Wien, 2019-2024
 */

 public class LeedIntegerArray {
    private int size;
    private int[] data;

    /** Creates a new expandable array with an initial capacity of 100. */
    public LeedIntegerArray() {
        this(100);
    }

    /** Creates a new expandable array with specified initial capacity.
     * @throws IllegalArgumentException if the specified initial capacity is less than zero */
    public LeedIntegerArray(int initialCapacity) {
        if (initialCapacity < 0) throw new IllegalArgumentException("Illegal IntArray Capacity: "+initialCapacity);
        data = new int[initialCapacity];
    }

    /** Returns the number of elements in the array. */
    public int size() {
        return size;
    }

    /** Removes all elements form this IntArray. */
    public void clear() {
        size = 0;
    }

    /** Returns a new int array containing all elements of this LeedIntegerArray. */
    public int[] toArray() {
        int[] out = new int[size];
        System.arraycopy(data, 0, out, 0, size);
        return out;
    }

    /** Returns an int array containing all elements of this LeedIntegerArray.
     *  This may be the array used internally or a new array.
     *  Use this method only if you are not going to overwrite the values
     *  in the output array (or if you are sure you are not going to use
     *  this LeedIntegerArray any more). Oterwise use the toArray() method
     *  instead. */
    public int[] getArray() {
        if (size == data.length)
            return data;
        else
            return toArray();
    }

    /** Returns the element at the specified position in this LeedIntegerArray.
     *  @throws IndexOutOfBoundsException - if index is out of range (<code>index < 0 || index >= size()</code>). */
    public int get(int index) {
        if (index<0 || index>= size) throw new IndexOutOfBoundsException("IntArray Index out of Bounds: "+index);
        return data[index];
    }

    /** Returns the last element of this LeedIntegerArray.
     *  @throws IndexOutOfBoundsException - if this LeedIntegerArray is empty */
    public int getLast() {
        return get(size-1);
    }

    /** Replaces the element at the specified position with the value given.
     *  @return the value previously at the specified position.
     *  @throws IndexOutOfBoundsException - if index is out of range (<code>index < 0 || index >= size()</code>). */
    public int set(int index, int value) {
        if (index<0 || index>= size) throw new IndexOutOfBoundsException("IntArray Index out of Bounds: "+index);
        int previousValue = data[index];
        data[index] = value;
        return previousValue;
    }

	/** Increments the value at the given index. Expands the array
	 *  to the necessary size with values of 0 if the element at index does not exist yet. */
    public void increment(int index) {
		if (index >= size) {
			for (int i=size; i<=Math.min(index, data.length-1); i++)
				data[i] = 0;
			size = index + 1;
		}
		if (size > data.length) {
            int[] newData = new int[size*2 + 50];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
		}
		data[index]++;
	}

    /** Appends the specified value to the end of this LeedIntegerArray.
     *  Returns the number of elements after adding. */
    public int add(int value) {
        if (size >= data.length) {
            int[] newData = new int[size*2 + 50];
            System.arraycopy(data, 0, newData, 0, size);
            data = newData;
        }
        data[size++] = value;
        return size;
    }

    /** Appends the first n values from the specified array to this LeedIntegerArray.
     *  Returns the number of elements after adding. */
    public int add(int[] a, int n) {
        if (size + n > data.length) {
            int[] newData = new int[size*2 + n + 50];
            System.arraycopy(data, 0, newData, 0, size);
            data = newData;
        }
        System.arraycopy(a, 0, data, size, n);
        size += n;
        return size;
    }

    /** Deletes the last <code>n</code> element from this LeedIntegerArray.
     *  <code>n</code> may be larger than the number of elements; in that
     *  case, all elements are removed. */
    public void removeLast(int n) {
        size -= Math.min(n, size);
    }

    /** Deletes the value at the given index */
    public void remove( int index) {
        if (index < 0 || index >= size)
            throw new IndexOutOfBoundsException(index+" not in range 0-"+(size-1));
        System.arraycopy(data, index+1, data, index, size-index-1);
        size--;
    }

    /** Returns the index of the first occurrence of the specified value, or -1 if
     *  the value is not contained in the array */
    public int indexOf(int value) {
        for (int i=0; i<size; i++)
            if (data[i] == value) return i;
        return -1;
    }

    /** Trims the capacity of this LeedIntegerArray instance to be its current size,
     *  minimizing the storage of the LeedIntegerArray instance. */
    public void trimToSize() {
        data = toArray();
    }
}
