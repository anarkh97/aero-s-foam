#ifndef __TOREDUCEDROWECHELONFORM__
#define __TOREDUCEDROWECHELONFORM__

/**
 * Row echelon form with full pivoting. Note that columns and
 * rows might be swapped due to the full pivoting, returned in rowmap and colmap if not NULL
 *  
 * @param mx            the matr to reduce
 * @param reduced       if true, the reduced row echelon form is returned, that
 *                      is, zeros above the diagonal, too
 * @param rowmap        the row mapping to reestablish the original row ordering. 
 * @param colmap        the column mapping to reestablish the original column ordering.
 * @return              the rank of the matr
 */

#include <limits>

template <typename T, typename Matrix>
int rowEchelon(Matrix& mx, bool reduced, int* rowmap, int* colmap) 
{
  int rows = mx.rows();
  int cols = mx.cols();
  int pivs = std::min<int>(rows, cols);
        
  //find pivot row/column
  int prow = -1;
  int pcol = -1;
  T pval = -std::numeric_limits<T>::max();
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
      T val = std::abs<T>(mx(row, col));
      if (val > pval) {
        pval = val;
        prow = row;
        pcol = col;
      }
    }
  }
        
  //precondition (each iteration): prow/pcol/pval are set
  for (int pivot = 0; pivot < pivs; pivot++) {
    if (pval <= 10*std::numeric_limits<T>::epsilon()) return pivot;
          
    //swap rows / columns
    if (prow != pivot) {
      mx.row(prow).swap(mx.row(pivot));
      if (rowmap != NULL) {
        int tmp = rowmap[prow]; 
        rowmap[prow] = rowmap[pivot]; 
        rowmap[pivot] = tmp;
      }
    }

    if (pcol != pivot) {
      mx.col(pcol).swap(mx.col(pivot));
      if (colmap != NULL) {
        int tmp = colmap[pcol];
        colmap[pcol] = colmap[pivot];
        colmap[pivot] = tmp;
      }
    }
    
    //divide pivot row
    pval = mx(pivot, pivot);
    for (int col = pivot + 1; col < cols; col++) {
      mx(pivot,col) /= pval;
    }
    mx(pivot,pivot) = 1.0;

    //subtract pivot row from other rows
    //find next pivot at the same time
    pval = -std::numeric_limits<T>::max();
    for (int row = pivot + 1; row < rows; row++) {
      T rpiv = mx(row, pivot);                          
      mx(row, pivot) = 0.0;
      for (int col = pivot + 1; col < cols; col++) {
        T val = mx(row, col);
        T sub = mx(pivot, col);
        val -= sub * rpiv;
        mx(row, col) = val;
        // is this our new pivot?
        if (row < rows && col < cols) {
          if (val < 0.0) val = -val;
          if (val > pval) {
            pval = val;
            prow = row;
            pcol = col;
          }
        }
      }
    }
    if (reduced) {
      //subtract pivot from rows above pivot row, too
      for (int row = 0; row < pivot; row++) {
        T rpiv = mx(row, pivot);                          
        mx(row, pivot) = 0.0;
        for (int col = pivot + 1; col < cols; col++) {
          T val = mx(row, col);
          T sub = mx(pivot, col);
          mx(row, col) = val - sub * rpiv;
        }
      }
    }
  }
  return pivs;
}

template <typename T, typename Matrix>
void ToReducedRowEchelonForm(Matrix& m, int *rowmap) 
{
  int lead = 0;
  int rowCount = m.rows();
  int colCount = m.cols();

  int i;
  T lv;
 
  for(int r = 0; r < rowCount; r++) {
    if(lead >= colCount)
      return;
    i = r;
    //while(m(i,lead) == 0) {
    while(abs<T>(m(i,lead)) <= std::numeric_limits<T>::epsilon()) {
      i++;
      if(i == rowCount) {
        i = r;
        lead++;
        if(lead == colCount)
          return;
      }
    }

    m.row(i).swap(m.row(r));
    if (rowmap != NULL) {
      int tmp = rowmap[i];
      rowmap[i] = rowmap[r];
      rowmap[r] = tmp;
    }

    lv = m(r,lead);
    m.row(r) /= lv;
    for(i = 0; i < rowCount; i++) {
      if(i != r) {
        lv = m(i, lead);
        m.row(i) -= m.row(r)*lv;
      }
    }
    lead++;
  }
}

#endif
