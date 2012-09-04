/*
 * OperatorSet.h
 *
 *  Created on: Apr 1, 2009
 *      Author: Michel Lesoinne
 */

#ifndef OPERATORSET_H_
#define OPERATORSET_H_

/** \brief Structure containing the mapping between the natural DOFs of the element and the constrained system DOFs
 */
struct MPCMap {
     int nCoefs;
     int *dofs; //!< Dofs to which we map
     double *coefs;
};


/** @brief Container of operators for which a mapping must be applied before matrices are added
 *
 */
template <typename Scalar>
class MappedOperatorSet {
    MPCMap **dofMaps;
  public:
    MappedOperatorSet();
    ~MappedOperatorSet();
    /** \brief Maps any matrix given with a set of dofs to a matrix with the mapped dofs before adding
     * to underlying operators
     *
     */
    void add(GenFullSquareMatrix<Scalar> &elMat, int *dofs) {
        int nDofs = elMat.dim();
        bool hasMapping = false;
        for(int i = 0; i < nDofs; ++i)
          if(dofMaps[dofs[i]]!= 0) {
            hasMapping = true;
            break;
          }
        if(hasMapping) {
          std::map<int, int> dofNumbering;
          for(int i = 0; i < nDofs; ++i) {
            std::list<int>::iterator it = dofNumbering.find()
          }
        } else {
          opSet.add(elMat, dofs);
        }
    }
};

#endif /* OPERATORSET_H_ */
