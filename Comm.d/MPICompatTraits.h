//
// Created by Michel Lesoinne on 2/13/18.
//

#ifndef FEM_MPICOMPATTRAITS_H
#define FEM_MPICOMPATTRAITS_H

template <>
struct CommTypeCompatibility<MPI_Win, HandleType::Window> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Comm, HandleType::Communicator> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Datatype , HandleType::Type> {
	static const bool isCompatible = true;
};

template <>
struct CommTypeCompatibility<MPI_Op , HandleType::Op> {
	static const bool isCompatible = true;
};
#endif //FEM_MPICOMPATTRAITS_H
