#include "msr.hpp"

MsrMatrix::MsrMatrix(quint32 size, quint32 nzcount):
	size(size),
	indices(new quint32[size+nzcount+1]),
	elements(new double[size+nzcount+1]),
	rightCol(new double[size])
{}

MsrMatrix::~MsrMatrix() {
	delete[] elements;
	delete[] indices;
	delete[] rightCol;
}

double MsrMatrix::getElement(quint32 i, quint32 j) const {
	quint32 index;
	if (i == j)
		return elements[i];
	else
		for (index = indices[i]; index < indices[i+1]; ++index)
			if (indices[index] == j)
				return elements[index];
	return 0;
}

void MsrMatrix::applyToVector(double *vector, double *newVector) const {
	quint32 i, ind;
	for (i = 0; i < size; ++i) {
		newVector[i] = elements[i] * vector[i];
		for (ind = indices[i]; ind < indices[i+1]; ++ind)
			newVector[i] += elements[ind] * vector[indices[ind]];
	}
}
double MsrMatrix::getResidual(double *vector) const {
	double result = 0, c;
	quint32 i, ind;
	for (i = 0; i < size; ++i) {
		c = elements[i] * vector[i];
		for (ind = indices[i]; ind < indices[i+1]; ++ind)
			c += elements[ind] * vector[indices[ind]];
		result += (rightCol[i] - c) * (rightCol[i] - c);
	}
	return result;
}

double MsrMatrix::scalarProduct(double *vector1, double *vector2) const {
	double result = 0;
	for (quint32 i = 0; i < size; ++i)
		result += vector1[i] * vector2[i];
	return result;
}

quint32 MsrMatrix::solve(double *result) const {
	double *r = new double[size];
	double *p = new double[size];
	double alpha, residual, prevResidual = 0;
	quint32 i;
	for (i = 0; i < size; ++i) {
		result[i] = 0;
		r[i] = rightCol[i];
	}
	applyToVector(r, p);
	residual = getResidual(result);
	quint32 iterations = 0;
	while (qAbs(prevResidual - residual) > 1e-14) {
		++iterations;
		alpha = scalarProduct(p, r) / scalarProduct(p, p);
		//qDebug() << "Alpha is" << alpha;
		for (i = 0; i < size; ++i) {
			result[i] += alpha * r[i];
			r[i] -= alpha * p[i];
		}
		applyToVector(r, p);
		prevResidual = residual;
		residual = getResidual(result);
		//qDebug() << "Residual is" << residual;
	}
	delete[] r;
	delete[] p;
	return iterations;
}
