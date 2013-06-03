#include "threads.hpp"
#include <QDebug>

void synchronize(Args *a) {
	pthread_mutex_lock(a->mutex);
	++(*(a->donethreads));
	if (*(a->donethreads) == a->threads)
		pthread_cond_signal(a->condvar_tomf);
	pthread_cond_wait(a->condvar_totf, a->mutex);
	pthread_mutex_unlock(a->mutex);
}

void taskGetNzCount(Args *a) {
	quint32 i, j;
	TriangleList triangleList;
	for (i = a->thread; i < POINTS(a->layers, a->size); i += a->threads) {
		a->locnzc[i] = 0;
		for (j = 0; j < POINTS(a->layers, a->size); ++j) {
			if (i == j)
				continue;
			triangleList = getCommonTriangles(a->points, a->size, a->layers,
				getVertex(a->size, i), getVertex(a->size, j), a->centerPtr);
			if (!(triangleList.isEmpty()))
				++(a->locnzc[i]);
		}
	}
}

void taskFillMatrix(Args *a) {
	quint32 i, j, l, ind, msize = POINTS(a->layers, a->size);
	int k;
	double c;
	TriangleList triangleList;
	Polynom polynom, interpolationPolynom;
	QPointF vertex1, vertex2;
	QPointF intPoints[3];
	ind = msize + 1;
	for (i = 0; i < msize; ++i) {
		if (i % a->threads != a->thread) {
			ind += a->locnzc[i];
			continue;
		}
		vertex1 = getPoint(a->points, a->size, a->layers,
			getVertex(a->size, i), a->centerPtr);
		a->matrix->indices[i] = ind;
		for (j = 0; j < msize; ++j) {
			if (i == j) {
				triangleList = getSurroundingTriangles(a->points, a->size,
					a->layers, getVertex(a->size, i), a->centerPtr);
				// fill right part while we have that list
				a->matrix->rightCol[i] = 0;
				for (k = 0; k < triangleList.size(); ++k) {
					Triangle t = triangleList.at(k);
					for (l = 0; l < 4; ++l) {
						intPoints[0] = (l < 3) ? t[l] : (t[1] + t[2]) / 2;
						intPoints[1] = (t[l % 3] + t[(l + 1) % 3]) / 2;
						intPoints[2] = (t[l % 3] + t[(l + 2) % 3]) / 2;
						interpolationPolynom = getLinearInterpolation(intPoints);
						polynom = polynomMultiply(
							phiFunctionPolynom(t, vertex1),
							interpolationPolynom);
						a->matrix->rightCol[i] += getIntegral(intPoints, polynom);
					}
				}
			} else
				triangleList = getCommonTriangles(a->points, a->size, a->layers,
				getVertex(a->size, i), getVertex(a->size, j), a->centerPtr);
			if (triangleList.isEmpty())
				continue;
			c = 0;
			vertex2 = getPoint(a->points, a->size, a->layers,
				getVertex(a->size, j), a->centerPtr);
			for (k = 0; k < triangleList.size(); ++k) {
				polynom = polynomMultiply(
					phiFunctionPolynom(triangleList.at(k), vertex1),
					phiFunctionPolynom(triangleList.at(k), vertex2));
				c += getIntegral(triangleList.at(k), polynom);
			}
			if (i == j)
				a->matrix->elements[i] = c;
			else {
				a->matrix->elements[ind] = c;
				a->matrix->indices[ind] = j;
				++ind;
			}
		}
	}
	a->matrix->indices[msize] = ind;
}

void taskFillValues(Args *a) {
	quint32 i, j;
	double phi;
	QPointF point;
	for (j = a->thread; j < POINTS(a->drawLayers, a->size); j += a->threads) {
		a->values[j] = 0;
		point = getPoint(a->points, a->size, a->drawLayers,
			getVertex(a->size, j), a->centerPtr);
		for (i = 0; i < POINTS(a->layers, a->size); ++i) {
			phi = phiFunction(a->points, a->size, a->layers, point, i,
				a->centerPtr);
			a->values[j] += a->alphas[i] * phi;
		}
	}
}

void *tf(void *args) {
	Args *a = (Args *)args;
	synchronize(a);
	while (a->task) {
		if (a->task == GetNzCount)
			taskGetNzCount(a);
		if (a->task == FillAndSolveMatrix)
			taskFillMatrix(a);
		if (a->task == FillValues)
			taskFillValues(a);
		synchronize(a);
	}
	return NULL;
}
