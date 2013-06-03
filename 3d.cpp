#include <QFile>
#include <QLineF>
#include <QTextStream>
#include <qmath.h>

#include "function.hpp"
#include "3d.hpp"

#define NORMALIZE_ANGLE(phi) (phi > 180) ? (360 - phi) : phi
#define EPS .02

quint32 getPoints(QPointF **points) {
	quint32 size;
	QFile coordinatesFile("coordinates.txt");
	coordinatesFile.open(QFile::ReadOnly);
	QTextStream stream(&coordinatesFile);
	stream >> size;
	*points = new QPointF[size];
	qreal x, y;
	for (quint32 i = 0; i < size; ++i) {
		stream >> x >> y;
		(*points)[i] = QPointF(x, y);
	}
	coordinatesFile.close();
	return size;
}

QPointF getCenter(QPointF *points, quint32 size) {
	QPointF center(0, 0);
	for (quint32 i = 0; i < size; ++i)
		center += points[i] / size;
	return center;
}

QPointF getPoint(QPointF *points, quint32 size, quint32 layer, quint32 sector,
quint32 layers, quint32 index, QPointF *centerPtr) {
	QPointF center = centerPtr ? *centerPtr : getCenter(points, size);
	if (!layer)
		return center;
	QPointF result = (points[sector] * (layer - index) + points[(sector+1)%size] * index) / layer;
	return (result * layer + center * (layers - layer)) / layers;
}

QPointF getPoint(QPointF *points, quint32 size, quint32 layers,
Vertex vertex, QPointF *centerPtr) {
	return getPoint(points, size, vertex.layer, vertex.sector, layers,
		vertex.index, centerPtr);
}

quint32 getVertexIndex(Vertex v, quint32 size) {
	return POINT(v.layer, v.sector, v.index);
}

Vertex getVertex(quint32 layer, quint32 sector, quint32 index, quint32 size) {
	Vertex v;
	if (layer && index >= layer) {
		sector = (sector + 1) % size;
		index -= layer;
	}
	v.layer = layer;
	v.sector = sector;
	v.index = index;
	return v;
}

Vertex getVertex(quint32 size, quint32 number) {
	quint32 layer = 1, sector, index;
	if (!number)
		layer = sector = index = 0;
	else {
		layer = 1;
		quint32 layerStart = 1;
		while (number >= layerStart + layer * size) {
			layerStart += layer * size;
			++layer;
		}
		sector = (number - layerStart) / layer;
		index  = (number - layerStart) % layer;
	}
	return getVertex(layer, sector, index, 0);
}

Triangle getTriangle(QPointF *points, quint32 size, quint32 layers,
quint32 layer, quint32 sector, quint32 index, QPointF *centerPtr) {
	Vertex a, b, c;
	if (index & 1) {
		a = getVertex(layer,     sector, index / 2,     size);
		b = getVertex(layer + 1, sector, index / 2 + 1, size);
		c = getVertex(layer,     sector, index / 2 + 1, size);
	} else {
		a = getVertex(layer,     sector, index / 2,     size);
		b = getVertex(layer + 1, sector, index / 2,     size);
		c = getVertex(layer + 1, sector, index / 2 + 1, size);
	}
	return Triangle(points, size, layers, centerPtr, a, b, c);
}

Triangle getTriangle(QPointF *points, quint32 size, quint32 layers, quint32 number,
QPointF *centerPtr) {
	quint32 layer = 0;
	quint32 layerStart = 0;
	while (number >= layerStart + (layer * 2 + 1) * size) {
		layerStart += (layer * 2 + 1) * size;
		++layer;
	}
	quint32 trianglesInSector = layer * 2 + 1;
	quint32 sector = (number - layerStart) / trianglesInSector;
	quint32 index  = (number - layerStart) % trianglesInSector;
	return getTriangle(points, size, layers, layer, sector, index, centerPtr);
}

Polynom psiFunctionPolynom(const Triangle &triangle, short n) {
	double u1, u2, v1, v2;
	QPointF a = triangle[0], b = triangle[1], c = triangle[2];
	Polynom result;
	if (n == 0) {
		u1 = b.x();
		v1 = b.y();
	} else {
		u1 = a.x();
		v1 = a.y();
	}
	if (n == 2) {
		u2 = b.x();
		v2 = b.y();
	} else {
		u2 = c.x();
		v2 = c.y();
	}
	addToPolynom(result, 1, 0, v2 - v1);
	addToPolynom(result, 0, 1, u1 - u2);
	addToPolynom(result, 0, 0, u2 * v1 - v2 * u1);
	return result;
}

double psiFunction(QPointF point, const Triangle &triangle, short n) {
	double u1, u2, v1, v2;
	QPointF a = triangle[0], b = triangle[1], c = triangle[2];
	if (n == 0) {
		u1 = b.x();
		v1 = b.y();
	} else {
		u1 = a.x();
		v1 = a.y();
	}
	if (n == 2) {
		u2 = b.x();
		v2 = b.y();
	} else {
		u2 = c.x();
		v2 = c.y();
	}
	return (point.x() - u1) * (v2 - v1) - (point.y() - v1) * (u2 - u1);
}

bool pointInTriangle(QPointF point, const Triangle &triangle) {
	QPointF a = triangle[0], b = triangle[1], c = triangle[2];
	if (point == a || point == b || point == c)
		return true;
	QLineF line0(point, a), line1(point, b), line2(point, c);
	qreal angle0 = line0.angleTo(line1);
	angle0 = NORMALIZE_ANGLE(angle0);
	qreal angle1 = line1.angleTo(line2);
	angle1 = NORMALIZE_ANGLE(angle1);
	qreal angle2 = line2.angleTo(line0);
	angle2 = NORMALIZE_ANGLE(angle2);
	return qAbs(angle0 + angle1 + angle2 - 360) < EPS;
}

double phiFunction(QPointF point, const Triangle &triangle, QPointF vertex) {
	int n = 0;
	while (triangle[n] != vertex)
		++n;
	Q_ASSERT(qAbs(psiFunction(vertex, triangle, n)) > 1e-10);
	Q_ASSERT(qAbs(psiFunction(point, triangle, n)) < 1e10);
	return psiFunction(point, triangle, n) / psiFunction(vertex, triangle, n);
}

Polynom phiFunctionPolynom(const Triangle &triangle, QPointF vertex) {
	int n = 0;
	while (triangle[n] != vertex)
		++n;
	double psiValue = psiFunction(vertex, triangle, n);
	Polynom result = psiFunctionPolynom(triangle, n);
	for (n = 0; n < result.size(); ++n)
		result[n].c /= psiValue;
	return result;
}

TriangleList getCommonTriangles(QPointF *points, quint32 size, quint32 layers,
Vertex v1, Vertex v2, QPointF *centerPtr) {
	TriangleList result;
	if (v1.layer > v2.layer)
		qSwap(v1, v2);
	if (!v2.layer)
		return result;
	if (v1.layer == 0 && v2.layer == 1) {
		result.append(TRIANGLE(0, v2.sector, 0));
		result.append(TRIANGLE(0, (v2.sector + size - 1) % size, 0));
	} else if (v1.sector == v2.sector) {
		if (v1.layer == v2.layer) {
			if (v2.index == v1.index + 1) {
				if (v1.layer < layers)
					result.append(TRIANGLE(v1.layer, v1.sector, v1.index * 2 + 1));
				result.append(TRIANGLE(v1.layer - 1, v1.sector, v1.index * 2));
			} else if (v2.index == v1.index - 1) {
				if (v1.layer < layers)
					result.append(TRIANGLE(v1.layer, v1.sector, v2.index * 2 + 1));
				result.append(TRIANGLE(v1.layer - 1, v1.sector, v2.index * 2));
			}
		}
		else if (v2.layer == v1.layer + 1) {
			if (v2.index == v1.index) {
				if (v1.index)
					result.append(TRIANGLE(v1.layer, v1.sector, v1.index * 2 - 1));
				else
					result.append(TRIANGLE(v1.layer, (v1.sector + size - 1) % size, v1.layer * 2));
				result.append(TRIANGLE(v1.layer, v1.sector, v1.index * 2));
			} else if (v2.index == v1.index + 1) {
				result.append(TRIANGLE(v1.layer, v1.sector, v1.index * 2));
				result.append(TRIANGLE(v1.layer, v1.sector, v1.index * 2 + 1));
			}
		}
	} else if ((v2.sector == (v1.sector + 1) % size)
	&& v1.index == v1.layer - 1 && v2.index == 0) {
		if (v1.layer == v2.layer) {
			result.append(TRIANGLE(v1.layer - 1, v1.sector, v1.index * 2));
			if (v1.layer < layers)
				result.append(TRIANGLE(v1.layer, v1.sector, v1.index * 2 + 1));
		} else if (v1.layer == v2.layer + 1) {
			result.append(TRIANGLE(v2.layer, v1.sector, v1.index * 2 - 1));
			result.append(TRIANGLE(v2.layer, v1.sector, v1.index * 2));
		}
	} else if ((v1.sector == (v2.sector + 1) % size)
	&& v2.index == v2.layer - 1 && v1.index == 0) {
		if (v1.layer == v2.layer) {
			result.append(TRIANGLE(v1.layer - 1, v2.sector, v2.index * 2));
			if (v1.layer < layers)
				result.append(TRIANGLE(v1.layer, v2.sector, v2.index * 2 + 1));
		} else if (v1.layer == v2.layer - 1) {
			result.append(TRIANGLE(v1.layer, v2.sector, v2.index * 2 - 1));
			result.append(TRIANGLE(v1.layer, v2.sector, v2.index * 2));
		}
	}
	return result;
}

TriangleList getSurroundingTriangles(QPointF *points, quint32 size, quint32 layers,
Vertex v, QPointF *centerPtr) {
	TriangleList result;
	if (!v.layer) {
		for (quint32 i = 0; i < size; ++i)
			result.append(TRIANGLE(0, i, 0));
		return result;
	}
	// inner part
	result.append(TRIANGLE(v.layer - 1, v.sector, 2 * v.index));
	if (v.index) {
		result.append(TRIANGLE(v.layer - 1, v.sector, 2 * v.index - 1));
		result.append(TRIANGLE(v.layer - 1, v.sector, 2 * v.index - 2));
	} else
		result.append(TRIANGLE(v.layer - 1, (v.sector + size - 1) % size, 2 * v.layer - 2));
	// outer part
	if (v.layer < layers) {
		result.append(TRIANGLE(v.layer, v.sector, 2 * v.index));
		result.append(TRIANGLE(v.layer, v.sector, 2 * v.index + 1));
		if (v.index) {
			result.append(TRIANGLE(v.layer, v.sector, 2 * v.index - 1));
		} else {
			result.append(TRIANGLE(v.layer, (v.sector + size - 1) % size, 2 * v.layer - 1));
			result.append(TRIANGLE(v.layer, (v.sector + size - 1) % size, 2 * v.layer));
		}
	}
	return result;
}

double phiFunction(QPointF *points, quint32 size, quint32 layers,
QPointF point, quint32 vertexNumber, QPointF *centerPtr) {
	Vertex vertex = getVertex(size, vertexNumber);
	TriangleList triangles = getSurroundingTriangles(points, size,
		layers, vertex, centerPtr);
	for (int i = 0; i < triangles.size(); ++i) {
		if (pointInTriangle(point, triangles.at(i)))
			return phiFunction(point, triangles.at(i),
			getPoint(points, size, layers, vertex, centerPtr));
	}
	return 0;
}

double determinant(double a11, double a12, double a13, double a21,
double a22, double a23, double a31, double a32, double a33) {
	return a11 * a22 * a33 + a21 * a32 * a13 + a12 * a23 * a31
	     - a31 * a22 * a13 - a11 * a23 * a32 - a21 * a12 * a33;
}

Polynom getLinearInterpolation(QPointF *p) {
	Polynom result;
	double a, b, c, d;
	d = determinant(
		p[0].x(), p[0].y(), 1,
		p[1].x(), p[1].y(), 1,
		p[2].x(), p[2].y(), 1
	);
	a = determinant(
		function(p[0].x(), p[0].y()), p[0].y(), 1,
		function(p[1].x(), p[1].y()), p[1].y(), 1,
		function(p[2].x(), p[2].y()), p[2].y(), 1
	) / d;
	b = determinant(
		p[0].x(), function(p[0].x(), p[0].y()), 1,
		p[1].x(), function(p[1].x(), p[1].y()), 1,
		p[2].x(), function(p[2].x(), p[2].y()), 1
	) / d;
	c = determinant(
		p[0].x(), p[0].y(), function(p[0].x(), p[0].y()),
		p[1].x(), p[1].y(), function(p[1].x(), p[1].y()),
		p[2].x(), p[2].y(), function(p[2].x(), p[2].y())
	) / d;
	addToPolynom(result, 1, 0, a);
	addToPolynom(result, 0, 1, b);
	addToPolynom(result, 0, 0, c);
	return result;
}
