#include <qmath.h>
#include "3d.hpp"

void transformateTriangle(QPointF *triangle, double *transform) {
	QPointF a = triangle[0], b = triangle[1], c = triangle[2];
	transform[0] = 1 / (b.x() - a.x());
	transform[2] = -transform[0] * a.x();
	transform[1] = 1 / (c.y() - a.y());
	transform[3] = -transform[1] * a.y();
}

void addToPolynom(Polynom &p, quint32 xPow, quint32 yPow, double c) {
	if (qFuzzyIsNull(c))
		return;
	for (int i = 0; i < p.size(); ++i)
		if (p.at(i).xPow == xPow && p.at(i).yPow == yPow) {
			p[i].c += c;
			return;
		}
	Monom m;
	m.xPow = xPow;
	m.yPow = yPow;
	m.c = c;
	p.append(m);
}

void addToPolynom(Polynom &p, const Monom &m) {
	for (int i = 0; i < p.size(); ++i)
		if (p.at(i).xPow == m.xPow && p.at(i).yPow == m.yPow) {
			p[i].c += m.c;
			return;
		}
	if (!qFuzzyIsNull(m.c))
		p.append(m);
}

Polynom polynomAdd(const Polynom &p1, const Polynom &p2) {
	Polynom result = p1;
	for (int i = 0; i < p2.size(); ++i)
		addToPolynom(result, p2.at(i));
	return result;
}

Polynom polynomMultiply(const Polynom &p1, const Polynom &p2) {
	Polynom result;
	int i, j;
	Monom m;
	for (i = 0; i < p1.size(); ++i)
		for (j = 0; j < p2.size(); ++j) {
			m.xPow = p1.at(i).xPow + p2.at(j).xPow;
			m.yPow = p1.at(i).yPow + p2.at(j).yPow;
			m.c = p1.at(i).c * p2.at(j).c;
			addToPolynom(result, m);
		}
	return result;
}

double polynomIntegral(const Polynom &p) {
	// all xPow's should be zero
	double result = 0;
	for (int i = 0; i < p.size(); ++i)
		result += p.at(i).c / (1 + p.at(i).yPow);
	return result;
}

Polynom polynomPower(const Polynom &p, quint32 power) {
	Polynom result;
	if (!power) {
		addToPolynom(result, 0, 0, 1);
		return result;
	}
	result = p;
	for (quint32 i = 0; i < power - 1; ++i)
		result = polynomMultiply(result, p);
	return result;
}

double getIntegralForGoodTriangle(QPointF *triangle, const Polynom &polynom) {
	double transform[4];
	transformateTriangle(triangle, transform);
	double newcx = transform[0] * triangle[2].x() + transform[2];
	Polynom newPolynom, polynomPart1, polynomPart2, polynomBase;
	int i, j;
	for (i = 0; i < polynom.size(); ++i) {
		polynomPart1.clear();
		polynomBase.clear();
		addToPolynom(polynomBase, 1, 0, 1 / transform[0]);
		addToPolynom(polynomBase, 0, 0, -transform[2] / transform[0]);
		polynomPart1 = polynomPower(polynomBase, polynom.at(i).xPow);
		polynomPart2.clear();
		polynomBase.clear();
		addToPolynom(polynomBase, 0, 1, 1 / transform[1]);
		addToPolynom(polynomBase, 0, 0, -transform[3] / transform[1]);
		polynomPart2 = polynomPower(polynomBase, polynom.at(i).yPow);
		for (j = 0; j < polynomPart2.size(); ++j)
			polynomPart2[j].c *= polynom.at(i).c;
		newPolynom = polynomAdd(newPolynom,
			polynomMultiply(polynomPart1, polynomPart2));
	}
	double result = 0;
	for (i = 0; i < newPolynom.size(); ++i) {
		polynomBase.clear();
		addToPolynom(polynomBase, 0, 0, 1);
		addToPolynom(polynomBase, 0, 1, newcx - 1);
		polynomBase = polynomPower(polynomBase, newPolynom.at(i).xPow + 1);
		addToPolynom(polynomBase, 0,
			newPolynom.at(i).xPow + 1, qPow(newcx, newPolynom.at(i).xPow + 1));
		for (j = 0; j < polynomBase.size(); ++j) {
			polynomBase[j].c /= newPolynom.at(i).xPow + 1;
			polynomBase[j].yPow += newPolynom.at(i).yPow;
		}
		result += newPolynom.at(i).c * polynomIntegral(polynomBase) /
			qAbs(transform[0] * transform[1]);
	}
	return result;
}

double getIntegral(QPointF *points, const Polynom &polynom) {
	if (points[0].y() > points[1].y())
		qSwap(points[0], points[1]);
	if (points[0].y() > points[2].y())
		qSwap(points[0], points[2]);
	if (points[1].y() > points[2].y())
		qSwap(points[1], points[2]);
	if (qFuzzyIsNull(points[0].y() - points[1].y()))
		return getIntegralForGoodTriangle(points, polynom);
	if (qFuzzyIsNull(points[1].y() - points[2].y())) {
		qSwap(points[0], points[2]);
		return getIntegralForGoodTriangle(points, polynom);
	}
	double result = 0;
	double h1 = points[1].y() - points[0].y(),
	       h2 = points[2].y() - points[1].y();
	QPointF triangle[3];
	triangle[0] = points[1];
	triangle[1] = (points[0] * h2 + points[2] * h1) / (h1 + h2);
	triangle[2] = points[0];
	result = getIntegralForGoodTriangle(triangle, polynom);
	triangle[2] = points[2];
	return result + getIntegralForGoodTriangle(triangle, polynom);
}

double getIntegral(const Triangle &triangle, const Polynom &polynom) {
	QPointF points[3];
	points[0] = triangle[0];
	points[1] = triangle[1];
	points[2] = triangle[2];
	return getIntegral(points, polynom);
}
