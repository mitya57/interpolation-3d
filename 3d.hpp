#include <QPointF>
#include <QList>
#include "msr.hpp"

#define POINT(l, s, i) (1+size*l*(l-1)/2+s*l+i)
#define TRIANGLE(l, s, i) getTriangle(points, size, layers, l, s, i, centerPtr)

#define POINTS(ls, sz) (sz*(ls+1)*ls/2+1)
#define TRIANGLES(ls, sz) (sz*ls*ls)

quint32 getPoints(QPointF **points);

QPointF getPoint(QPointF *points, quint32 size, quint32 layer, quint32 sector,
quint32 layers, quint32 index, QPointF *centerPtr);

struct Vertex {
	quint32 layer;
	quint32 sector;
	quint32 index;
};

class Triangle {
private:
	QPointF *points;
	QPointF *centerPtr;
	quint32 size;
	quint32 layers;

public:
	Vertex a;
	Vertex b;
	Vertex c;
	Triangle(QPointF *points, quint32 size, quint32 layers, QPointF *centerPtr,
	Vertex a, Vertex b, Vertex c):
		points(points), centerPtr(centerPtr), size(size), layers(layers),
		a(a), b(b), c(c)
	{}
	QPointF operator[](quint32 i) const {
		Vertex v;
		switch(i) {
			case 0: v = a; break;
			case 1: v = b; break;
			case 2: v = c; break;
			default: return QPointF(); // to fix a warning
		}
		return getPoint(points, size, v.layer, v.sector, layers, v.index, centerPtr);
	}
};

struct Monom {
	quint32 xPow;
	quint32 yPow;
	double c;
};

typedef QList<Triangle> TriangleList;
typedef QList<Monom> Polynom;

QPointF getPoint(QPointF *points, quint32 size, quint32 layers,
Vertex vertex, QPointF *centerPtr);

QPointF getCenter(QPointF *points, quint32 size);

Triangle getTriangle(QPointF *points, quint32 size, quint32 layers, quint32 index,
QPointF *centerPtr);

Triangle getTriangle(QPointF *points, quint32 size, quint32 layers,
quint32 layer, quint32 sector, quint32 index, QPointF *centerPtr);

Vertex getVertex(quint32 size, quint32 number);

quint32 getVertexIndex(Vertex v, quint32 size);

double phiFunction(QPointF *points, quint32 size, quint32 layers,
QPointF point, quint32 vertexNumber, QPointF *centerPtr);

TriangleList getCommonTriangles(QPointF *points, quint32 size, quint32 layers,
Vertex v1, Vertex v2, QPointF *centerPtr);

TriangleList getSurroundingTriangles(QPointF *points, quint32 size, quint32 layers,
Vertex v, QPointF *centerPtr);

Polynom getLinearInterpolation(QPointF *p);

Polynom phiFunctionPolynom(const Triangle &triangle, QPointF vertex);

void addToPolynom(Polynom &p, quint32 xPow, quint32 yPow, double c);

Polynom polynomMultiply(const Polynom &p1, const Polynom &p2);

double getIntegral(QPointF *triangle, const Polynom &polynom);

double getIntegral(const Triangle &triangle, const Polynom &polynom);
