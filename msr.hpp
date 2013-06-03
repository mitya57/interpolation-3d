#include <QtGlobal>

class MsrMatrix {
public:
	quint32 size;
	quint32 *indices;
	double *elements;
	double *rightCol;

	MsrMatrix(quint32 size, quint32 nzcount);
	~MsrMatrix();
	double getElement(quint32 i, quint32 j) const;
	quint32 solve(double *result) const;

private:
	double scalarProduct(double *vector1, double *vector2) const;
	void applyToVector(double *vector, double *newVector) const;
	double getResidual(double *vector) const;
};
