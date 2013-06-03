#include <pthread.h>
#include <QtGlobal>
#include "3d.hpp"

enum Task {
	FinishWork         = 0,
	GetNzCount         = 1,
	FillAndSolveMatrix = 2,
	FillValues         = 4
};

struct Args {
	Task task;
	pthread_mutex_t *mutex;
	pthread_cond_t *condvar_tomf;
	pthread_cond_t *condvar_totf;
	quint32 thread;
	quint32 threads;
	quint32 *donethreads;
	MsrMatrix *matrix;
	QPointF *points;
	QPointF *centerPtr;
	quint32 size;
	quint32 layers;
	quint32 drawLayers;
	quint32 *locnzc;
	double *alphas;
	double *values;
};

void *tf(void *args);
