#include <QApplication>
#include <cfenv>
#include "window.hpp"

int main(int argc, char **argv) {
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	QApplication app(argc, argv);
	MyMainWindow window;
	window.show();
	return app.exec();
}
