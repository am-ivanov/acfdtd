#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char** argv) {
	if (argc != 2) {
		cout << "Usage: " << argv[0] << " config.txt" << endl;
		cout << "Config format:" << endl;
		cout << "filename val x y z" << endl;
		cout << "val lx ly lz rx ry rz" << endl;
		exit(0);
	}

	fstream fs(argv[1], ios_base::in);
	if (!fs.is_open()) throw logic_error("Can't open file");

	double val;
	string filename;
	long nx, ny, nz;
	fs >> filename >> val >> nx >> ny >> nz;

	vector<double> data(nx * ny * nz, val);

	while (true) {
		long lx, ly, lz;
		long rx, ry, rz;
		fs >> val >> lx >> ly >> lz >> rx >> ry >> rz;
		if (fs.eof()) break;
		for (int k = lz; k < rz; ++k)
		for (int j = ly; j < ry; ++j)
		for (int i = lx; i < rx; ++i) {
			data[i + j * nx + k * nx * ny] = val;
		}
	}

	fs.close();

	fs.open(filename.c_str(), ios_base::out | ios_base::trunc);
	for (int k = 0; k < nz; ++k)
	for (int j = 0; j < ny; ++j)
	for (int i = 0; i < nx; ++i) {
		fs.write((char*)&data[i + j * nx + k * nx * ny], sizeof(double));
	}
	fs.close();
}
