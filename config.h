#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>

#include "common.h"
#include "rgrid/darrayscatter.h"

class Config {
public:
	void readConfig(std::string file);
	void readRhoX(std::string file);
	void readRhoY(std::string file);
	void readRhoZ(std::string file);
	void readK(std::string file);
	void readSources(std::string file);
	void readReceivers(std::string file);

	long dims;
	long nx, ny, nz, steps;
	double dx, dy, dz, dt;
	long lx, ly, lz;
	long gx, gy, gz;
	long saveStep;

	bool isPml[3][2];
	long pml_len;
	double max_pml;

	std::string rcvsOut;

// Ready to work data
	std::vector<Source> src;
	std::vector<Receiver> rcv;

	rgrid::DArrayScatter<double, int> K;
	rgrid::DArrayScatter<double, int> rhox, rhoy, rhoz;
};

#endif // CONFIG_H
