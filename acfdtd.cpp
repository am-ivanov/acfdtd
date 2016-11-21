#include "config.h"

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "rgrid/vtksaver.h"
#include "rgrid/darrayscatter.h"

//#include <fenv.h>

using namespace std;
using namespace rgrid;

inline double source(double w, double t0, double t) {
	double tmp = sqr(M_PI * w * (t - t0));
	return (1 - 2 * tmp) * exp(-tmp);
}

// first and second derivatives in pml layer
enum {
	PHI1X,
	PHI1Y,
	PHI1Z,
	PHI2X,
	PHI2Y,
	PHI2Z,
	PHINUM
};

struct PMLParams {
	double pmlVal;
	double a;
	double b;
};

int main(int argc, char** argv) {
	//feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
	rgmpi::init(&argc, &argv);
	{
	if (argc != 2) {
		cout << "Usage: " << argv[0] << " config" << endl;
		return 0;
	}
	Config cfg;
	cfg.readConfig(argv[1]);

	DArrayScatter<double, int> das;
	das.setSizes(
		Dim3D<int>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<int>(cfg.gx, cfg.gy, cfg.gz),
		Dim3D<int>(cfg.lx, cfg.ly, cfg.lz),
		Dim3D<int>(1, 1, cfg.dims==3 ? 1 : 0),
		1);
	DArrayScatter<double, int> dasNext;
	dasNext.setSizes(
		Dim3D<int>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<int>(cfg.gx, cfg.gy, cfg.gz),
		Dim3D<int>(cfg.lx, cfg.ly, cfg.lz),
		Dim3D<int>(1, 1, cfg.dims==3 ? 1 : 0),
		1);
	DArrayScatter<float, int> dasSave;
	dasSave.setSizes(
		Dim3D<int>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<int>(cfg.gx, cfg.gy, cfg.gz),
		Dim3D<int>(1, 1, 1),
		Dim3D<int>(0, 0, 0),
		1);

	VTKSaver<float, int> vs;
	vs.setHeaderStructPoints(
		"Acoustic FDTD with CPML",
		Dim3D<int>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<float>(0, 0, 0),
		Dim3D<float>(cfg.dx, cfg.dy, cfg.dz));
	vs.appendData("u", dasSave, VTKSaver<float, int>::POINT_DATA);

	DArrayScatter<double, int> x1das, y1das, z1das;
	x1das.setSizes(
		cfg.rhox.getWidth(),
		cfg.rhox.getLocalContainer().getWidth(),
		Dim3D<int>(1, 0, 0),
		1);
	y1das.setSizes(
		cfg.rhoy.getWidth(),
		cfg.rhoy.getLocalContainer().getWidth(),
		Dim3D<int>(0, 1, 0),
		1);
	if (cfg.dims == 3) {
		z1das.setSizes(
			cfg.rhoz.getWidth(),
			cfg.rhoz.getLocalContainer().getWidth(),
			Dim3D<int>(0, 0, 1),
			1);
	} else {
		z1das.setSizes(
			Dim3D<int>(cfg.gx * cfg.lx, cfg.gy*cfg.ly, cfg.gz*cfg.lz),
			Dim3D<int>(cfg.gx, cfg.gy, cfg.gz),
			Dim3D<int>(cfg.lx, cfg.ly, cfg.lz),
			Dim3D<int>(0,0,0),
			1);
	}

	DArrayScatter<double, int>* u = &das;
	DArrayScatter<double, int>* un = &dasNext;

	DArray<double, int> pml[3][2];

	pml[0][0].resize(cfg.pml_len, cfg.ny, cfg.nz, cfg.pml_len, cfg.ny, cfg.nz, 0, 0, 0, 0, 0, 0, PHINUM);
	pml[0][1].resize(cfg.pml_len, cfg.ny, cfg.nz, cfg.pml_len, cfg.ny, cfg.nz, 0, 0, 0, 0, 0, 0, PHINUM);
	pml[1][0].resize(cfg.nx, cfg.pml_len, cfg.nz, cfg.nx, cfg.pml_len, cfg.nz, 0, 0, 0, 0, 0, 0, PHINUM);
	pml[1][1].resize(cfg.nx, cfg.pml_len, cfg.nz, cfg.nx, cfg.pml_len, cfg.nz, 0, 0, 0, 0, 0, 0, PHINUM);
	if (cfg.dims == 3) {
		pml[2][0].resize(cfg.nx, cfg.ny, cfg.pml_len, cfg.nx, cfg.ny, cfg.pml_len, 0, 0, 0, 0, 0, 0, PHINUM);
		pml[2][1].resize(cfg.nx, cfg.ny, cfg.pml_len, cfg.nx, cfg.ny, cfg.pml_len, 0, 0, 0, 0, 0, 0, PHINUM);
	}

	vector<PMLParams> pmlParams1; // first derivative
	vector<PMLParams> pmlParams2; // second derivative
	for (int i = 0; i != cfg.pml_len; ++i) {
		PMLParams p;
		p.pmlVal = cfg.max_pml * sqr(1 - (i + 0.5)/cfg.pml_len);
		p.b = exp(-p.pmlVal*cfg.dt);
		p.a = p.b - 1;
		pmlParams1.push_back(p);
	}
	for (int i = 0; i != cfg.pml_len; ++i) {
		PMLParams p;
		p.pmlVal = cfg.max_pml * sqr(1 - (i + 1.0)/cfg.pml_len);
		p.b = exp(-p.pmlVal * cfg.dt);
		p.a = p.b - 1;
		pmlParams2.push_back(p);
	}

	for (int i = 0; i != 3; ++i)
		for (int j = 0; j != 2; ++j)
			pml[i][j].fill(0);

	//fstream recvsout;
	//if (rgmpi::worldRank() == 0) {
	//	recvsout.open("rcvs_out", ios_base::out | ios_base::trunc);
	//	if (!recvsout.is_open()) throw logic_error("Can't open receivers file");
	//}

	double t = 0.0;
	for (int step = 0; step != cfg.steps; ++step) {
		if (step % cfg.saveStep == 0) {
			// save to VTK
			if (rgmpi::worldRank() == 0)
				cout << "Step " << step << endl;
			char name[100];
			sprintf(name, "out-%06d.vtk", step);
			DArray<float, int>& ds = dasSave.getDArrayPart(0, 0, 0);
			DArrayContainer<double, int>& dac = das.getLocalContainer();
			for (int k = 0; k != dac.numParts(Z); ++k)
			for (int j = 0; j != dac.numParts(Y); ++j)
			for (int i = 0; i != dac.numParts(X); ++i) {
				DArray<double, int>& d = dac.getDArrayPart(i, j, k);
				for (int k2 = 0; k2 != d.localSize(Z); ++k2)
				for (int j2 = 0; j2 != d.localSize(Y); ++j2)
				for (int i2 = 0; i2 != d.localSize(X); ++i2) {
					ds(i2 + dac.partOrigin(X, i), j2 + dac.partOrigin(Y, j), k2 + dac.partOrigin(Z, k), 0) = d(i2, j2, k2, 0);
				}
			}
			ds.inverseBytes();
			vs.save(std::string(name));
		}

		DArrayContainer<double, int>& dac = u->getLocalContainer();
		DArrayContainer<double, int>& dacNext = un->getLocalContainer();

		u->externalSyncStart();
		u->internalSync();
		u->externalSyncEnd();

		for (int gk = 0; gk != dac.numParts(Z); ++gk)
			for (int gj = 0; gj != dac.numParts(Y); ++gj)
				for (int gi = 0; gi != dac.numParts(X); ++gi) {
					DArray<double, int>& p = dac.getDArrayPart(gi, gj, gk);
					DArray<double, int>& x1 = x1das.getDArrayPart(gi, gj, gk);
					DArray<double, int>& y1 = y1das.getDArrayPart(gi, gj, gk);
					DArray<double, int>& z1 = z1das.getDArrayPart(gi, gj, gk);

					// save receivers
					//recvsout << t;
					//for (vector<Receiver>::size_type i = 0; i != cfg.rcv.size(); ++i) {
					//	recvsout << " " << u->val(cfg.rcv.at(i).x, cfg.rcv.at(i).y, cfg.rcv.at(i).z, 0);
					//}
					//recvsout << endl;

					// inner area
					for (int k = 0; k != p.localSize(Z); ++k) {
						for (int j = 0; j != p.localSize(Y); ++j) {
							for (int i = 0; i != p.localSize(X)+1; ++i) {
								int i2 = i + p.origin(X);
								int j2 = j + p.origin(Y);
								int k2 = k + p.origin(Z);
								x1(i, j, k, 0) = (p(i, j, k, 0) - p(i-1, j, k, 0)) / cfg.dx;
								if (i2 < cfg.pml_len && cfg.isPml[0][0]) {
									pml[0][0](i2,j2,k2,PHI1X) = pmlParams1.at(i2).b * pml[0][0](i2,j2,k2,PHI1X) + pmlParams1.at(i).a * x1(i,j,k,0);
									x1(i,j,k,0) = x1(i,j,k,0) + pml[0][0](i2,j2,k2,PHI1X);
								} else if (i2 > cfg.nx - cfg.pml_len && cfg.isPml[0][1]) {
									pml[0][1](cfg.nx-i2,j2,k2,PHI1X) = pmlParams1.at(cfg.nx-i2).b * pml[0][1](cfg.nx-i2,j2,k2,PHI1X) + pmlParams1.at(cfg.nx-i2).a * x1(i,j,k,0);
									x1(i,j,k,0) = x1(i,j,k,0) + pml[0][1](cfg.nx-i2,j2,k2,PHI1X);
								}
							}
						}
					}

					for (int k = 0; k != p.localSize(Z); ++k) {
						for (int j = 0; j != p.localSize(Y)+1; ++j) {
							for (int i = 0; i != p.localSize(X); ++i) {
								int i2 = i + p.origin(X);
								int j2 = j + p.origin(Y);
								int k2 = k + p.origin(Z);
								y1(i,j,k,0) = (p(i,j,k,0) - p(i,j-1,k,0)) / cfg.dy;
								if (j2 < cfg.pml_len && cfg.isPml[1][0]) {
									pml[1][0](i2,j2,k2,PHI1Y) = pmlParams1.at(j2).b * pml[1][0](i2,j2,k2,PHI1Y) + pmlParams1.at(j2).a * y1(i,j,k,0);
									y1(i,j,k,0) = y1(i,j,k,0) + pml[1][0](i2,j2,k2,PHI1Y);
								} else if (j2 > cfg.ny - cfg.pml_len && cfg.isPml[1][1]) {
									pml[1][1](i2,cfg.ny-j2,k2,PHI1Y) = pmlParams1.at(cfg.ny-j2).b * pml[1][1](i2,cfg.ny-j2,k2,PHI1Y) + pmlParams1.at(cfg.ny-j2).a * y1(i,j,k,0);
									y1(i,j,k,0) = y1(i,j,k,0) + pml[1][1](i2,cfg.ny-j2,k2,PHI1Y);
								}
							}
						}
					}

					if (cfg.dims == 3)
						for (int k = 0; k != p.localSize(Z)+1; ++k) {
							for (int j = 0; j != p.localSize(Y); ++j) {
								for (int i = 0; i != p.localSize(X); ++i) {
									int i2 = i + p.origin(X);
									int j2 = j + p.origin(Y);
									int k2 = k + p.origin(Z);
									z1(i,j,k,0) = (p(i,j,k,0) - p(i,j,k-1,0)) / cfg.dz;
									if (k2 < cfg.pml_len && cfg.isPml[2][0]) {
										pml[2][0](i2,j2,k2,PHI1Z) = pmlParams1.at(k2).b * pml[2][0](i2,j2,k2,PHI1Z) + pmlParams1.at(k2).a * z1(i,j,k,0);
										z1(i,j,k,0) = z1(i,j,k,0) + pml[2][0](i2,j2,k2,PHI1Z);
									} else if (k2 > cfg.nz - cfg.pml_len && cfg.isPml[2][1]) {
										pml[2][1](i2,j2,cfg.nz-k2,PHI1Z) = pmlParams1.at(cfg.nz-k2).b * pml[2][1](i2,j2,cfg.nz-k2,PHI1Z) + pmlParams1.at(cfg.nz-k2).a * z1(i,j,k,0);
										z1(i,j,k,0) = z1(i,j,k,0) + pml[2][1](i2,j2,cfg.nz-k2,PHI1Z);
									}
								}
							}
						}
		//			}

		/*
		 * Next synchronizations are not necessary, because
		 * we calculate up to nx+1, ny+1, nz+1 indexes in previous block
		 */
		//x1das.externalSyncStart();
		//y1das.externalSyncStart();
		//if (cfg.dims == 3) z1das.externalSyncStart();

		//x1das.internalSync();
		//y1das.internalSync();
		//if (cfg.dims == 3) z1das.internalSync();

		//x1das.externalSyncEnd();
		//y1das.externalSyncEnd();
		//if (cfg.dims == 3) z1das.externalSyncEnd();

		//for (int gk = 0; gk != dac.numParts(Z); ++gk)
		//	for (int gj = 0; gj != dac.numParts(Y); ++gj)
		//		for (int gi = 0; gi != dac.numParts(X); ++gi) {
		//			DArray<double, int>& p = dac.getDArrayPart(gi, gj, gk);
					DArray<double, int>& pn = dacNext.getDArrayPart(gi, gj, gk);
					DArray<double, int>& K = cfg.K.getDArrayPart(gi, gj, gk);
					DArray<double, int>& rhox = cfg.rhox.getDArrayPart(gi, gj, gk);
					DArray<double, int>& rhoy = cfg.rhoy.getDArrayPart(gi, gj, gk);
					DArray<double, int>& rhoz = cfg.rhoz.getDArrayPart(gi, gj, gk);
		//			DArray<double, int>& x1 = x1das.getDArrayPart(gi, gj, gk);
		//			DArray<double, int>& y1 = y1das.getDArrayPart(gi, gj, gk);
		//			DArray<double, int>& z1 = z1das.getDArrayPart(gi, gj, gk);

					for (int k = 0; k != p.localSize(Z); ++k) {
						for (int j = 0; j != p.localSize(Y); ++j) {
							for (int i = 0; i != p.localSize(X); ++i) {
								double x2 = 0, y2 = 0, z2 = 0;
								int i2 = i + p.origin(X);
								int j2 = j + p.origin(Y);
								int k2 = k + p.origin(Z);

								x2 = (x1(i+1,j,k,0)/rhox(i+1,j,k,0) - x1(i,j,k,0)/rhox(i,j,k,0)) / cfg.dx;
								if (i2 < cfg.pml_len && cfg.isPml[0][0]) {
									pml[0][0](i2,j2,k2,PHI2X) = pmlParams2.at(i2).b * pml[0][0](i2,j2,k2,PHI2X) + pmlParams2.at(i2).a * x2;
									x2 = x2 + pml[0][0](i2,j2,k2,PHI2X);
								} else if (i2 > cfg.nx - cfg.pml_len - 1 && cfg.isPml[0][1]) {
									pml[0][1](cfg.nx-i2-1,j2,k2,PHI2X) = pmlParams2.at(cfg.nx-i2-1).b * pml[0][1](cfg.nx-i2-1,j2,k2,PHI2X) + pmlParams2.at(cfg.nx-i2-1).a * x2;
									x2 = x2 + pml[0][1](cfg.nx-i2-1,j2,k2,PHI2X);
								}

								y2 = (y1(i,j+1,k,0)/rhoy(i,j+1,k,0) - y1(i,j,k,0)/rhoy(i,j,k,0)) / cfg.dy;
								if (j2 < cfg.pml_len && cfg.isPml[1][0]) {
									pml[1][0](i2,j2,k2,PHI2Y) = pmlParams2.at(j2).b * pml[1][0](i2,j2,k2,PHI2Y) + pmlParams2.at(j).a * y2;
									y2 = y2 + pml[1][0](i2,j2,k2,PHI2Y);
								} else if (j2 > cfg.ny - cfg.pml_len - 1 && cfg.isPml[1][1]) {
									pml[1][1](i2,cfg.ny-j2-1,k2,PHI2Y) = pmlParams2.at(cfg.ny-j2-1).b * pml[1][1](i2,cfg.ny-j2-1,k2,PHI2Y) + pmlParams2.at(cfg.ny-j2-1).a * y2;
									y2 = y2 + pml[1][1](i2,cfg.ny-j2-1,k2,PHI2Y);
								}

								if (cfg.dims == 3) {
									z2 = (z1(i,j,k+1,0)/rhoz(i,j,k+1,0) - z1(i,j,k,0)/rhoz(i,j,k,0)) / cfg.dz;
									if (k2 < cfg.pml_len && cfg.isPml[2][0]) {
										pml[2][0](i2,j2,k2,PHI2Z) = pmlParams2.at(k2).b * pml[2][0](i2,j2,k2,PHI2Z) + pmlParams2.at(k2).a * z2;
										z2 = z2 + pml[2][0](i2,j2,k2,PHI2Z);
									} else if (k2 > cfg.nz - cfg.pml_len - 1 && cfg.isPml[2][1]) {
										pml[2][1](i2,j2,cfg.nz-k2-1,PHI2Z) = pmlParams2.at(cfg.nz-k2-1).b * pml[2][1](i2,j2,cfg.nz-k2-1,PHI2Z) + pmlParams2.at(cfg.nz-k2-1).a * z2;
										z2 = z2 + pml[2][1](i2,j2,cfg.nz-k2-1,PHI2Z);
									}
								}

								pn(i, j, k, 0) =
									2.0 * p.val(i, j, k, 0) - pn(i, j, k, 0)
									+ K(i, j, k, 0) * sqr(cfg.dt) * (x2 + y2 + z2);
							}
						}
					}

					// insert source
					for (vector<Source>::iterator it = cfg.src.begin(); it != cfg.src.end(); ++it)
					{
						if (pn.check(it->x, it->y, it->z)) {
							int i = it->x-pn.origin(X);
							int j = it->y-pn.origin(Y);
							int k = it->z-pn.origin(Z);
							pn(i, j, k, 0) +=
								sqr(cfg.dt) * K(i,j,k,0) * it->val.at(step);
						}
					}
				}

		swap(un, u);

		t += cfg.dt;
	}

	//if (rgmpi::worldRank() == 0) {
	//	recvsout.close();
	//}
	}
	rgmpi::forceFinalize();
}
