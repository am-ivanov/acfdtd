#include "config.h"

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "rgrid/vtksaver.h"
#include "rgrid/darrayscatter.h"

#include <fenv.h>

using namespace std;
using namespace rgrid;

inline real_t source(real_t w, real_t t0, real_t t) {
	real_t tmp = sqr(M_PI * w * (t - t0));
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
	real_t pmlVal;
	real_t a;
	real_t b;
};



int main(int argc, char** argv) {
	feenableexcept(FE_INVALID | FE_OVERFLOW);  // Enable all floating point exceptions but FE_INEXACT
	rgmpi::init(&argc, &argv);
	{
	if (argc != 2) {
		cout << "Usage: " << argv[0] << " config" << endl;
		return 0;
	}
	Config cfg;
	cfg.readConfig(argv[1]);

	for (int i = -cfg.ho; i <= cfg.ho; ++i) {
		if (i != 0)
			cout << cfg.fdc.sc1(i) << " ";
	}
	cout << endl;

	DArrayScatter<real_t, int_t> das;
	das.setSizes(
		Dim3D<int_t>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<int_t>(cfg.gx, cfg.gy, cfg.gz),
		Dim3D<int_t>(cfg.lx, cfg.ly, cfg.lz),
		Dim3D<int_t>(cfg.ho, cfg.ho, cfg.dims==3 ? cfg.ho : 0),
		1);
	DArrayScatter<real_t, int_t> dasNext;
	dasNext.setSizes(
		Dim3D<int_t>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<int_t>(cfg.gx, cfg.gy, cfg.gz),
		Dim3D<int_t>(cfg.lx, cfg.ly, cfg.lz),
		Dim3D<int_t>(cfg.ho, cfg.ho, cfg.dims==3 ? cfg.ho : 0),
		1);
	DArrayScatter<float, int_t> dasSave;
	dasSave.setSizes(
		Dim3D<int_t>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<int_t>(cfg.gx, cfg.gy, cfg.gz),
		Dim3D<int_t>(1, 1, 1),
		Dim3D<int_t>(0, 0, 0),
		1);

	VTKSaver<float, int_t> vs;
	vs.setHeaderStructPoints(
		"Acoustic FDTD with CPML",
		Dim3D<int_t>(cfg.nx, cfg.ny, cfg.nz),
		Dim3D<float>(0, 0, 0),
		Dim3D<float>(cfg.dx, cfg.dy, cfg.dz));
	vs.appendData("u", dasSave, VTKSaver<float, int_t>::POINT_DATA);

	DArrayScatter<real_t, int_t> x1das, y1das, z1das;
	x1das.setSizes(
		cfg.rhox.getWidth(),
		cfg.rhox.getLocalContainer().getWidth(),
		Dim3D<int_t>(cfg.ho, 0, 0),
		1);
	y1das.setSizes(
		cfg.rhoy.getWidth(),
		cfg.rhoy.getLocalContainer().getWidth(),
		Dim3D<int_t>(0, cfg.ho, 0),
		1);
	if (cfg.dims == 3) {
		z1das.setSizes(
			cfg.rhoz.getWidth(),
			cfg.rhoz.getLocalContainer().getWidth(),
			Dim3D<int_t>(0, 0, cfg.ho),
			1);
	} else {
		z1das.setSizes(
			Dim3D<int_t>(cfg.gx * cfg.lx, cfg.gy*cfg.ly, cfg.gz*cfg.lz),
			Dim3D<int_t>(cfg.gx, cfg.gy, cfg.gz),
			Dim3D<int_t>(cfg.lx, cfg.ly, cfg.lz),
			Dim3D<int_t>(0, 0, 0),
			1);
	}

	DArrayScatter<real_t, int_t>* u = &das;
	DArrayScatter<real_t, int_t>* un = &dasNext;

	DArray<real_t, int_t> pml[3][2];

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
	for (int_t i = 0; i != cfg.pml_len; ++i) {
		PMLParams p;
		p.pmlVal = cfg.max_pml * sqr(1 - (i + 0.5)/cfg.pml_len);
		p.b = exp(-p.pmlVal*cfg.dt);
		p.a = p.b - 1;
		pmlParams1.push_back(p);
	}
	for (int_t i = 0; i != cfg.pml_len; ++i) {
		PMLParams p;
		p.pmlVal = cfg.max_pml * sqr(1 - (i * 1.0)/cfg.pml_len);
		p.b = exp(-p.pmlVal * cfg.dt);
		p.a = p.b - 1;
		pmlParams2.push_back(p);
	}

	for (int_t i = 0; i != 3; ++i)
		for (int_t j = 0; j != 2; ++j)
			pml[i][j].fill(0);

#ifdef USE_MPI
	MPI_File fh;
	if (0 != MPI_File_open(MPI_COMM_WORLD, cfg.rcvsOut.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh))
		throw logic_error("Can't open receivers file for writing");
	if (0 != MPI_File_set_size(fh, 0))
		throw logic_error("Can't truncate receivers file");
#else
	fstream recvsout;
	if (rgmpi::worldRank() == 0) {
		recvsout.open(cfg.rcvsOut.c_str(), ios_base::out | ios_base::trunc);
		if (!recvsout.is_open()) throw logic_error("Can't open receivers file for writing");
	}
#endif

	real_t t = 0.0;
	for (int_t step = 0; step != cfg.steps; ++step) {

		DArrayContainer<real_t, int_t>& dac = u->getLocalContainer();
		DArrayContainer<real_t, int_t>& dacNext = un->getLocalContainer();

		if (step % cfg.saveStep == 0) {
			// save to VTK
			if (rgmpi::worldRank() == 0)
				cout << "Step " << step << endl;
			char name[100];
			sprintf(name, "out-%06ld.vtk", (long)step);
			DArray<float, int_t>& ds = dasSave.getDArrayPart(0, 0, 0);
			for (int_t k = 0; k != dac.numParts(Z); ++k)
			for (int_t j = 0; j != dac.numParts(Y); ++j)
			for (int_t i = 0; i != dac.numParts(X); ++i) {
				DArray<real_t, int_t>& d = dac.getDArrayPart(i, j, k);
				for (int_t k2 = 0; k2 != d.localSize(Z); ++k2)
				for (int_t j2 = 0; j2 != d.localSize(Y); ++j2)
				for (int_t i2 = 0; i2 != d.localSize(X); ++i2) {
					ds(i2 + dac.partOrigin(X, i), j2 + dac.partOrigin(Y, j), k2 + dac.partOrigin(Z, k), 0) = d(i2, j2, k2, 0);
				}
			}
			ds.inverseBytes();
			vs.save(std::string(name));
		}

		u->externalSyncStart();
		u->internalSync();
		u->externalSyncEnd();

		// save receivers
		for (vector<Receiver>::size_type i = 0; i != cfg.rcv.size(); ++i) {
			const Dim3D<int_t> rind(cfg.rcv.at(i).i, cfg.rcv.at(i).j, cfg.rcv.at(i).k);
			if (u->isPresentGhostGlobal(rind)) {
				float val[LR][LR][LR] = {{{0}}};
				val[L][L][L] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X],rind[Y],rind[Z]),0) * cfg.rcv.at(i).c[L][L][L];
				val[R][L][L] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X]+1,rind[Y],rind[Z]),0) * cfg.rcv.at(i).c[R][L][L];
				val[L][R][L] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X],rind[Y]+1,rind[Z]),0) * cfg.rcv.at(i).c[L][R][L];
				val[R][R][L] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X]+1,rind[Y]+1,rind[Z]),0) * cfg.rcv.at(i).c[R][R][L];
				if (cfg.dims == 3) {
					val[L][L][R] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X],rind[Y],rind[Z]+1),0) * cfg.rcv.at(i).c[L][L][R];
					val[R][L][R] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X]+1,rind[Y],rind[Z]+1),0) * cfg.rcv.at(i).c[R][L][R];
					val[L][R][R] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X],rind[Y]+1,rind[Z]+1),0) * cfg.rcv.at(i).c[L][R][R];
					val[R][R][R] = u->getNodeGhostLocal(Dim3D<int_t>(rind[X]+1,rind[Y]+1,rind[Z]+1),0) * cfg.rcv.at(i).c[R][R][R];
				}
				float rval =
					val[L][L][L] +
					val[R][L][L] +
					val[L][R][L] +
					val[R][R][L] +
					val[L][L][R] +
					val[R][L][R] +
					val[L][R][R] +
					val[R][R][R];
#ifdef USE_MPI
				if (0 != MPI_File_write_at(fh, (step * cfg.rcv.size() + i) * sizeof(float), &rval, 1, MPI_FLOAT, MPI_STATUS_IGNORE))
					throw logic_error("Can't write in receivers file");
#else
				recvsout.write(reinterpret_cast<const char*>(&rval), sizeof(float));
#endif
			}
		}

		// recalculate own DArrays
		for (int_t gk = 0; gk != dac.numParts(Z); ++gk)
			for (int_t gj = 0; gj != dac.numParts(Y); ++gj)
				for (int_t gi = 0; gi != dac.numParts(X); ++gi) {
					DArray<real_t, int_t>& p = dac.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& x1 = x1das.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& y1 = y1das.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& z1 = z1das.getDArrayPart(gi, gj, gk);

					DArray<real_t, int_t>& rhox = cfg.rhox.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& rhoy = cfg.rhoy.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& rhoz = cfg.rhoz.getDArrayPart(gi, gj, gk);

					// inner area
					for (int_t k = 0; k != rhox.localSize(Z); ++k) {
						for (int_t j = 0; j != rhox.localSize(Y); ++j) {
							for (int_t i = 0; i != rhox.localSize(X); ++i) {
								int_t i2 = i + rhox.origin(X);
								int_t j2 = j + rhox.origin(Y);
								int_t k2 = k + rhox.origin(Z);
								x1(i, j, k, 0) = 0;
								for (int_t n = 0; n < cfg.ho; ++n)
									x1(i, j, k, 0) += cfg.fdc.sc1(n+1) * (p(i+1+n, j, k, 0) - p(i-n, j, k, 0)) / cfg.dx;
								if (i2 < cfg.pml_len && cfg.isPml[0][0]) {
									pml[0][0](i2,j2,k2,PHI1X) = pmlParams1.at(i2).b * pml[0][0](i2,j2,k2,PHI1X) + pmlParams1.at(i).a * x1(i,j,k,0);
									x1(i,j,k,0) = x1(i,j,k,0) + pml[0][0](i2,j2,k2,PHI1X);
								} else if (i2 > cfg.nx-cfg.pml_len-2 && cfg.isPml[0][1]) {
									pml[0][1](cfg.nx-i2-2,j2,k2,PHI1X) = pmlParams1.at(cfg.nx-i2-2).b * pml[0][1](cfg.nx-i2-2,j2,k2,PHI1X) + pmlParams1.at(cfg.nx-i2-2).a * x1(i,j,k,0);
									x1(i,j,k,0) = x1(i,j,k,0) + pml[0][1](cfg.nx-i2-2,j2,k2,PHI1X);
								}
							}
						}
					}

					for (int_t k = 0; k != rhoy.localSize(Z); ++k) {
						for (int_t j = 0; j != rhoy.localSize(Y); ++j) {
							for (int_t i = 0; i != rhoy.localSize(X); ++i) {
								int_t i2 = i + rhoy.origin(X);
								int_t j2 = j + rhoy.origin(Y);
								int_t k2 = k + rhoy.origin(Z);
								y1(i,j,k,0) = 0;
								for (int_t n = 0; n < cfg.ho; ++n)
									y1(i,j,k,0) += cfg.fdc.sc1(n+1) * (p(i,j+1+n,k,0) - p(i,j-n,k,0)) / cfg.dy;
								if (j2 < cfg.pml_len && cfg.isPml[1][0]) {
									pml[1][0](i2,j2,k2,PHI1Y) = pmlParams1.at(j2).b * pml[1][0](i2,j2,k2,PHI1Y) + pmlParams1.at(j2).a * y1(i,j,k,0);
									y1(i,j,k,0) = y1(i,j,k,0) + pml[1][0](i2,j2,k2,PHI1Y);
								} else if (j2 > cfg.ny-cfg.pml_len-2 && cfg.isPml[1][1]) {
									pml[1][1](i2,cfg.ny-j2-2,k2,PHI1Y) = pmlParams1.at(cfg.ny-j2-2).b * pml[1][1](i2,cfg.ny-j2-2,k2,PHI1Y) + pmlParams1.at(cfg.ny-j2-2).a * y1(i,j,k,0);
									y1(i,j,k,0) = y1(i,j,k,0) + pml[1][1](i2,cfg.ny-j2-2,k2,PHI1Y);
								}
							}
						}
					}

					if (cfg.dims == 3)
						for (int_t k = 0; k != rhoz.localSize(Z); ++k) {
							for (int_t j = 0; j != rhoz.localSize(Y); ++j) {
								for (int_t i = 0; i != rhoz.localSize(X); ++i) {
									int_t i2 = i + rhoz.origin(X);
									int_t j2 = j + rhoz.origin(Y);
									int_t k2 = k + rhoz.origin(Z);
									z1(i,j,k,0);
									for (int_t n = 0; n < cfg.ho; ++n)
										z1(i,j,k,0) += cfg.fdc.sc1(n+1) * (p(i,j,k+1+n,0) - p(i,j,k-n,0)) / cfg.dz;
									if (k2 < cfg.pml_len && cfg.isPml[2][0]) {
										pml[2][0](i2,j2,k2,PHI1Z) = pmlParams1.at(k2).b * pml[2][0](i2,j2,k2,PHI1Z) + pmlParams1.at(k2).a * z1(i,j,k,0);
										z1(i,j,k,0) = z1(i,j,k,0) + pml[2][0](i2,j2,k2,PHI1Z);
									} else if (k2 > cfg.nz-cfg.pml_len-2 && cfg.isPml[2][1]) {
										pml[2][1](i2,j2,cfg.nz-k2-2,PHI1Z) = pmlParams1.at(cfg.nz-k2-2).b * pml[2][1](i2,j2,cfg.nz-k2-2,PHI1Z) + pmlParams1.at(cfg.nz-k2-2).a * z1(i,j,k,0);
										z1(i,j,k,0) = z1(i,j,k,0) + pml[2][1](i2,j2,cfg.nz-k2-2,PHI1Z);
									}
								}
							}
						}
				}

		x1das.externalSyncStart();
		y1das.externalSyncStart();
		z1das.externalSyncStart();

		x1das.internalSync();
		y1das.internalSync();
		z1das.internalSync();

		x1das.externalSyncEnd();
		y1das.externalSyncEnd();
		z1das.externalSyncEnd();

		for (int_t gk = 0; gk != dac.numParts(Z); ++gk)
			for (int_t gj = 0; gj != dac.numParts(Y); ++gj)
				for (int_t gi = 0; gi != dac.numParts(X); ++gi) {

					DArray<real_t, int_t>& p = dac.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& pn = dacNext.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& K = cfg.K.getDArrayPart(gi, gj, gk);

					DArray<real_t, int_t>& x1 = x1das.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& y1 = y1das.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& z1 = z1das.getDArrayPart(gi, gj, gk);

					DArray<real_t, int_t>& rhox = cfg.rhox.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& rhoy = cfg.rhoy.getDArrayPart(gi, gj, gk);
					DArray<real_t, int_t>& rhoz = cfg.rhoz.getDArrayPart(gi, gj, gk);

					for (int_t k = 0; k != p.localSize(Z); ++k) {
						for (int_t j = 0; j != p.localSize(Y); ++j) {
							for (int_t i = 0; i != p.localSize(X); ++i) {
								real_t x2 = 0, y2 = 0, z2 = 0;
								int_t i2 = i + p.origin(X);
								int_t j2 = j + p.origin(Y);
								int_t k2 = k + p.origin(Z);

								for (int_t n = 0; n < cfg.ho; ++n)
									x2 += cfg.fdc.sc1(n+1) * (x1(i+n,j,k,0)/rhox(i+n,j,k,0) - x1(i-1-n,j,k,0)/rhox(i-1-n,j,k,0)) / cfg.dx;
								if (i2 < cfg.pml_len && cfg.isPml[0][0]) {
									pml[0][0](i2,j2,k2,PHI2X) = pmlParams2.at(i2).b * pml[0][0](i2,j2,k2,PHI2X) + pmlParams2.at(i2).a * x2;
									x2 = x2 + pml[0][0](i2,j2,k2,PHI2X);
								} else if (i2 > cfg.nx - cfg.pml_len - 1 && cfg.isPml[0][1]) {
									pml[0][1](cfg.nx-i2-1,j2,k2,PHI2X) = pmlParams2.at(cfg.nx-i2-1).b * pml[0][1](cfg.nx-i2-1,j2,k2,PHI2X) + pmlParams2.at(cfg.nx-i2-1).a * x2;
									x2 = x2 + pml[0][1](cfg.nx-i2-1,j2,k2,PHI2X);
								}

								for (int_t n = 0; n < cfg.ho; ++n)
									y2 += cfg.fdc.sc1(n+1) * (y1(i,j+n,k,0)/rhoy(i,j+n,k,0) - y1(i,j-1-n,k,0)/rhoy(i,j-1-n,k,0)) / cfg.dy;
								if (j2 < cfg.pml_len && cfg.isPml[1][0]) {
									pml[1][0](i2,j2,k2,PHI2Y) = pmlParams2.at(j2).b * pml[1][0](i2,j2,k2,PHI2Y) + pmlParams2.at(j).a * y2;
									y2 = y2 + pml[1][0](i2,j2,k2,PHI2Y);
								} else if (j2 > cfg.ny - cfg.pml_len - 1 && cfg.isPml[1][1]) {
									pml[1][1](i2,cfg.ny-j2-1,k2,PHI2Y) = pmlParams2.at(cfg.ny-j2-1).b * pml[1][1](i2,cfg.ny-j2-1,k2,PHI2Y) + pmlParams2.at(cfg.ny-j2-1).a * y2;
									y2 = y2 + pml[1][1](i2,cfg.ny-j2-1,k2,PHI2Y);
								}

								if (cfg.dims == 3) {
									for (int_t n = 0; n < cfg.ho; ++n)
										z2 += cfg.fdc.sc1(n+1) * (z1(i,j,k+n,0)/rhoz(i,j,k+n,0) - z1(i,j,k-1-n,0)/rhoz(i,j,k-1-n,0)) / cfg.dz;
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

				}

		// insert source
		for (vector<Source>::iterator it = cfg.src.begin(); it != cfg.src.end(); ++it) {
			// put source only in real no ghost nodes
			if (un->isPresent(Dim3D<int_t>(it->i  ,it->j  ,it->k  )))
				un->getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k  ),0)
					+= it->c[L][L][L] * sqr(cfg.dt) * it->val.at(step)
					* cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k  ),0);

			if (un->isPresent(Dim3D<int_t>(it->i+1,it->j  ,it->k  )))
				un->getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k  ),0)
					+= it->c[R][L][L] * sqr(cfg.dt) * it->val.at(step)
					* cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k  ),0);

			if (un->isPresent(Dim3D<int_t>(it->i  ,it->j+1,it->k  )))
				un->getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k  ),0)
					+= it->c[L][R][L] * sqr(cfg.dt) * it->val.at(step)
					* cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k  ),0);

			if (un->isPresent(Dim3D<int_t>(it->i+1,it->j+1,it->k  )))
				un->getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k  ),0)
					+= it->c[R][R][L] * sqr(cfg.dt) * it->val.at(step)
					* cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k  ),0);

			if (cfg.dims == 3) {
				if (un->isPresent(Dim3D<int_t>(it->i  ,it->j  ,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k+1),0)
						+= it->c[L][L][R] * sqr(cfg.dt) * it->val.at(step)
						* cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k+1),0);

				if (un->isPresent(Dim3D<int_t>(it->i+1,it->j  ,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k+1),0)
						+= it->c[R][L][R] * sqr(cfg.dt) * it->val.at(step)
						* cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k+1),0);

				if (un->isPresent(Dim3D<int_t>(it->i  ,it->j+1,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k+1),0)
						+= it->c[L][R][R] * sqr(cfg.dt) * it->val.at(step)
						* cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k+1),0);

				if (un->isPresent(Dim3D<int_t>(it->i+1,it->j+1,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k+1),0)
						+= it->c[R][R][R] * sqr(cfg.dt) * it->val.at(step)
						* cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k+1),0);
			}
		}

		swap(un, u);

		t += cfg.dt;
	}

#ifdef USE_MPI
	MPI_File_close(&fh);
#else
	if (rgmpi::worldRank() == 0) {
		recvsout.close();
	}
#endif
	}
	rgmpi::forceFinalize();
}
