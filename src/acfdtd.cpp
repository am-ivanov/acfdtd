#include "config.h"

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "rgrid/vtksaver.h"
#include "rgrid/darrayscatter.h"

#include "rgrid/gridoperator.h"
#include "rgrid/nodeoperator.h"
#include "rgrid/range.h"

#include <fenv.h>

using namespace std;
using namespace rgrid;
using namespace operators;

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

	//for (int i = -cfg.ho; i <= cfg.ho; ++i) {
	//	if (i != 0)
	//		cout << cfg.fdc.sc1(i) << " ";
	//}
	//cout << endl;

	if (rgmpi::worldRank() == 0) {
		cout << "dimensions = " << cfg.dims << endl;
		cout << "order = " << cfg.ho * 2 << endl;
		cout << "nodes = " << cfg.nx << ", " << cfg.ny << ", " << cfg.nz << endl;
		cout << "origin = " << cfg.ox << ", " << cfg.oy << ", " << cfg.oz << endl;
		cout << "space_step = " << cfg.dx << ", " << cfg.dy << ", " << cfg.dz << endl;
		cout << "time_steps = " << cfg.steps << endl;
		cout << "time_step = " << cfg.dt << endl;
		cout << "save_every = " << cfg.saveStep << endl;
		cout << "pml_max = " << cfg.max_pml << endl;
		cout << "pml_nodes = " << cfg.pml_len << endl;
		cout << "left_boundaries = " << cfg.isPml[0][0] << ", " << cfg.isPml[1][0] << ", " << cfg.isPml[2][0] << endl;
		cout << "right_boundaries = " << cfg.isPml[0][1] << ", " << cfg.isPml[1][1] << ", " << cfg.isPml[2][1] << endl;
		cout << "global_parts = " << cfg.gx << ", " << cfg.gy << ", " << cfg.gz << endl;
		cout << "local_parts = " << cfg.lx << ", " << cfg.ly << ", " << cfg.lz << endl;

		cout << "Sources: " << cfg.src.size() << endl;
		for (vector<Source>::size_type i = 0; i < cfg.src.size(); ++i) {
			cout << "source " << i << " position:" << endl;
			cout
				<< "index: "
				<< cfg.src.at(i).i << " " << cfg.src.at(i).j << " " << cfg.src.at(i).k
				<< ", coord: "
				<< cfg.src.at(i).x << " " << cfg.src.at(i).y << " " << cfg.src.at(i).z
				<< endl;
		}

		cout << "Receivers: " << cfg.src.size() << endl;
		for (vector<Receiver>::size_type i = 0; i < cfg.src.size(); ++i) {
			cout << "receiver " << i << " position:" << endl;
			cout
				<< "index: "
				<< cfg.rcv.at(i).i << " " << cfg.rcv.at(i).j << " " << cfg.rcv.at(i).k
				<< ", coord: "
				<< cfg.rcv.at(i).x << " " << cfg.rcv.at(i).y << " " << cfg.rcv.at(i).z
				<< endl;
		}
	}

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

	Range<int_t> r_a2(
	cfg.pml_len, cfg.pml_len, 0, 
	cfg.nx-cfg.pml_len, cfg.ny-cfg.pml_len, cfg.nz);
	
	Range<int_t> r_a(0, 0, 0, cfg.nx, cfg.ny, cfg.nz);
	
	Range<int_t> r_x(cfg.pml_len, 0, 0, cfg.nx-cfg.pml_len, cfg.ny, cfg.nz);
	Range<int_t> r_x_l(0, 0, 0, cfg.pml_len, cfg.ny, cfg.nz);
	Range<int_t> r_x_r(cfg.nx-cfg.pml_len, 0, 0, cfg.nx, cfg.ny, cfg.nz);
	
	Range<int_t> r_y(0, cfg.pml_len, 0, cfg.nx, cfg.ny-cfg.pml_len, cfg.nz);
	Range<int_t> r_y_l(0, 0, 0, cfg.nx, cfg.pml_len, cfg.nz);
	Range<int_t> r_y_r(0, cfg.ny-cfg.pml_len, 0, cfg.nx, cfg.ny, cfg.nz);
	
	Range<int_t> r_z(0, 0, cfg.pml_len, cfg.nx, cfg.ny, cfg.nz-cfg.pml_len);
	Range<int_t> r_z_l(0, 0, 0, cfg.nx, cfg.ny, cfg.pml_len);
	Range<int_t> r_z_r(0, 0, cfg.nz-cfg.pml_len, cfg.nx, cfg.ny, cfg.nz);

	// operator expressions
	GridOp<real_t, int_t> gop;
	GridOp<real_t, int_t>::VarType var_ux = gop.setVar("ux", x1das);
	GridOp<real_t, int_t>::VarType var_uy = gop.setVar("uy", y1das);
	GridOp<real_t, int_t>::VarType var_k = gop.setVar("k", cfg.K);
	GridOp<real_t, int_t>::VarType var_rhox = gop.setVar("rhox", cfg.rhox);
	GridOp<real_t, int_t>::VarType var_rhoy = gop.setVar("rhoy", cfg.rhoy);
	GridOp<real_t, int_t>::VarType var_rhoz = gop.setVar("rhoz", cfg.rhoz);
	GridOp<real_t, int_t>::ConstType con_t2(cfg.dt*cfg.dt);
	GridOp<real_t, int_t>::ConstType con_2(2);

	real_t t = 0.0;
	for (int_t step = 0; step != cfg.steps; ++step) {
		// rename DArrays because of swapping variables
		GridOp<real_t, int_t>::VarType var_u = gop.setVar("u", *u);
		GridOp<real_t, int_t>::VarType var_un = gop.setVar("un", *un);

		DArrayContainer<real_t, int_t>& dac = u->getLocalContainer();
		//DArrayContainer<real_t, int_t>& dacNext = un->getLocalContainer();

		if (step % cfg.saveStep == 0) {
			if (rgmpi::worldRank() == 0)
				cout << "Step " << step << endl;

			char name[100];

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

			if (cfg.format == 0) {
				// save to VTK
				ds.inverseBytes();
				sprintf(name, "out-%06ld.vtk", (long)step);
				vs.save(std::string(name));
			}
			else {
				// save to binary
				sprintf(name, "out-%06ld.bin", (long)step);
				fstream trunc_fs(name, fstream::out | fstream::trunc);
				trunc_fs.close();
				dasSave.appendData(std::string(name));
			}
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
		
		gop.apply(r_a2, fd1<X>(cfg.dx, cfg.ho, var_u) / var_rhox, x1das);
		gop.apply(r_a2, fd1<Y>(cfg.dy, cfg.ho, var_u) / var_rhoy, y1das);

		x1das.externalSyncStart();
		y1das.externalSyncStart();
		z1das.externalSyncStart();

		x1das.internalSync();
		y1das.internalSync();
		z1das.internalSync();

		x1das.externalSyncEnd();
		y1das.externalSyncEnd();
		z1das.externalSyncEnd();
		
		gop.apply(
			r_a2,
			con_2 * var_u - var_un + var_k * con_t2 * (bd1<X>(cfg.dx, cfg.ho, var_ux) + bd1<Y>(cfg.dy, cfg.ho, var_uy)),
			*un);

		// insert source
		for (vector<Source>::iterator it = cfg.src.begin(); it != cfg.src.end(); ++it) {
			real_t val = sqr(cfg.dt) * it->val.at(step) / (cfg.dx * cfg.dy * cfg.dz);
			// put source only in real no ghost nodes
			if (un->isPresent(Dim3D<int_t>(it->i  ,it->j  ,it->k  )))
				un->getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k  ),0)
					+= it->c[L][L][L] * val * cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k  ),0);

			if (un->isPresent(Dim3D<int_t>(it->i+1,it->j  ,it->k  )))
				un->getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k  ),0)
					+= it->c[R][L][L] * val * cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k  ),0);

			if (un->isPresent(Dim3D<int_t>(it->i  ,it->j+1,it->k  )))
				un->getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k  ),0)
					+= it->c[L][R][L] * val * cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k  ),0);

			if (un->isPresent(Dim3D<int_t>(it->i+1,it->j+1,it->k  )))
				un->getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k  ),0)
					+= it->c[R][R][L] * val * cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k  ),0);

			if (cfg.dims == 3) {
				if (un->isPresent(Dim3D<int_t>(it->i  ,it->j  ,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k+1),0)
						+= it->c[L][L][R] * val * cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j  ,it->k+1),0);

				if (un->isPresent(Dim3D<int_t>(it->i+1,it->j  ,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k+1),0)
						+= it->c[R][L][R] * val * cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j  ,it->k+1),0);

				if (un->isPresent(Dim3D<int_t>(it->i  ,it->j+1,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k+1),0)
						+= it->c[L][R][R] * val * cfg.K.getNode(Dim3D<int_t>(it->i  ,it->j+1,it->k+1),0);

				if (un->isPresent(Dim3D<int_t>(it->i+1,it->j+1,it->k+1)))
					un->getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k+1),0)
						+= it->c[R][R][R] * val * cfg.K.getNode(Dim3D<int_t>(it->i+1,it->j+1,it->k+1),0);
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
