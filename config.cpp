#include "config.h"
#include "rgrid/darrayscatter.h"

#include <fstream>

using namespace rgrid;

using namespace std;

void Config::readConfig(string file) {
#ifdef USE_MPI
	MPI_File fh;
	MPI_Offset filesize;
	if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
		throw logic_error("(MPI) Couldn't open \"" + file + "\"");
	if (0 != MPI_File_get_size(fh, &filesize))
		throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
	vector<char> buf(filesize);
	MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);

	std::istringstream ss(std::string(&buf[0],static_cast<size_t>(filesize)));
	istream& is = ss;
#else
	fstream fs;
	fs.open(file.c_str(), ios_base::in);
	if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
	istream& is = fs;
#endif
	is >> dims;

	if (dims == 2) {
		is >> nx >> ny;
		is >> dx >> dy;
		nz = 1;
		dz = 1;
	} else if (dims == 3) {
		is >> nx >> ny >> nz;
		is >> dx >> dy >> dz;
	} else throw logic_error("Wrong dimensions");

	is >> steps;
	is >> dt;

	is >> saveStep;

	is >> max_pml >> pml_len;

	is >> isPml[0][0] >> isPml[0][1] >> isPml[1][0] >> isPml[1][1];
	if (dims == 3) {
		is >> isPml[2][0] >> isPml[2][1];
	} else {
		isPml[2][0] = false;
		isPml[2][1] = false;
	}

	is >> gx >> gy;
	if (dims == 3)
		is >> gz;
	else
		gz = 1;
	is >> lx >> ly;
	if (dims == 3)
		is >> lz;
	else
		lz = 1;

	string sourcesFile;
	is >> sourcesFile;

	string receiversFile;
	is >> receiversFile;

	string KFile;
	is >> KFile;

	string rhoXFile;
	is >> rhoXFile;

	string rhoYFile;
	is >> rhoYFile;

	string rhoZFile;
	if (dims == 3) {
		is >> rhoZFile;
	}
#ifndef USE_MPI
	fs.close();
#endif
	readSources(sourcesFile);
	readReceivers(receiversFile);
	readK(KFile);
	readRhoX(rhoXFile);
	readRhoY(rhoYFile);
	if (dims == 3) {
		readRhoZ(rhoZFile);
	} else {
		rhoz.setSizes(
			Dim3D<int>(nx,ny,nz),
			Dim3D<int>(gx,gy,gz),
			Dim3D<int>(lx,ly,lz),
			Dim3D<int>(0,0,0),
			1);
	}
}

void Config::readSources(string file) {
	std::istringstream ss;
#ifdef USE_MPI
	{
		MPI_File fh;
		MPI_Offset filesize;
		if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
			throw logic_error("(MPI) Couldn't open \"" + file + "\"");
		if (0 != MPI_File_get_size(fh, &filesize))
			throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
		vector<char> buf(filesize);
		MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		ss.str(std::string(&buf[0],static_cast<size_t>(filesize)));
	}
#else
	{
		fstream fs;
		fs.open(file.c_str(), ios_base::in);
		if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
		std::string contents((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
		fs.close();
		ss.str(contents);
	}
#endif
	while (1) {
		Source e;
		string file_src;
		ss >> file_src >> e.x >> e.y;
		if (dims == 3)
			ss >> e.z;
		else
			e.z = 0;
		if (ss.eof()) break;
		istringstream ss_src;
#ifdef USE_MPI
	{
		string& file = file_src;
		MPI_File fh;
		MPI_Offset filesize;
		if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
			throw logic_error("(MPI) Couldn't open \"" + file + "\"");
		if (0 != MPI_File_get_size(fh, &filesize))
			throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
		vector<char> buf(filesize);
		MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		ss_src.str(std::string(&buf[0],static_cast<size_t>(filesize)));
	}
#else
	{
		string& file = file_src;
		fstream fs;
		fs.open(file.c_str(), ios_base::in);
		if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
		std::string contents((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
		fs.close();
		ss_src.str(contents);
	}
#endif
		for (int i = 0; i != steps; ++i) {
			double val;
			ss_src >> val;
			e.val.push_back(val);
		}
		src.push_back(e);
	}
}

void Config::readReceivers(string file) {
	std::istringstream ss;
#ifdef USE_MPI
	{
		MPI_File fh;
		MPI_Offset filesize;
		if (0 != MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
			throw logic_error("(MPI) Couldn't open \"" + file + "\"");
		if (0 != MPI_File_get_size(fh, &filesize))
			throw logic_error("(MPI) Couldn't get \"" + file + "\" size");
		vector<char> buf(filesize);
		MPI_File_read_all(fh, &buf[0], filesize, MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);
		ss.str(std::string(&buf[0],static_cast<size_t>(filesize)));
	}
#else
	{
		fstream fs;
		fs.open(file.c_str(), ios_base::in);
		if (!fs.is_open()) throw logic_error("Couldn't open \"" + file + "\"");
		std::string contents((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
		fs.close();
		ss.str(contents);
	}
#endif
	while (1) {
		Receiver e;
		ss >> e.x >> e.y;
		if (dims == 3)
			ss >> e.z;
		else
			e.z = 0;
		if (ss.eof()) break;
		rcv.push_back(e);
	}
}

void Config::readK(string file) {
	K.setSizes(
		Dim3D<int>(nx, ny, nz),
		Dim3D<int>(gx, gy, gz),
		Dim3D<int>(lx, ly, lz),
		Dim3D<int>(0, 0, 0),
		1);
	K.loadData(file);
}

void Config::readRhoX(string file) {
	// last block in global partitioning takes one additional node
	Dim3D<vector<int> > globalWidth = K.getWidth();
	globalWidth[X].at(globalWidth[X].size()-1) += 1;
	Dim3D<vector<int> > localWidth = K.getLocalContainer().getWidth();
	// so if last block on current process we add node here
	if (K.getInternalPos(X) == gx - 1) {
		localWidth[X].at(localWidth[X].size()-1) += 1;
	}
	rhox.setSizes(
		globalWidth,
		localWidth,
		Dim3D<int>(1, 1, dims == 3 ? 1 : 0),
		1);
	rhox.loadData(file);
	rhox.externalSyncStart();
	rhox.internalSync();
	rhox.externalSyncEnd();
}

void Config::readRhoY(string file) {
	Dim3D<vector<int> > globalWidth = K.getWidth();
	globalWidth[Y].at(globalWidth[Y].size()-1) += 1;
	Dim3D<vector<int> > localWidth = K.getLocalContainer().getWidth();
	if (K.getInternalPos(Y) == gy - 1) {
		localWidth[Y].at(localWidth[Y].size()-1) += 1;
	}
	rhoy.setSizes(
		globalWidth,
		localWidth,
		Dim3D<int>(1, 1, dims == 3 ? 1 : 0),
		1);
	rhoy.loadData(file);
	rhoy.externalSyncStart();
	rhoy.internalSync();
	rhoy.externalSyncEnd();
}

void Config::readRhoZ(string file) {
	Dim3D<vector<int> > globalWidth = K.getWidth();
	globalWidth[Z].at(globalWidth[Z].size()-1) += 1;
	Dim3D<vector<int> > localWidth = K.getLocalContainer().getWidth();
	if (K.getInternalPos(Z) == gz - 1) {
		localWidth[Z].at(localWidth[Z].size()-1) += 1;
	}
	rhoz.setSizes(
		globalWidth,
		localWidth,
		Dim3D<int>(1, 1, dims == 3 ? 1 : 0),
		1);
	rhoz.loadData(file);
	rhoz.externalSyncStart();
	rhoz.internalSync();
	rhoz.externalSyncEnd();
}
