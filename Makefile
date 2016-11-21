simple: acfdtd.cpp config.cpp common.h config.h
	g++ -fno-inline-small-functions acfdtd.cpp config.cpp -g3 ~/code/rgrid/build/libRGridStatic.a -I ~/code/rgrid -W -Wall -o acfdtd

mpi: acfdtd.cpp config.cpp common.h config.h
	mpicxx -fno-inline-small-functions -DUSE_MPI acfdtd.cpp config.cpp -g3 ~/code/rgrid/build_mpi/libRGridStatic.a -I ~/code/rgrid -W -Wall -o acfdtd

simple_opt: acfdtd.cpp config.cpp common.h config.h
	g++ -O2 acfdtd.cpp config.cpp -g ~/code/rgrid/build/libRGridStatic.a -I ~/code/rgrid -W -Wall -o acfdtd

mpi_opt: acfdtd.cpp config.cpp common.h config.h
	mpicxx -O2 -DUSE_MPI acfdtd.cpp config.cpp -g ~/code/rgrid/build_mpi/libRGridStatic.a -I ~/code/rgrid -W -Wall -o acfdtd
