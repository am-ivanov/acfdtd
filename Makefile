CXXFLAGS = -I . -I ~/code/rgrid -W -Wall -Wextra -Werror=narrowing

simple: acfdtd.cpp config.cpp common.h config.h
	g++ $(CXXFLAGS) -fno-inline-small-functions $(wildcard kutils/*.cpp) acfdtd.cpp config.cpp -g3 ~/code/rgrid/build/libRGridStatic.a -o acfdtd

mpi: acfdtd.cpp config.cpp common.h config.h
	mpicxx $(CXXFLAGS) -fno-inline-small-functions -DUSE_MPI $(wildcard kutils/*.cpp) acfdtd.cpp config.cpp -g3 ~/code/rgrid/build_mpi/libRGridStatic.a -o acfdtd -fmax-errors=1

simple_opt: acfdtd.cpp config.cpp common.h config.h
	g++ $(CXXFLAGS) -O2 $(wildcard kutils/*.cpp) acfdtd.cpp config.cpp -g ~/code/rgrid/build/libRGridStatic.a -o acfdtd

mpi_opt: acfdtd.cpp config.cpp common.h config.h
	mpicxx $(CXXFLAGS) -O2 -DUSE_MPI $(wildcard kutils/*.cpp) acfdtd.cpp config.cpp -g ~/code/rgrid/build_mpi/libRGridStatic.a -o acfdtd

tools:
	g++ -O2 -o datagen datagen.cpp
	g++ -O2 -o btoa bintoascii.cpp
