#!/bin/bash

source mad-source.sh

mkdir -p test
cd test

cat > rec.txt <<END
50 100 100
END

cat > rec-mad.txt <<END
320 200
360 200
400 200
440 200
480 200
520 200
END

rm -f rcvs_mad

sfspike > Fvel.rsf mag=1500 nsp=1 k1=1 l1=200 d1=4 d2=4 \
	label1=z label2=x n1=200 n2=200 o1=2 o2=2 unit1=m unit2=m
sfspike > Fden.rsf mag=1 nsp=1 k1=1 l1=200 d1=4 d2=4 label1=z \
	label2=x n1=200 n2=200 o1=2 o2=2 unit1=m unit2=m
sfspike n1=2 nsp=2 k1=1,2 mag=400,400 o1=0 o2=0 > Fsou.rsf
echo in=rec-mad.txt n1=2 n2=6 data_format=ascii_float | sfdd form=native > Frec.rsf
#sfspike n1=2 nsp=2 k1=1,2 mag=402,202 o1=0 o2=0 > Frec.rsf
sfspike nsp=1 n1=2000 d1=0.0005 k1=400 | sfricker1 frequency=20 | sftransp > Fwav.rsf
sfawefd2d < Fwav.rsf vel=Fvel.rsf sou=Fsou.rsf rec=Frec.rsf wfl=Fwfl.rsf \
	den=Fden.rsf > Fdat.rsf verb=y free=n fdorder=2 expl=y snap=y dabc=n jdata=1 jsnap=10
echo 'label1=z unit1=m label2=x unit2=m' >> Fwfl.rsf

cat > config.txt <<END
dimensions = 3
order = 2
nodes = 50, 50, 50
origin = 0.0, 0.0, 0.0
space_step = 4.0, 4.0, 4.0
time_steps = 2000
time_step = 0.0005
save_every = 10
save_format = bin
pml_max = 400 # значение PML на границе (растет квадратично)
pml_nodes = 10 # количество узлов, отведенных под PML
left_boundaries = free, free, free # граничные условия на левых границах (free/free)
right_boundaries = free, free, free # граничные условия на правых границах
global_parts = 1, 1, 1
local_parts = 1, 1, 1
sources_position = srcs.txt
receivers_position = rec.txt
bulk_modulus = K.txt
density = rhox.txt, rhoy.txt, rhoz.txt
receivers_output = rcvs.txt
END

cat > srcs.txt <<END
src_func.txt 100 100 100
END

rm -f K.txt
cat > datagenconf.txt <<END
K.txt 2250000 50 50 50
END
../datagen datagenconf.txt

rm -f rhox.txt rhoy.txt rhoz.txt

cat > datagenconf.txt <<END
rhox.txt 1 49 50 50
END
../datagen datagenconf.txt

cat > datagenconf.txt <<END
rhoy.txt 1 50 49 50
END
../datagen datagenconf.txt

cat > datagenconf.txt <<END
rhoz.txt 1 50 50 49
END
../datagen datagenconf.txt

sfdd form=ascii out=src_func.txt < Fwav.rsf > /dev/null
sfdisfil > rcvs_mad col=6 format="%e " number=n < Fdat.rsf

tr -d '\n' < K.txt > __tmp
tr ' ' '\n' < __tmp > K.txt

tr -d '\n' < rho.txt > __tmp
tr ' ' '\n' < __tmp > rho.txt

tr -d '\n' < src_func.txt > __tmp
tr ' ' '\n' < __tmp > src_func.txt

rm __tmp

rm -f rcvs.asc
rm -f rcvs.txt

#time mpirun -np 1 ../acfdtd config.txt #> /dev/null
../acfdtd config.txt
echo OK

../btoa <rcvs.txt >rcvs.asc 1

../plot
