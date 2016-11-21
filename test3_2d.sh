#!/bin/bash

# Свободная граница с одной стороны (free=y в madagascar).
# !!! Тут нет центральной точки, т.к. количество узлов четное. (100,100) чуть дальше от начала координат.
# В мадагаскаре просто зануляют первую не ghost линию узлов, а мы - нет. Поэтому у нас на 1 узел с той стороны меньше (из-за этого сдвигаются все координаты источников)

source mad-source.sh

mkdir -p test
cd test

sfspike > Fvel.rsf mag=1500 nsp=1 k1=1 l1=200 d1=4 d2=4 \
	label1=z label2=x n1=200 n2=200 o1=2 o2=2 unit1=m unit2=m
sfspike > Fden.rsf mag=1 nsp=1 k1=1 l1=200 d1=4 d2=4 label1=z \
	label2=x n1=200 n2=200 o1=2 o2=2 unit1=m unit2=m
sfspike n1=2 nsp=2 k1=1,2 mag=402,402 o1=0 o2=0 > Fsou.rsf
sfspike n1=2 nsp=2 k1=1,2 mag=402,202 o1=0 o2=0 > Frec.rsf
sfspike nsp=1 n1=2000 d1=0.0005 k1=400 | sfricker1 frequency=20 | sftransp > Fwav.rsf
sfawefd2d < Fwav.rsf vel=Fvel.rsf sou=Fsou.rsf rec=Frec.rsf wfl=Fwfl.rsf \
	den=Fden.rsf > Fdat.rsf verb=y free=y fdorder=2 expl=y snap=y dabc=y jdata=1 jsnap=10
echo 'label1=z unit1=m label2=x unit2=m' >> Fwfl.rsf

cat > config.txt <<END
2
200 199
4 4
2000
0.0005
10
400 10
1 1 0 1
4 2
3 5
srcs.txt
rcvs.txt
K.txt
rhox.txt rhoy.txt
END

cat > srcs.txt <<END
src_func.txt 100 100
END

cat > rcvs.txt <<END
100 50
END

rm -f K.txt
cat > datagenconf.txt <<END
K.txt 2250000 200 200 1
END
../datagen datagenconf.txt

rm -f rhox.txt rhoy.txt

cat > datagenconf.txt <<END
rhox.txt 1 201 200 1
END
../datagen datagenconf.txt

cat > datagenconf.txt <<END
rhoy.txt 1 200 201 1
END
../datagen datagenconf.txt

sfdd form=ascii out=src_func.txt < Fwav.rsf > /dev/null
sfdd form=ascii out=rcvs_mad < Fdat.rsf > /dev/null

tr -d '\n' < K.txt > __tmp
tr ' ' '\n' < __tmp > K.txt

tr -d '\n' < rho.txt > __tmp
tr ' ' '\n' < __tmp > rho.txt

tr -d '\n' < src_func.txt > __tmp
tr ' ' '\n' < __tmp > src_func.txt

tr -d '\n' < rcvs_mad > __tmp
tr ' ' '\n' < __tmp > rcvs_mad

rm __tmp

time mpirun -np 8 ../acfdtd config.txt #> /dev/null
echo OK
