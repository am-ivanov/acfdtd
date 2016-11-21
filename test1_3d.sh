#!/bin/bash

# Отключены все поглощающие граничные условия (dabc=n в madagascar и 6 нулей в конфиге у меня).
# Центральная точка - 25, 25, 25

source mad-source.sh

mkdir -p test
cd test

sfspike > Fvel.rsf mag=1500 nsp=1 k1=1 l1=51 d1=4 d2=4 d3=4\
	label1=z label2=y label3=x n1=51 n2=51 n3=51 o1=2 o2=2 o3=2 unit1=m unit2=m
sfspike > Fden.rsf mag=1 nsp=1 k1=1 l1=51 d1=4 d2=4 d3=4 \
	label1=z label2=y label3=x n1=51 n2=51 n3=51 o1=2 o2=2 o3=2 unit1=m unit2=m
sfspike n1=3 nsp=3 k1=1,2,3 mag=102,102,102 o1=0 o2=0 o3=0 > Fsou.rsf
sfspike n1=3 nsp=3 k1=1,2,3 mag=102,102,62 o1=0 o2=0 o3=0 > Frec.rsf
sfspike nsp=1 n1=2000 d1=0.0005 k1=400 | sfricker1 frequency=20 | sftransp > Fwav.rsf
sfawefd3d < Fwav.rsf vel=Fvel.rsf sou=Fsou.rsf rec=Frec.rsf wfl=Fwfl.rsf \
	den=Fden.rsf > Fdat.rsf verb=y free=n fdorder=2 expl=y snap=y dabc=n jdata=1 jsnap=10
echo 'label1=z unit1=m label2=y unit2=m label3=x unit3=m' >> Fwfl.rsf

cat > config.txt <<END
3
51 51 51
4 4 4
2000
0.0005
10
400 10
0 0 0 0 0 0
2 1 2
1 1 1
srcs.txt
rcvs.txt
K.txt
rhox.txt rhoy.txt rhoz.txt
END

cat > srcs.txt <<END
src_func.txt 25 25 25
END

cat > rcvs.txt <<END
25 25 15
END

rm -f K.txt
cat > datagenconf.txt <<END
K.txt 2250000 51 51 51
END
../datagen datagenconf.txt

rm -f rhox.txt rhoy.txt rhoz.txt

cat > datagenconf.txt <<END
rhox.txt 1 52 51 51
END
../datagen datagenconf.txt

cat > datagenconf.txt <<END
rhoy.txt 1 51 52 51
END
../datagen datagenconf.txt

cat > datagenconf.txt <<END
rhoz.txt 1 51 51 52
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

time mpirun -np 4 ../acfdtd config.txt #> /dev/null
echo OK

../plot
