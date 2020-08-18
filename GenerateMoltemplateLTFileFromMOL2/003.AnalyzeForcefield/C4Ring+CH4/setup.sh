prep(){
    if [ ! -d $1 ]
    then
	mkdir $1
    fi
    cd $1
    echo $1
    cp $2 ./system.lt
    pwd
    moltemplate.sh -atomstyle full system.lt -xyz $3
    
    rm system.in
    rm -rf output_ttree

    cp ../run.in ./
    cd ..
}

# Run them
lmp(){
    mpirun -np 4 ~/Applications/lammps-3Mar20/src/lmp_mpi -in $1
}

run(){
    cd 00.c4r
    lmp run.in
    cd ..
    cd 01.ch4
    lmp run.in
    cd ..
    
    for i in $(ls configs)
    do
	cd $i
	echo $i
	lmp run.in
	cd ..
    done
}

prepall(){
    
    prep 00.c4r ../only.c4r.lt ../c4r.xyz
    prep 01.ch4 ../only.ch4.lt ../ch4.xyz
    
    for i in $(ls configs)
    do
	prep $i ../system.lt ../configs/$i/fs.xyz
    done;
}

#prepall

# Run them
run
