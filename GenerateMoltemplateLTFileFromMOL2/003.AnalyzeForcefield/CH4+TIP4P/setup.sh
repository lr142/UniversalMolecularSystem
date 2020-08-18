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

    for i in $(ls configs)
    do
	cd $i
	echo $i
	lmp run.in
	cd ..
    done
}

prepall(){
    
    for i in $(ls configs)
    do
	prep $i ../system.lt ../configs/$i/ch4h2o.xyz
    done;
}

prepall

# Run them
run
