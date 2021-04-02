#!/bin/sh

n_thr=4
n_itr=1000

if [ -n "$1" ]; then
    n_thr=$1
fi

if [ -n "$2" ]; then
    n_itr=$2
fi

build_dir=../build
debug_dir=$build_dir/debug

if [ ! -d $build_dir ]
then
    echo "No Cmake environment for current project, set one:"
    echo "Make ../build directory and run from there cmake -G [-your environment] .."
    exit
fi

if [ ! -d $debug_dir ]
then
    cd $build_dir
    echo "Binaries not built, building ... "
    cmake --build . --config Debug --target ALL_BUILD
    echo "Binaries were built!"
    cd ../src
fi

cd $debug_dir/bin
echo "started, check logs for optimization details..."
OMP_NUM_THREADS=$n_thr ./to_optimize $n_itr
echo "optimization finished! Preparing plot data..."
./to_plot
echo "Plot data Completed! Executing python script..."
cd ../../../src

pipenv run python3 graph.py || python3 graph.py

echo "Program finished! Check output files in current directory!"

