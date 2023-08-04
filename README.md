# Streaming-Dtruss
This is the code of the paper titled as "Truss-based Community Search over Streaming Directed Graphs". 
The input is the edge list, each row is in the form of "from-node to-node time-stamp". 
To compile the code, go to the folder and run the compile.sh.

The running commands for running the algorithms are shown below, note that the command should be run in the corresponding folder, e.g., to run the repeel without OPT, you should be in the "repeel" folder

For the repeel without OPT: repeeling filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf

For the repeel with OPT1: repeel-hindex filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf

For the repeel with OPT2: bfsbased filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf

For the repeel with OP3: prediction filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf

For the repeel with all OPTs: repeel-all filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf

For the order-based method: order filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf

Welcome to the readme file of the code for **Truss-based Community Search over Streaming Directed Graphs**. In this work, several efficient algorithms are proposed to answer the community search query over streaming directed graphs. All of our work are implemented in c++.

## Dependencies
The algorithms should also work on Unix-like distributions. Building our algorithms requires the following software installed as dependencies.

- A modern C++ compiler compliant with C++-11 standard. (g++ >= 4.8.1 or clang++ >= 3.3)


## The input

The input is a .txt file that has three columns, in the form of "u v t", an example is shown below:

```bash
1 2 1
3 4 2
5 2 3
6 7 4
8 7 5
9 10 6
9 11 7
12 13 8
9 14 9
9 15 10
```

## Compiling the algorithms

Go to the directory of corresponding algorithms and simply run ./compile.sh


## How to run the algorithms

If you want to run the algorithm, you should go to the corresponding directory first, e.g., to run the repeel without OPT, you should be in the "repeel" folder

```bash
For the repeel without OPT: repeeling filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf
For the repeel with OPT1: repeel-hindex filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf
For the repeel with OPT2: bfsbased filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf
For the repeel with OP3: prediction filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf
For the repeel with all OPTs: repeel-all filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf
For the order-based method: order filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf
```



Here is the readme file of the implementation for the algorithms in **Truss-based Community Search over Streaming Directed Graphs**. All the codes are implemented in c++:

**libgrape-lite** is a C++ library from Alibaba for parallel graph processing. It differs from prior systems in its ability to parallelize sequential graph algorithms as a whole by following the *PIE* programming model from [GRAPE](https://dl.acm.org/doi/10.1145/3035918.3035942). Sequential algorithms can be easily ["plugged into"](examples/analytical_apps/sssp/sssp_auto.h) libgrape-lite with only minor changes and get parallelized to handle large graphs efficiently. In addition to the ease of programming, libgrape-lite is designed to be highly [efficient](Performance.md) and [flexible](examples/gnn_sampler), to cope the scale, variety and complexity from real-life graph applications.

## Building **libgrape-lite**

### Dependencies
**libgrape-lite** is developed and tested on CentOS 7. It should also work on other unix-like distributions. Building libgrape-lite requires the following softwares installed as dependencies.

- [CMake](https://cmake.org/) (>=2.8)
- A modern C++ compiler compliant with C++-11 standard. (g++ >= 4.8.1 or clang++ >= 3.3)
- [MPICH](https://www.mpich.org/) (>= 2.1.4) or [OpenMPI](https://www.open-mpi.org/) (>= 3.0.0)
- [glog](https://github.com/google/glog) (>= 0.3.4)


Here are the dependencies for optional features:
- [jemalloc](http://jemalloc.net/) (>= 5.0.0) for better memory allocation;
- [Doxygen](https://www.doxygen.nl/index.html) (>= 1.8) for generating documentation;
- Linux [HUGE_PAGES](http://www.kernel.org/doc/Documentation/vm/hugetlbpage.txt) support, for better performance.

Extra dependencies are required by examples:
- [gflags](https://github.com/gflags/gflags) (>= 2.2.0);
- [Apache Kafka](https://github.com/apache/kafka) (>= 2.3.0);
- [librdkafka](https://github.com/edenhill/librdkafka)(>= 0.11.3);


### Building libgrape-lite and examples

Once the required dependencies have been installed, go to the root directory of libgrape-lite and do a out-of-source build using CMake.

```bash
mkdir build && cd build
cmake ..
make -j
```

The building targets include a shared/static library, and two sets of examples: [analytical_apps](./examples/analytical_apps) and a [gnn_sampler](./examples/gnn_sampler).

Alternatively, you can build a particular target with command:

```bash
make libgrape-lite # or
make analytical_apps # or
make gnn_sampler
```

## Running libgrape-lite applications

### Graph format

The input of libgrape-lite is formatted following the [LDBC Graph Analytics](http://graphalytics.org) benchmark, with two files for each graph, a `.v` file for vertices with 1 or 2 columns, which are a vertex_id and optionally followed by the data assigned to the vertex; and a `.e` file for edges with 2 or 3 columns, representing source, destination and optionally the data on the edge, correspondingly. See sample files `p2p-31.v` and `p2p-31.e` under the [dataset](dataset/) directory. 

### Example applications
To run the our algorithms, users may use command like this:

```bash
#All the scripts below should add the item to specify partition strategy as follows
#使用metis strategy，在GRAPE来处理图(rfile 可以是其他partition 算法的处理结果)：
mpirun  -n 4 ./run_app --vfile ../dataset/message.v --efile ../dataset/message.e --rfile ../dataset/message.r4 --application dcorefirst --partition metis --out_prefix ./output_message_step1 --directed

#使用grape自带的partition处理图(seg/hash)：
mpirun  -n 4 ./run_app --vfile ../dataset/message.v --efile ../dataset/message.e --application dcorefirst --partition seg/hash --out_prefix ./output_message_step1 --directed



# run phase1 (vertex-centric implementation) of Distributed Anchored Coreness Computation.
mpirun -hostfile ../hosts -n 4 ./run_app --vfile ../dataset/sample_dcore.v --efile ../dataset/sample_dcore.e --application dcorefirst --out_prefix ./output_sample_step1 --directed
#cat the output of phase1 to sample_step1.v, which will used as the input of phase2 
cat ./output_sample_step1/*frag* &> ./output_sample_step1/sample_step1.v"
# run phase2 (vertex-centric implementation) of Distributed Anchored Coreness Computation.
mpirun -hostfile ../hosts -n 4 ./run_app --vfile ./output_sample_step1/sample_step1.v --efile ../dataset/sample_dcore.e --application dcoresecond --out_prefix ./output_sample_step2 --directed
#cat the output of phase2 to sample_step2.v, which will used as the input of phase3 
cat ./output_sample_step2/*frag* &> ./output_sample_step1/sample_step2.v"
# run phase3 (vertex-centric implementation) of Distributed Anchored Coreness Computation.
mpirun -hostfile ../hosts -n 4 ./run_app --vfile ./output_sample_step2/sample_step1.2 --efile ../dataset/sample_dcore.e --application dcorethird --out_prefix ./output_sample_step3 --directed

#To run the distributed skyline coreness computation algorithm, you have to initilize the input as k_{max}.l_{max}(k_{max},v)
#To get the l_{max}(k_{max},v):
mpirun -hostfile ../hosts -n 4 ./run_app --vfile ../dataset/sample_dcore.v --efile ../dataset/sample_dcore.e --application getoutcoreness --out_prefix ./output_outcoreness --directed
#use the output of dcorefirst and getoutcoreness to  initialize the input as sample_SH.v:
mpirun -hostfile ../hosts -n 4 ./run_app --vfile ../dataset/sample_SH.v --efile ../dataset/sample_dcore.e --application dcoreoptimized --out_prefix ./output_sample_optimized --directed

#for the block implementation, you may use similar scripts to run them.
# HOSTFILE provides a list of hosts where MPI processes are launched. 


# see more flags info.
./run_app --help
```


## Documentation

Documentation is generated using Doxygen. Users can build doxygen documentation in the build directory using:

```bash
cd build
make doc
# open docs/index.html
```

The latest version of online documentation can be found at [https://alibaba.github.io/libgrape-lite](https://alibaba.github.io/libgrape-lite)

## License

**libgrape-lite** is distributed under [Apache License 2.0](./LICENSE). Please note that third-party libraries may not have the same license as libgrape-lite.


