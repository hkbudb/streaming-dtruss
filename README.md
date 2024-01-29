
## Description
Welcome to the readme file of the code for the paper titled **Truss-based Community Search over Streaming Directed Graphs**. In this work, we provide a formal definition of the problem of D-truss-based community search over streaming directed graphs. Using the sliding window model as a basis, we initially introduce a peeling-based algorithm along with three optimizations to reduce the peeling input. Additionally, we propose an order-based algorithm to avoid retrieving the community from scratch.

## Contributions
- Novel problem of D-truss-based community search over streaming graphs
- A peeling-based algorithm that incorporates three optimizations
- An order-based algorithm

## Environment
All of our algorithms are implemented in C++ 11 and compiled with the g++ compiler at -O3 optimization level. The development environment used for implementing and testing is:

-Linux version: Oracle Linux 8.6

-g++ version: 11.2.0


## Dataset format

The input is a .txt file that has three columns in the form of "u v t". An example is shown below:

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

## How to compile the algorithms

To compile the corresponding code, first go to the directory of corresponding algorithms, then:

```bash
./compile.sh
```
After compilation, executable files named ```order```/```repeel-all```/```repeel-bfs```/```repeel-hindex```/```repeel-prediction```/```repeeling```/ will be generated.


## How to run the algorithms

The running commands for running the algorithms are shown below. Note that if you want to run the algorithm, you should go to the corresponding directory first. E.g., to run the repeel without OPT, you should be in the "repeel" folder. Please note that the query vertices Q are in the form of string: "v1 v2 v3 ...".

```bash
For the repeel without OPT: repeeling filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with OPT1: repeel-hindex filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with OPT2: repeel-bfs filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with OP3: repeel-prediction filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with all OPTs: repeel-all filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the order-based method: order filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
```

## 



