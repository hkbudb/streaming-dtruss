Welcome to the readme file of the code for the paper titled as **Truss-based Community Search over Streaming Directed Graphs**. In this work, several efficient algorithms are proposed to answer the community search query over streaming directed graphs. All of our algorithms are implemented in C++.

## Dependencies
The algorithms should also work on Unix-like systems. Building our algorithms requires the following software installed as dependencies.

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

The running commands for running the algorithms are shown below. Note that if you want to run the algorithm, you should go to the corresponding directory first, e.g., to run the repeel without OPT, you should be in the "repeel" folder. Please note that the query vertices Q is in the form of string: "v1 v2 v3 ...".

```bash
For the repeel without OPT: repeeling filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with OPT1: repeel-hindex filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with OPT2: bfsbased filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with OP3: prediction filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the repeel with all OPTs: repeel-all filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
For the order-based method: order filepath (e.g. ~/data/msg.txt) windowsize stridesize kc kf Q
```

## 



