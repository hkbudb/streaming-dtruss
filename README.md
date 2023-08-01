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
