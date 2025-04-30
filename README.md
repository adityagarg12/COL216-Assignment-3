# COL216-Assignment-3
## How to Compile and Run

### Compilation
To compile the code, use the `make` command in the terminal:
```bash
make
```
This will generate an executable file named `L1simulate`.

### Generating the Report
To generate the `report.pdf` document, use the `make report` command:
```bash
make report
```
Ensure that you have `texlive` installed on your system. You can install it using the following commands:
```bash
sudo apt update
sudo apt install texlive-latex-base
sudo apt install texlive-latex-extra
```

### Cleaning Up
To remove the executable and the generated `report.pdf` document, use the `make clean` command:
```bash
make clean
```

### Running the Executable
You can run the executable using the following command:
```bash
./L1simulate [OPTIONS]
```

### Usage
The program accepts the following options:
- `-t <tracefile>`: Name of the parallel application (e.g., `app1`) whose 4 traces are to be used.
- `-s <s>`: Number of set index bits (number of sets = 2^s).
- `-E <E>`: Associativity (number of cache lines per set).
- `-b <b>`: Number of block bits (block size = 2^b bytes).
- `-o <outfile>`: Logs output in a file for plotting, etc.
- `-h`: Prints this help message.

### Example
Here is an example of how to run the program:
```bash
./L1simulate -t app1 -s 5 -E 2 -b 5 -o output.txt
```