# Sinussum

This project, developed as part of the ICC course at EPFL, investigates the impact of using a finite number of terms to approximate signals that are theoretically represented as infinite sums of sinusoids. The project focuses on three signal types: sawtooth, square, and triangular waves.[Projet_Sinussum_V1.01.pdf](https://github.com/user-attachments/files/16339506/Projet_Sinussum_V1.01.pdf)


## Features

- **Signal Approximation**: Approximates sawtooth, square, and triangular signals using a specified number of sinusoidal terms.
- **Parameter Input**: Accepts user input for the type of signal, number of terms, and display parameters.
- **Graphical Representation**: Displays both the theoretical and approximated signals graphically in the terminal.
- **Maximum Value Calculation**: Efficiently calculates the maximum value of the approximated signal over a period using a dichotomic search algorithm.

## Getting Started

### Prerequisites

- A C++11 compatible compiler.

### Running the Program

1. Clone the repository:
    ```bash
    git clone https://github.com/Bahey-shalash/Sinussum.git
    ```
    ```
    cd Sinussum
    ```

2. Compile the source code:
    ```bash
    g++ -std=c++11 Sinussum.cc -o proj
    ```
    or
   ```bash
   clang++ -std=c++11 Sinussum.cc -o proj
   ```

4. Create a test input file (e.g., `t01.txt`), which includes the parameters for the program:
    ```
    SQUARE
    1
    0. 1.
    -1.3 1.3
    5
    ```

5. Run the program with the test file:
    ```bash
    ./proj < t01.txt
    ```

6. To redirect the output to a file:
    ```bash
    ./proj < t01.txt > out01.txt
    ```

## Example

Here's an example of how to execute the program and test its output:

- **Input file (`t01.txt`):**
    ```
    SQUARE
    1
    0. 1.
    -1.3 1.3
    5
    ```

- **Running the program:**
    ```bash
    ./proj < t01.txt > out01.txt
    ```

- **Output file (`out01.txt`):**
    ```
    ---------
     +*+     
     * *     
    *...*...*
         * * 
         +*+ 
    ---------
    1.27323954
    ```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
