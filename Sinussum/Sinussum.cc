#include <cmath>
#include <iomanip>  // for setprecision
#include <iostream>
#include <string>
#include <vector>
using namespace std;

typedef vector<vector<char>> CharGrid;

// Constants and error messages from useful_stuff.txt
const double EPSIL_DICHO(1e-9);
const double EPSIL_T(1e-13);
const string BAD_SIGNAL("Error: signal type must be SAWTOOTH, SQUARE or TRIANGLE");
const string NBN_TOO_SMALL("Error: the number of sine terms must be greater than 0");
const string TIME_MIN_MAX("Error: time max is not bigger than min");
const string WRONG_TIME_VAL("Error: both time values must belong to [0., 1.]");
const string SIGNAL_MIN_MAX("Error: signal max is not bigger than min");
const string NBL_TOO_SMALL("Error: the number of lines must be greater than 2");
const string NBL_NOT_ODD("Error: the number of lines must be odd");

//
void readAndValidateParameters(string& signalType, int& nbN, double& tmin,
                               double& tmax, double& signalMin, double& signalMax,
                               int& nbL, int& nbC);
double sawtoothWave(double time);
double squareWave(double time);
double triangularWave(double t);
double calculateTheoreticalValue(const string& signalType, double time);
double approximateSawtooth(double t, int nbN);
double approximateSquare(double t, int nbN);
double approximateTriangle(double t, int nbN);
double approximateSignalValue(double t, int nbN, const string& signalType);
void populateGrid(CharGrid& grid, const string& signalType, int nbN, double tmin,
                  double tmax, double signalMin, double signalMax, int nbC);
CharGrid initializeGrid(int nbL, int nbC);
void printDashes(int nbC);
void displayGrid(const CharGrid& grid);
void setTimeInterval(const string& signalType, int nbN, double& tmin, double& tmax);
double findMaximum(const string& signalType, int nbN, double tmin, double tmax);

void print_error(const string& message) {
    cout << message << endl;
    exit(0);
}

int main() {
    string signalType;
    int nbN;
    double tmin, tmax;
    double signalMin, signalMax;
    int nbL, nbC;

    readAndValidateParameters(signalType, nbN, tmin, tmax, signalMin, signalMax, nbL,
                              nbC);

    CharGrid grid = initializeGrid(nbL, nbC);
    populateGrid(grid, signalType, nbN, tmin, tmax, signalMin, signalMax, nbC);
    displayGrid(grid);

    // tache 3
    cout << setprecision(8) << fixed;
    double maxSignal = findMaximum(signalType, nbN, tmin, tmax);
    // cout << "Maximum value of the approximated signal: " << maxSignal << endl;
    cout << maxSignal << endl;

    return 0;
}

void readAndValidateParameters(string& signalType, int& nbN, double& tmin,
                               double& tmax, double& signalMin, double& signalMax,
                               int& nbL, int& nbC) {
    // Read and validate signal type
    cin >> signalType;
    if (signalType != "SAWTOOTH" && signalType != "SQUARE" &&
        signalType != "TRIANGLE") {
        print_error(BAD_SIGNAL);
    }
    // Read and validate number of terms
    cin >> nbN;
    if (nbN <= 0) {
        print_error(NBN_TOO_SMALL);
    }
    // Read and validate time interval
    cin >> tmin >> tmax;
    if (tmax <= tmin) {
        print_error(TIME_MIN_MAX);
    }
    if (tmin < 0.0 || tmax > 1.0) {
        print_error(WRONG_TIME_VAL);
    }

    // Read and validate signal amplitude interval
    cin >> signalMin >> signalMax;
    if (signalMax <= signalMin) {
        print_error(SIGNAL_MIN_MAX);
    }
    // Read and validate number of display lines
    cin >> nbL;
    // cout<<endl;
    if (nbL <= 2) {
        print_error(NBL_TOO_SMALL);
    }
    if (nbL % 2 == 0) {
        print_error(NBL_NOT_ODD);
    }
    // Calculate number of display columns
    nbC = 2 * nbL - 1;
}

double sawtoothWave(double time) {
    const double period = 2.0;
    double normalized_time = fmod(time, period);
    const double slope = 2.0;
    // At the discontinuity, return the midpoint value
    if (normalized_time == 0.0 || normalized_time == 1.0) {
        return 0.0;  // Discontinuity at t=0 and t=1
    } else if (normalized_time < 1.0) {
        return slope * normalized_time - 1.0;  // Upward slope
    } else {
        return slope * (normalized_time - 1.0) - 1.0;  // Continue the upward slope
    }
}

double squareWave(double time) {
    const double period = 1.0;  // Period of the square wave
    double normalized_time = fmod(time, period);

    // Handle the discontinuities explicitly
    if (fabs(normalized_time - 0.0) < EPSIL_T ||
        fabs(normalized_time - 0.5) < EPSIL_T) {
        return 0.0;  // Return 0 at the points of discontinuity
    } else if (normalized_time < 0.5) {
        // High from just after 0 to just before 0.5
        return 1.0;
    } else {
        // Low from just after 0.5 to just before the next discontinuity
        return -1.0;
    }
}

double triangularWave(double time) {
    const double period = 2.0;  // Full period of the wave
    double normalized_time = fmod(time, period);

    if (normalized_time < 1.0) {
        // In the first half of the period, we have an upward slope followed by a
        // downward slope
        if (normalized_time < 0.5) {
            // Upward slope from -1 to 1 for t in [0, 0.5)
            return 4.0 * normalized_time - 1.0;
        } else {
            // Downward slope from 1 to -1 for t in [0.5, 1)
            return -4.0 * normalized_time + 3.0;
        }
    } else {
        // In the second half of the period, the pattern repeats
        if (normalized_time < 1.5) {
            // Upward slope from -1 to 1 for t in [1, 1.5)
            return 4.0 * (normalized_time - 1.0) - 1.0;
        } else {
            // Downward slope from 1 to -1 for t in [1.5, 2)
            return -4.0 * (normalized_time - 1.0) + 3.0;
        }
    }
}

double calculateTheoreticalValue(const string& signalType, double time) {
    if (signalType == "SAWTOOTH") {
        return sawtoothWave(time);
    } else if (signalType == "SQUARE") {
        return squareWave(time);
    } else if (signalType == "TRIANGLE") {
        return triangularWave(time);
    }
    return 0.0;
}

double approximateSawtooth(double t, int nbN) {
    double sum = 0.0;
    for (int k = 1; k <= nbN; ++k) {
        sum += pow(-1, k) * sin(2 * M_PI * k * (t - 0.5)) / k;
    }
    return -2 / M_PI * sum;
}

double approximateSquare(double t, int nbN) {
    double sum = 0.0;
    for (int k = 1; k <= nbN; ++k) {
        sum += sin(2 * M_PI * (2 * k - 1) * t) / (2 * k - 1);
    }
    return 4 / M_PI * sum;
}

double approximateTriangle(double t, int nbN) {
    double sum = 0.0;
    for (int k = 1; k <= nbN; ++k) {
        sum +=
            pow(-1, k) * sin(2 * M_PI * (2 * k - 1) * (t - 0.25)) / pow(2 * k - 1, 2);
    }
    return -8 / (M_PI * M_PI) * sum;
}

// Function to handle the signal approximation based on the type of signal
double approximateSignalValue(double t, int nbN, const string& signalType) {
    if (signalType == "SAWTOOTH") {
        return approximateSawtooth(t, nbN);
    } else if (signalType == "SQUARE") {
        return approximateSquare(t, nbN);
    } else if (signalType == "TRIANGLE") {
        return approximateTriangle(t, nbN);
    } else {
        cerr << "Unknown signal type!" << endl;
        return 0.0;
    }
}

void populateGrid(CharGrid& grid, const string& signalType, int nbN, double tmin,
                  double tmax, double signalMin, double signalMax, int nbC) {
    // Calculate delta_t and delta_s based on provided specifications
    double delta_t = (tmax - tmin) / (nbC - 1);
    double delta_s = (signalMax - signalMin) / (grid.size() - 1);

    // Determine the row index for the time axis
    int timeAxisRow =
        (signalMin < 0.0 && signalMax > 0.0)
            ? grid.size() - 1 - static_cast<int>(round(-signalMin / delta_s))
            : -1;

    for (int j = 0; j < nbC; ++j) {
        double t = tmin + j * delta_t;  // Calculate the time for this column
        double theoreticalValue = calculateTheoreticalValue(signalType, t);
        double approxValue = approximateSignalValue(t, nbN, signalType);
        // Calculate the row index for the theoretical value
        int rowIndexTheoretical =
            grid.size() - 1 -
            static_cast<int>(round((theoreticalValue - signalMin) / delta_s));
        // Calculate the row index for the approximated value
        int rowIndexApprox =
            grid.size() - 1 -
            static_cast<int>(round((approxValue - signalMin) / delta_s));

        // Place the theoretical value marker
        if (rowIndexTheoretical >= 0 && rowIndexTheoretical < grid.size() &&
            grid[rowIndexTheoretical][j] != '*') {
            grid[rowIndexTheoretical][j] = '+';
        }  // Place the approximated value marker
        if (rowIndexApprox >= 0 && rowIndexApprox < grid.size()) {
            grid[rowIndexApprox][j] = '*';
        }  // Place the time axis marker
        if (timeAxisRow != -1 &&
            grid[timeAxisRow][j] == ' ') {  // Only place '.' if the cell is empty and
                                            // the time axis row is valid
            grid[timeAxisRow][j] = '.';
        }
    }
}

CharGrid initializeGrid(int nbL, int nbC) {
    return CharGrid(nbL, vector<char>(nbC, ' '));  // Initialize grid with spaces
}

void printDashes(int nbC) {
    for (int i = 0; i < nbC; ++i) {
        cout << '-';
    }
    cout << endl;
}

void displayGrid(const CharGrid& grid) {
    int nbC = grid[0].size();
    printDashes(nbC);
    for (const auto& row : grid) {
        for (char c : row) {
            cout << c;
        }
        cout << endl;
    }
    printDashes(nbC);
}

double findMaximum(const string& signalType, int nbN, double tmin, double tmax) {
    double searchMin, searchMax;

    if (signalType == "SAWTOOTH") {
        searchMin = 1.0 - (1.0 / (2.0 * nbN + 1.0));
        searchMax = 1.0;
    } else if (signalType == "SQUARE") {
        searchMin = 0.0;
        searchMax = 1.0 / (2.0 * nbN + 1.0);
    } else if (signalType == "TRIANGLE") {
        searchMin = 0.25 - (1.0 / (2.0 * (2.0 * nbN + 1.0)));
        searchMax = 0.25 + (1.0 / (2.0 * (2.0 * nbN + 1.0)));
    }
    if (signalType == "TRIANGLE") {
        // For a symmetric triangle wave, the maximum is at t = 0.5
        return approximateTriangle(0.5, nbN);
    }
    // Dichotomous search for the maximum value
    double left = searchMin;
    double right = searchMax;
    double epsilon = EPSIL_DICHO;
    double mid, fmid, fleft, fright;
    while ((right - left) > epsilon) {
        mid = (left + right) / 2;
        fmid = approximateSignalValue(mid, nbN, signalType);
        fleft = approximateSignalValue(left, nbN, signalType);
        fright = approximateSignalValue(right, nbN, signalType);
        // Adjust the search interval based on where the maximum lies
        if (fmid > fleft && fmid > fright) {  // The maximum is within (left, right)
            left = left + (right - left) / 4;
            right = right - (right - left) / 4;
        } else if (fleft > fmid) {  // The maximum is within (left, mid)
            right = mid;
        } else {
            // The maximum is within (mid, right)
            left = mid;
        }
    }
    return fmid;  // At this point, the maximum value should be approximately at 'mid'
}
