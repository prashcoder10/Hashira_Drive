#include <iostream>
#include <fstream>
#include <vector>
#include <json/json.h> 
#include <cmath>
#include <cctype>

using namespace std;

long long baseToDecimal(const string& value, int base) {
    long long num = 0;
    for (char c : value) {
        int digit = 0;
        if (c >= '0' && c <= '9')
            digit = c - '0';
        else if (c >= 'a' && c <= 'z')
            digit = 10 + (c - 'a');
        else if (c >= 'A' && c <= 'Z')
            digit = 10 + (c - 'A');
        else {
            cerr << "Invalid digit " << c << " in input\n";
            exit(1);
        }
        if (digit >= base) {
            cerr << "Invalid Digit " << "\n";
            exit(1);
        }
        num = num * base + digit;
    }
    return num;
}


void gaussianElimination(vector<vector<double>>& A, vector<double>& Y, vector<double>& C) {
    int n = Y.size();
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int j = i + 1; j < n; j++)
            if (fabs(A[j][i]) > fabs(A[maxRow][i]))
                maxRow = j;
        swap(A[i], A[maxRow]);
        swap(Y[i], Y[maxRow]);

        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++)
                A[j][k] -= factor * A[i][k];
            Y[j] -= factor * Y[i];
        }
    }

    C.resize(n);
    for (int i = n - 1; i >= 0; i--) {
        C[i] = Y[i];
        for (int j = i + 1; j < n; j++)
            C[i] -= A[i][j] * C[j];
        C[i] /= A[i][i];
    }
}

bool isNumber(const string& s) {
    for (char c : s)
        if (!isdigit(c))
            return false;
    return !s.empty();
}

int main() {
    ifstream file("input.json");
    if (!file.is_open()) {
        cerr << "Failed to open input.json\n";
        return 1;
    }

    Json::Value root;
    file >> root;

    int n = root["keys"]["n"].asInt();
    int k = root["keys"]["k"].asInt();
    int m = k - 1;

    vector<double> Xvals;
    vector<double> Yvals;

    for (const auto& key : root.getMemberNames()) {
        if (key == "keys")
            continue;

        if (!isNumber(key))
            continue;  

        int x = stoi(key);
        int base = stoi(root[key]["base"].asString());
        string val = root[key]["value"].asString();

        long long yDecoded = baseToDecimal(val, base);
        Xvals.push_back((double)x);
        Yvals.push_back((double)yDecoded);
    }

    if (Xvals.size() < (size_t)k) {
        cerr << "Insufficient roots to solve polynomial\n";
        return 1;
    }

    vector<vector<double>> A(k, vector<double>(k));
    vector<double> Y(k);

    for (int i = 0; i < k; i++) {
        double x = Xvals[i];
        double y = Yvals[i];
        Y[i] = y;
        A[i][0] = 1.0;
        for (int j = 1; j < k; j++)
            A[i][j] = pow(x, j);
    }

    vector<double> coeffs;
    gaussianElimination(A, Y, coeffs);

    cout << "Coefficients c (polynomial degree " << m << "):" << endl;
    for (int i = 0; i < (int)coeffs.size(); i++)
        cout << "c[" << i << "] = " << coeffs[i] << endl;

    return 0;
}
