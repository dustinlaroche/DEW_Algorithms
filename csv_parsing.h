//
// Created by dlaroche on 6/2/24.
//

#ifndef CSV_PARSING_H
#define CSV_PARSING_H
#include <complex>
#include <fstream>
#include <vector>
#include <iostream>
#endif //CSV_PARSING_H

using namespace std;
typedef complex<double> Complex;


// Function to read CSV file and extract ADC sample values as real values
pair<vector<double>, vector<double>> retrieve_columns_as_real_voltages(const string& csv_file, int negative_column, int positive_column, int ADCMaxValue, double vref, bool fixInverted) {
    vector<double> negative_samples;
    vector<double> postive_samples;
    ifstream file(csv_file);

    if (!file.is_open()) {
        throw runtime_error("Failed to open CSV file");
    }

    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string str_value;
        int column = 0;

        // input voltage range = - (gain + v_offset)  to (v_ref - v_offset)
        // Adjusted Voltage = (adc_sample / adc_max) * (v_ref - v_offset)
        // Voltage = (adc_sample / (2^Resolution - 1)) * Vref
        while (getline(ss, str_value, ',')) {
            if (column == negative_column) {
                if(fixInverted == true) {
                    double gain = vref / ADCMaxValue;
                    double voltage_offset = gain * cos(M_PI);
                    double adjusted_voltage = (stod(str_value) / ADCMaxValue) * (vref - voltage_offset);
                    negative_samples.push_back(adjusted_voltage);
                } else {
                    double voltage = (stod(str_value) / ADCMaxValue) * vref;
                    negative_samples.push_back(voltage);
                }
            } else if (column == positive_column) {
                postive_samples.push_back(stod(str_value));
            }
            column++;
        }
    }

    file.close();

    return make_pair(negative_samples, postive_samples);
}
