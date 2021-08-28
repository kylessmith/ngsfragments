#ifndef RUNNING_MEAN_H
#define RUNNING_MEAN_H

double calculate_sum(const double values[], int half_window);
double rolling_mean(double values[], double sum, int i, int window_size);
void rolling_mean_array(const double values[], double means[], int length, int window_size);

#endif