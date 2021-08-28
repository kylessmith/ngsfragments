#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


double calculate_sum(const double values[], int half_window)
{
    int window_size = half_window * 2;
    double sum = 0;
    int i;
    for (i = half_window; i < half_window+window_size; i++)
    {
        sum += values[i];
    }

    return sum;
}


double rolling_mean(double values[], double sum, int i, int window_size)
{
    // Increment i
    i++;

    // Re-calculate sum
    sum = sum - values[i - window_size];
    sum = sum + values[i];
    double mean = sum / (double)window_size;

    return mean;
}


void rolling_mean_array(const double values[], double means[], int length, int window_size)
{

    // Initialize
    int half_window = window_size / 2;
    double sum = calculate_sum(values, half_window);

    // Iterate over values
    int i;
    for (i = half_window; i < length - half_window; i++)
    {
        
        means[i] = sum / (double)window_size;
        sum -= values[i - half_window];
        sum += values[i + half_window];
    }

    return;
}