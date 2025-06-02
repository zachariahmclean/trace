
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace Rcpp;

double calculateRSquared(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    if (n != y.size() || n < 2) return 0.0;
    
    double sum_x = 0.0, sum_y = 0.0;
    double sum_xy = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
    
    for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
    }
    
    double numerator = n * sum_xy - sum_x * sum_y;
    double denominator = std::sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y));
    
    if (denominator == 0.0) return 0.0;
    double correlation = numerator / denominator;
    return correlation * correlation;
}

// Structure to hold result combinations with their R-squared values
struct CombResult {
    std::vector<std::vector<double>> combinations;
    std::vector<double> rsq_values;
};

// Find the best combinations based on R-squared
CombResult findBestCombinations(const std::vector<std::vector<double>>& recombinations, 
                               const std::vector<double>& reference_sizes, 
                               int top_n) {
    std::vector<double> rsq_vector(recombinations.size(), 0.0);
    
    for (size_t i = 0; i < recombinations.size(); i++) {
        rsq_vector[i] = calculateRSquared(reference_sizes, recombinations[i]);
    }
    
    // Find indices of top N combinations
    std::vector<size_t> indices(rsq_vector.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    std::partial_sort(indices.begin(),
                    indices.begin() + std::min(top_n, (int)indices.size()),
                    indices.end(),
                    [&rsq_vector](size_t i1, size_t i2) { return rsq_vector[i1] > rsq_vector[i2]; });
    
    CombResult result;
    for (int i = 0; i < std::min(top_n, (int)indices.size()); i++) {
        result.combinations.push_back(recombinations[indices[i]]);
        result.rsq_values.push_back(rsq_vector[indices[i]]);
    }
    
    return result;
}

// Generate combinations of size k from vector v
std::vector<std::vector<double>> generateCombinations(const std::vector<double>& v, int k) {
    std::vector<std::vector<double>> result;
    std::vector<bool> bitmask(k, 1); // k leading 1's
    bitmask.resize(v.size(), 0); // v.size()-k trailing 0's
    
    // While there are still valid bitmasks
    do {
        std::vector<double> combination;
        for (size_t i = 0; i < v.size(); ++i) {
            if (bitmask[i]) combination.push_back(v[i]);
        }
        result.push_back(combination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    
    return result;
}

// Recursive function to explore assignments with backtracking
CombResult exploreAssignments(std::vector<double> ref, 
                             std::vector<double> obs,
                             int max_combinations, 
                             int branching_n) {
    
    // Calculate window size for this iteration
    int n_obs = obs.size();
    int n_ref = ref.size();
    
    // Check if we'll generate too many combinations
    // Calculate binomial coefficient C(n_obs, n_ref)
    long long n_recombinations = 1;
    int k = std::min(n_ref, n_obs - n_ref); // Use smaller k for efficiency
    for (int i = 1; i <= k; i++) {
        n_recombinations = n_recombinations * (n_obs - i + 1) / i;
    }
    
    if (n_recombinations > max_combinations) {
        stop("Recombinations (" + std::to_string(n_recombinations) + ") exceeded max_combinations");
    }
    
    // Generate combinations
    std::vector<std::vector<double>> all_combinations = generateCombinations(obs, n_ref);
    
    // Get top combinations
    CombResult top_combinations = findBestCombinations(all_combinations, ref, branching_n);
    
    return top_combinations;
}

// Main function exposed to R
// [[Rcpp::export]]
List ladder_assignment_cpp(NumericVector reference_sizes, 
                          NumericVector observed_sizes,
                          int max_combinations = 10000000,
                          int branching_n = 5) {
    
    // Convert R vectors to C++ vectors
    std::vector<double> ref_vec = as<std::vector<double>>(reference_sizes);
    std::vector<double> obs_vec = as<std::vector<double>>(observed_sizes);
    
    CombResult result = exploreAssignments(ref_vec, obs_vec, max_combinations, branching_n);
    
    // Convert combinations to R format
    List combinations_list(result.combinations.size());
    for (size_t i = 0; i < result.combinations.size(); i++) {
        combinations_list[i] = NumericVector(result.combinations[i].begin(), result.combinations[i].end());
    }
    
    // Create a list to return to R
    List result_list = List::create(
        Named("recombinations") = combinations_list,
        Named("rsq_values") = NumericVector(result.rsq_values.begin(), result.rsq_values.end())
    );
    
    return result_list;
}