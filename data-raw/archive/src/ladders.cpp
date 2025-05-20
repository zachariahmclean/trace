# function moved here since although it worked, it adds significant complexity to fix bugs
# revisit if faster ladder is required




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
    
    // something here to make denominator functionally zero. do else so not double return

    if (denominator == 0.0) return 0.0;
    double correlation = numerator / denominator;
    return correlation * correlation;
}

// Calculate mean R-squared over sliding windows
double meanRSq(const std::vector<double>& scan, const std::vector<double>& size) {
    if (scan.size() < 3 || size.size() < 3) return 0.0;
    
    std::vector<double> cors;
    for (size_t i = 0; i <= scan.size() - 3; i++) {
        std::vector<double> xi = {scan[i], scan[i+1], scan[i+2]};
        std::vector<double> yi = {size[i], size[i+1], size[i+2]};
        cors.push_back(calculateRSquared(yi, xi));
    }
    
    double sum = std::accumulate(cors.begin(), cors.end(), 0.0);
    return sum / cors.size();
}

// Structure to hold result combinations with their R-squared values
struct CombResult {
    std::vector<std::vector<double>> combinations;
    std::vector<double> rsq_values;
};

// Find the best combinations based on R-squared
CombResult findBestCombinations(const std::vector<std::vector<double>>& recombinations, 
                               const std::vector<double>& reference_sizes, int top_n = 3) {
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

// Structure to hold the assignment result
struct AssignmentResult {
    std::vector<double> assigned_observed;
    std::vector<double> assigned_reference;
    double final_rsq;
    bool valid;
    
    AssignmentResult() : final_rsq(-1.0), valid(false) {}
};

// Recursive function to explore assignments with backtracking
AssignmentResult exploreAssignments(std::vector<double> remaining_ref, 
                                  std::vector<double> remaining_obs,
                                  std::vector<double> current_assigned_ref, 
                                  std::vector<double> current_assigned_obs,
                                  double current_rsq, int depth = 1, 
                                  bool comb_reversed = false,
                                  int choose, 
                                  int max_combinations, 
                                  double r_squared_threshold) {
    
    AssignmentResult result;
    
    // If remaining reference and observed are same size, assign them all at once
    if (remaining_ref.size() == remaining_obs.size()) {
        current_assigned_obs.insert(current_assigned_obs.end(), remaining_obs.begin(), remaining_obs.end());
        current_assigned_ref.insert(current_assigned_ref.end(), remaining_ref.begin(), remaining_ref.end());
        remaining_ref.clear();
        remaining_obs.clear();
    }
    
    // Base case: all reference sizes assigned
    if (remaining_ref.empty()) {
        if (comb_reversed) {
            result.assigned_observed = current_assigned_ref;
            result.assigned_reference = current_assigned_obs;
        } else {
            result.assigned_observed = current_assigned_obs;
            result.assigned_reference = current_assigned_ref;
        }
        result.final_rsq = current_rsq;
        result.valid = true;
        return result;
    }
    
    // Calculate window size for this iteration
    int n_obs = remaining_obs.size();
    int n_ref = remaining_ref.size();
    int current_choose = std::min(choose, n_ref);
    
    // Calculate start window
    int remainder = n_ref - current_choose;
    int start_window = std::min(n_obs - remainder, n_obs);
    start_window = std::max(start_window, current_choose); // Ensure we have enough points
    
    // Check if we'll generate too many combinations
    // Calculate binomial coefficient
    long long n_recombinations = 1;
    for (int i = 1; i <= current_choose; i++) {
        n_recombinations = n_recombinations * (start_window - i + 1) / i;
    }
    
    if (n_recombinations > max_combinations) {
        // Fall back to simpler approach for large combinations
        start_window = std::min(start_window, 2 * current_choose);
        n_recombinations = 1;
        for (int i = 1; i <= current_choose; i++) {
            n_recombinations = n_recombinations * (start_window - i + 1) / i;
        }
        
        if (n_recombinations > max_combinations) {
            result.valid = false;
            return result;
        }
    }
    
    // Generate subset of observations for combinations
    std::vector<double> obs_subset(remaining_obs.begin(), remaining_obs.begin() + start_window);
    
    // Create reference subset
    std::vector<double> ref_subset(remaining_ref.begin(), remaining_ref.begin() + current_choose);
    
    // Generate combinations
    std::vector<std::vector<double>> all_combinations = generateCombinations(obs_subset, current_choose);
    
    // Get top combinations
    CombResult top_combinations = findBestCombinations(all_combinations, ref_subset, 3);
    
    AssignmentResult best_result;
    double best_rsq = -1.0;
    
    // Try each top combination
    for (size_t i = 0; i < top_combinations.combinations.size(); i++) {
        std::vector<double>& selected_comb = top_combinations.combinations[i];
        double current_window_rsq = top_combinations.rsq_values[i];
        
        // If R² is too low, skip this branch
        if (current_window_rsq < r_squared_threshold && depth > 1) {
            continue;
        }
        
        // Find positions in remaining observations
        std::vector<size_t> positions;
        for (double val : selected_comb) {
            auto it = std::find(remaining_obs.begin(), remaining_obs.end(), val);
            if (it != remaining_obs.end()) {
                positions.push_back(std::distance(remaining_obs.begin(), it));
            }
        }
        
        // Set up for next iteration
        std::vector<double> new_assigned_obs = current_assigned_obs;
        new_assigned_obs.insert(new_assigned_obs.end(), selected_comb.begin(), selected_comb.end());
        
        std::vector<double> new_assigned_ref = current_assigned_ref;
        new_assigned_ref.insert(new_assigned_ref.end(), ref_subset.begin(), ref_subset.end());
        
        // Calculate overall R² so far
        double overall_rsq = 1.0;
        if (new_assigned_obs.size() > 2) {
            overall_rsq = meanRSq(new_assigned_ref, new_assigned_obs);
        }
        
        // Only proceed if correlation is decent
        if (overall_rsq >= current_rsq * 0.8 || depth <= 2) {
            // Create new remaining vectors for next recursion
            std::vector<double> new_remaining_obs = remaining_obs;
            std::sort(positions.begin(), positions.end(), std::greater<size_t>());
            for (size_t pos : positions) {
                new_remaining_obs.erase(new_remaining_obs.begin() + pos);
            }
            
            std::vector<double> new_remaining_ref(remaining_ref.begin() + current_choose, remaining_ref.end());
            
            AssignmentResult recursive_result = exploreAssignments(
                new_remaining_ref, new_remaining_obs,
                new_assigned_ref, new_assigned_obs,
                overall_rsq, depth + 1, comb_reversed,
                choose, max_combinations, r_squared_threshold
            );
            
            if (recursive_result.valid && recursive_result.final_rsq > best_rsq) {
                best_result = recursive_result;
                best_rsq = recursive_result.final_rsq;
            }
        }
    }
    
    return best_result;
}

// Main function exposed to R
// [[Rcpp::export]]
DataFrame ladder_assignment_backtracking_cpp(NumericVector reference_sizes, 
                                           NumericVector observed_sizes,
                                           int choose = 5,
                                           int max_combinations = 10000,
                                           double r_squared_threshold = 0.9) {
    
    // Convert R vectors to C++ vectors
    std::vector<double> ref_vec = as<std::vector<double>>(reference_sizes);
    std::vector<double> obs_vec = as<std::vector<double>>(observed_sizes);
    
    AssignmentResult result;
    
    if (obs_vec.size() > ref_vec.size()) {
        result = exploreAssignments(
            ref_vec, obs_vec,
            std::vector<double>(), std::vector<double>(),
            1.0, 1, false,
            choose, max_combinations, r_squared_threshold
        );
    } else {
        result = exploreAssignments(
            obs_vec, ref_vec,
            std::vector<double>(), std::vector<double>(),
            1.0, 1, true,
            choose, max_combinations, r_squared_threshold
        );
    }
    
    if (!result.valid) {
        stop("Could not find a valid assignment");
    }
    
    // Create a data frame to return
    DataFrame df = DataFrame::create(
        Named("scan") = wrap(result.assigned_observed),
        Named("size") = wrap(result.assigned_reference)
    );
    
    return df;
}