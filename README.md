## Group Regularization Method Simulations
Provides functions to run simulations using lasso, ridge, group lasso, 
thresholded group lasso, elastic net, sparse group lasso, and group bridge methods.

- The function `testMethods` Generates design matrix with user-specified grouping structure and runs user-specified group regularization methods.
- Returns l2 difference between generated true-beta values and predicted beta-values based on what method was specified. User can also calculate the l2 diff using the `L2diff` function. Sign matching is used to take pairs of elements between two vectors and sees whether or not both elements are in the same category (abs value less than 0.2, abs value greater than 0.2) using the `SignMatch` function.
