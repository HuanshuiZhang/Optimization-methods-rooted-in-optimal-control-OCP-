# Nonlinear Function Minimization in MATLAB

This MATLAB function provides two algorithms for finding the minimum of a nonlinear function. Both algorithms involve first- and second-order information of the objective function and are designed to converge super-linearly.

## Algorithms

The function implements two methods:

1. **Method 1 (Algorithm [1])**
   - **Characteristics**: 
     - Faster convergence
     - Higher precision
     - More computationally expensive due to the need to invert a matrix with the same dimension as the Hessian matrix
   - **Reference**: [Zhang et al. (2023)](https://arxiv.org/abs/2312.01334)

2. **Method 2 (Algorithm [2])**
   - **Characteristics**: 
     - Slower convergence compared to Method 1
     - Less precise
     - Computationally cheaper
   - **Reference**: [Wang et al. (2024)](https://arxiv.org/abs/2403.11115)

## Comparison with Other Methods

- **Gradient Descent**: The provided methods converge faster than gradient descent.
- **Newton's Method**: The provided methods are more stable compared to Newton's method.
- Both methods (Method 1 and Method 2) can handle cases where the Hessian matrix is singular or indefinite at some points.

## Usage

To use the MATLAB function, select the desired method based on your requirements for calculation cost and accuracy:

- To implement **Method 1**, select the appropriate option in the function.
- To implement **Method 2**, select the corresponding option in the function.

Choose the method that best suits your needs for convergence speed and computational efficiency.
