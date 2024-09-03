# Non-convex Function Minimization in MATLAB

This MATLAB function **OMBOC** (Optimization Method Based on Optimal Control), implements two methods for finding the unconstrained minimum of a non-convex function. These methods are based on optimal control theory and are designed to efficiently converge to a solution even in complex optimization landscapes.

## Methods

`OMBOC` is unique in its integration of optimization with control theory, providing a robust approach to solving challenging optimization problems.

`OMBOC ` provides two methods for optimization:

1. **Method 1 (Algorithm in [1])**
   - **Characteristics**:
     - Relatively simple parameter selection.
     - Find different local minimum points.
     - Higher computational time due to the need for solving FBDEs.
   - **Reference**: [Xu et al. (2023)](https://arxiv.org/abs/2309.05280)
2. **Method 2 (Algorithm in [2])**
   - **Characteristics**:
     - Find different local minimum points.
     - Higher precision.
     - Reduced computational time compared to Method 1.
   - **Reference**: [Zhang et al. (2023)](https://arxiv.org/abs/2312.01334)

## Comparison with Other Methods

- **Gradient Descent**:  `OMBOC`  offers faster convergence than traditional gradient descent methods.
- **Newton's Method**:  `OMBOC` is more stable and versatile, especially in cases where the Hessian matrix is singular or indefinite.

## Usage

- Before using the MATLAB function, you should download the **CasADi** toolkit and add it to your MATLAB path. 
- The different local minimum points can be found by adjusting **R** in some cases.
- Choose the appropriate control weight matrix **R**  and the step size **Lambda** to best suit your needs for convergence speed and computational efficiency.
- To use the MATLAB function, select the desired method based on your requirements for calculation cost and efficiency:

  - To implement **Method 1**, select the appropriate option in the function.
  - To implement **Method 2**, select the corresponding option in the function.


### Example Usage

```matlab
% Define the objective function
myfun = @(x) x - 4*x^2 + 0.2*x^3 + 2*x^4;

% Initial value
x0 = 0.5;

% Call OMBOC with Method 1
[x, fval] = OMBOC(@(x)myfun(x), x0, eye(1), 'Method1', 100, 0.1, [], 1e-3, 'on');

% Display the result
disp(['Optimal X: ', num2str(x)]);
disp(['Function Value: ', num2str(fval)]);
```