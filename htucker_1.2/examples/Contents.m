% htucker toolbox
% Version 1.2
%
% Demos:
%   demo_basics           - Demonstration of basic htensor functionality.
%   demo_constructor      - Demonstration of htensor constructors.
%   demo_elem_reciprocal  - Demonstration of element-wise reciprocal.
%   demo_function         - Demonstration of htensor function approximation.
%   demo_invlaplace       - Demonstration of approximate inverse Laplace.
%   demo_operator         - Demonstration of operator-HTD format.
%
% Examples:
%   example_cancellation  - Cancellation in tan(x) + 1/x - tan(x).
%   example_cancellation2 - Cancellation in exp(-x^2) + sin(x)^2 + cos(x)^2.
%   example_cookies       - Apply CG method to a parametric PDE.
%   example_maxest        - Example for computing element of
%                           maximal absolute value.
%   example_spins         - Demonstration of operator-HTD for 1D spin system.
%   example_sum_reciproc  - Apply element-wise reciprocal method to sum of tensors.
%   example_truncation    - Comparison of speed for different
%                           truncation methods.
%
% Functions:
%   cg_tensor             - Truncated Conjugate Gradient method for htensor.
%   elem_reciprocal       - Iterative computation of elementwise
%                           reciprocal for htensor.
%   handle_inv_mat        - Function handle to Kronecker structured matrix multiplication.
%   handle_lin_mat        - Function handle to Kronecker structured matrix multiplication.
%   maxest                - Approximate element of maximal absolute value.
%
% MAT-files:
%   cookies_matrices_2x2       - FE discretization of parametric PDE
%                                with 4 parameters
%   cookies_matrices_3x3       - FE discretization of parametric PDE
%                                with 9 parameters
%   example_truncation_output  - Computation times of example_truncation


% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
