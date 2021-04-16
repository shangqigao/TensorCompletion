% htucker toolbox
% Version 1.2
%
% Toolbox for tensors in hierarchical Tucker decomposition (HTD)
%
% Functions.
%   truncate       - Truncate full tensor/htensor/CP to htensor.
%
% Example tensors.
%   gen_invlaplace - htensor for approx. inverse of Laplace-like matrix.
%   gen_laplace    - htensor for Laplace-like matrix.
%   gen_sin_cos    - Function-valued htensor for sine and cosine.
%   htenones       - htensor with all elements one.
%   htenrandn      - Random htensor.
%   laplace_core   - Core tensor for Laplace operator.
%   reciproc_sum   - Function-valued tensor for 1/(xi_1+ ... +xi_d)
%   
% Auxiliary functions for full tensors.
%   dematricize    - Determine (full) tensor from matricization.
%   diag3d         - Return third-order diagonal tensor.
%   isindexvector  - Check whether input is index vector.
%   khatrirao_aux  - Khatri-Rao product.
%   khatrirao_t    - Transposed Khatri-Rao product.
%   matricize      - Return matricization of (full) tensor.
%   spy3           - Plot sparsity pattern of order-3 tensor.
%   ttm            - N-mode multiplication of (full) tensor with matrix
%   ttt            - Tensor times tensor (full tensors).
%
% Test
%   test_functions - Simple test script to check for obvious bugs.
%
%
% --------------------------------------------------------------------------
%                                   @htensor/
% --------------------------------------------------------------------------
%
% Construction of htensor objects.
%   htensor            - Construct a tensor in HTD and return htensor object.
%   define_tree        - Define dimension tree.
%
% Basic functionality.
%   cat                - Concatenate two htensor objects.
%   change_dimtree     - Change dimension tree of htensor.
%   change_root        - Change root of the dimension tree.
%   check_htensor      - Check internal consistency of htensor.
%   conj               - Complex conjugate of htensor.
%   ctranspose         - Not defined for htensor.
%   disp               - Command window display of htensor.
%   display            - Command window display of htensor.
%   disp_all           - Command window display of dimension tree of htensor.
%   end                - Last index in one mode of htensor.
%   equal_dimtree      - Compare dimension trees of two htensor objects.
%   full               - Convert htensor to a (full) tensor.
%   full_block         - Return subblock of htensor as a (full) tensor.
%   full_leaves        - Convert leaf matrices U to dense matrices.
%   ipermute           - Inverse permute dimensions of htensor.
%   isequal            - Check whether two htensors are equal.
%   mrdivide           - Scalar division for htensor.
%   mtimes             - Scalar multiplication for htensor.
%   ndims              - Order (number of dimensions) of htensor.
%   ndofs              - Number of degrees of freedom in htensor.
%   norm               - Norm of htensor.
%   norm_diff          - Norm of difference between htensor and full tensor.
%   nvecs              - Dominant left singular vectors for matricization of htensor.
%   permute            - Permute dimensions of htensor.
%   plot_sv            - Plot singular value tree of htensor.
%   rank               - Hierarchical ranks of htensor.
%   singular_values    - Singular values for matricizations of htensor.
%   size               - Size of htensor.
%   sparse_leaves      - Convert leaf matrices U to sparse matrices.
%   spy                - Plot sparsity pattern of the nodes of htensor.
%   squeeze            - Remove singleton dimensions from htensor.
%   subsasgn           - Subscripted assignment for htensor.
%   subsref            - Subscripted reference for htensor.
%   subtree            - Return all nodes in the subtree of a node.
%   transpose          - Not defined for htensor.
%   uminus             - Unary minus (-) of htensor.
%   uplus              - Unary plus for htensor.
%
% Operations with htensor objects.
%   elem_mult          - Approximate element-by-element multiplication for htensor.
%   innerprod          - Inner product for htensor.
%   minus              - Binary subtraction for htensor.
%   plus               - Binary addition for htensor.
%   power              - Element-by-element square for htensor.
%   times              - Element-by-element multiplication for htensor.
%   ttm                - N-mode multiplication of htensor with matrix.
%   ttt                - Tensor-times-tensor for htensor.
%   ttv                - Tensor-times-vector for htensor.
%
% Orthogonalization and truncation.
%   gramians           - Reduced Gramians of htensor in orthogonalized HTD.
%   gramians_cp        - Reduced Gramians of CP tensor.
%   gramians_nonorthog - Reduced Gramians of htensor.
%   gramians_sum       - Reduced Gramians for sum of htensor objects.
%   left_svd_gramian   - Left singular vectors and values from Gramian.
%   left_svd_qr        - Left singular vectors and values.
%   orthog             - Orthogonalize HTD of htensor.
%   trunc_rank         - Return rank according to user-specified parameters.
%   truncate_cp        - Truncate CP tensor to lower-rank htensor.
%   truncate_ltr       - Truncate full tensor to htensor, leaves-to-root.
%   truncate_nonorthog - Truncate htensor to lower-rank htensor.
%   truncate_rtl       - Truncate full tensor to htensor, root-to-leaves.
%   truncate_std       - Truncate htensor to lower-rank htensor.
%   truncate_sum       - Truncate sum of htensor objects to lower-rank htensor.
%
% Linear Operators.
%   apply_mat_to_mat   - Applies an operator in HTD to another operator in HTD.
%   apply_mat_to_vec   - Applies an operator in HTD to htensor.
%   full_mat           - Full matrix represented by an operator in HTD.
%   innerprod_mat      - Weighted inner product for htensor.
% 
% Interface with Tensor Toolbox
%   ktensor_approx     - Approximation of htensor by ktensor.
%   mttkrp             - Building block for approximating htensor by ktensor.
%   ttensor            - Convert htensor into a Tensor Toolbox ttensor.
%
%
% --------------------------------------------------------------------------
%                                   examples/
% --------------------------------------------------------------------------
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
%   cookies_matrices_2x2  - FE discretization of parametric PDE
%                           with 4 parameters
%   cookies_matrices_3x3  - FE discretization of parametric PDE
%                           with 9 parameters


% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
