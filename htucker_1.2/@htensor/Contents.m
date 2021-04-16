% htucker toolbox
% Version 1.2
%
% @htensor:
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

% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
