# A collection of tensor completion algorithms
This is an official implementation of **S-LRTC, BS-TMac, SpBCD, EPT-BCD, MCP-BCD, and SCAD-BCD** via Matlab R2016. The implementation is based on [FaLRTC](https://ieeexplore.ieee.org/document/6138863), [TMac](https://link.springer.com/article/10.1007/s11464-012-0194-5) , and [geomCG](https://link.springer.com/article/10.1007/s10543-013-0455-z), we thanks the authors for sharing their code. One can quickly test these algorithms as follows.
### Quickly test "S-LRTC"
```S-LRTC``` is an implementation of "**A Mixture of Nuclear Norm and Matrix Factorization for Tensor Completion**", one can quickly test the algorithm by running 
```matlab
test_S_LRTC.m 
```
### Quickly test "BS-TMac"
```BS-TMac``` is an implementation of "**Robust balancing scheme-based approach for tensor completion**", one can quickly test the algorithm by running ```
```matlab
test_BS_TMac.m
```

### Quickly test "SpBCD"
```SpBCD``` is an implementation of "**Robust Schatten-p Norm Based Approach for Tensor Completion**", one can quickly test the algorithm by running
```matlab
test_SpBCD.m
```

### Quickly test "EPT-BCD, MCP-BCD, and SCAD-BCD"
```EPT-BCD```, ```MCP-BCD```, and ```SCAD-BCD``` is the implementation of "**Robust approximations of low-rank minimization for tensor completion**", one can quickly test the algorithms by running 
```matlab
test_EPT_MCP_SCAD.m
```

### Related works
[1] Liu, J., Musialski, P., Wonka, P., Ye, J.P.: Tensor completion for estimating missing values in visual data. IEEE Trans. Pattern Anal. Mach. Intell. 35(1), 1126–1153 (2013) [[link]](https://ieeexplore.ieee.org/document/6138863)

[2] Kressner, D., Steinlechner, M., Vandereycken, B.: Low-rank tensor completion by Riemannian optimization. BIT Numer. Math. 54(2), 447–468 (2014) [[link]](https://link.springer.com/article/10.1007/s10543-013-0455-z)

[3] Xu, Y., Hao, R., Yin, W., Su, Z.: Parallel matrix factorization for low-rank tensor completion. Inverse Probl. Imaging 9(2), 601–624 (2015) [[link]](https://link.springer.com/article/10.1007/s11464-012-0194-5)

### Citation
If our works are useful in your research or publication, please cite the works:

[1] S. Gao and Q. Fan. A Mixture of Nuclear Norm and Matrix Factorization for Tensor Completion. J Sci Comput, 75:43–64, 2018. [[link]](https://doi.org/10.1007/s10915-017-0521-9)

[2] S. Gao and Q. Fan. Robust balancing scheme-based approach for tensor completion. Neurocomputing, 330:328–336, 2019. [[link]](https://doi.org/10.1016/j.neucom.2018.11.033)

[3] S. Gao and Q. Fan. Robust Schatten-p Norm Based Approach for Tensor Completion. J Sci Comput, 82:11, 2020. [[link]](https://doi.org/10.1007/s10915-019-01108-9)

[4] S. Gao and X. Zhuang. Robust approximations of low-rank minimization for tensor completion. Neurocomputing, 379:319–333, 2020. [[link]](https://doi.org/10.1016/j.neucom.2019.10.086)

Don't hesitate to contact us via [shqgao@163.com](), if you have any questions.

