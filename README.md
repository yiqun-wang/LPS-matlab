# LPS-matlab

This code is a MATLAB implementation of **Local Point Signature** for 3D surface shape matching described in our CVPR 2019 paper:

["A Robust Local Spectral Descriptor for Matching Non-Rigid Shapes with Incompatible Shape Structures"](http://openaccess.thecvf.com/content_CVPR_2019/html/Wang_A_Robust_Local_Spectral_Descriptor_for_Matching_Non-Rigid_Shapes_With_CVPR_2019_paper.html) 

by Yiqun Wang, Jianwei Guo, Dong-Ming Yan, Kai Wang, Xiaopeng Zhang.

[Project Page](http://www.nlpr.ia.ac.cn/ivc/project/specmathcing/)

Please consider citing the above paper if you use the code/program (or part of it). 

## License

This program is free software; you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation; either version 2 of 
the License, or (at your option) any later version. 

## Usage	  

Please run "demo.m" using MATLAB for generating Local Point Signature(LPS).

1. input: "off", "obj", or "ply" format of mesh model in mesh_dir folder.

2. output: Local Point Signature on every vertex saved in ourput_dir folder. 
