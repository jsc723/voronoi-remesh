### How to Build:

```
mkdir build
cd build
cmake ..
make
```

### How to Run:

```
./final_proj voronoi [options] <inputfile_name> 
```

options:
```
-n INT    : set the number of clusters, default 400
-f PATH   : set the name of the output file, a default name will be generated if not provided
-r INT    : set the radius of the ring used by adaptive sampling, default 2
-seed INT : set the seed for srand(), default time(0)
-q        : use quadric based placement policy for vertex position
-o1       : enable adaptive remeshing
-o2       : enable adaptive remeshing with anisotropic metric
-o3       : shortcut for -o2 -q
-v        : enable verification for the clusters
-h        : show this help page
```
### Some examples:
#### The original mesh
<img src="img/block.png" alt="" width="600" />

#### Uniform remeshing
<img src="img/block-iso.png" alt="" width="600" />

#### Adaptive remeshing
<img src="img/block-adapt.png" alt="" width="600" />

#### Adaptive remeshing with anisotrophic metric
<img src="img/block-aniso.png" alt="" width="600" />

#### Adaptive remeshing with anisotrophic metric + quadric based placement policy
<img src="img/block-vp.png" alt="" width="600" />

#### Uniform remeshing for the dragon model
<img src="img/dragon-iso.png" alt="" width="600" />

#### Adaptive remeshing with anisotrophic metric + quadric based placement policy for the dragon model
<img src="img/dragon-vp.png" alt="" width="600" />


