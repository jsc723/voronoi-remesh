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

