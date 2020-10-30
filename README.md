# Test project for SIGMOD2014 PC / Q3

## Requirements

Use the `master` branches of:
- SuiteSparse:GraphBLAS (v4.0.0draft5)
- LAGraph

## Build and test

```bash
rm -rf cmake-build-release && mkdir cmake-build-release && pushd cmake-build-release && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. && make -j$(nproc) && ./main ; popd
```
