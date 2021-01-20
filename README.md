# Test project for SIGMOD2014 PC / Q1

## Data

```bash
wget http://mit.bme.hu/~szarnyas/spc2014/data/p100k.zip
unzip p100k.zip
```

## Build and test

```bash
rm -rf cmake-build-release && mkdir cmake-build-release && pushd cmake-build-release && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. && make -j$(nproc) && ./main ; popd
```
