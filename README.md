# Test project for SIGMOD2014 PC / Q1

## Data

```bash
wget https://www.dropbox.com/s/dktnm1h5d1syxkw/p100k.zip
unzip p100k.zip
```

## Build and test

```bash
rm -rf cmake-build-release && mkdir cmake-build-release && pushd cmake-build-release && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. && make -j$(nproc) && ./main ; popd
```
