hepd_dir=$PWD
mkdir build
cd build
echo "Configuration  with cmake!"
cmake -DSTATIC_BUILD=OFF -DHEPDSW_INSTALL_PREFIX=$hepd_dir $hepd_dir
echo "Compilation of the code with make!"
make -j4
echo "Installation of the code!"
make install
