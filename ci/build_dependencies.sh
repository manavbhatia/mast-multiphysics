#!/usr/bin/env bash

if [ "${TRAVIS_OS_NAME}" = linux ]; then # Ubuntu Linux

  if [ "${TRAVIS_DIST}" = xenial ]; then # Ubuntu 16.04 Xenial Xerus
    cd ${HOME} || exit

    # Python 3.6 apt repository (since its not neatly included in Ubuntu 16.04)
    sudo add-apt-repository -y ppa:deadsnakes/ppa

    # Regular libMesh/MAST dependencies.
    sudo apt-get -qq update
    sudo apt-get -qq install -y gfortran wget m4
    sudo apt-get -qq install -y openmpi-bin libopenmpi-dev
    sudo apt-get -qq install -y libpetsc3.6 libpetsc3.6.2-dev
    sudo apt-get -qq install -y libslepc3.6 libslepc3.6.1-dev libparpack2-dev
    sudo apt-get -qq install -y libboost-all-dev
    sudo apt-get -qq install -y libeigen3-dev
    sudo apt-get -qq install -y doxygen graphviz rsync
    sudo apt-get -qq install -y texlive-latex-base dvi2ps ghostscript
    sudo apt-get -qq install -y python3.6 python3.6-dev libpython3.6

    wget https://bootstrap.pypa.io/get-pip.py || exit
    sudo python3.6 get-pip.py || exit

    sudo which python
    sudo which python3
    sudo which python3.6
    python --version
    python3 --version
    python3.6 --version
    sudo which pip
    sudo which pip3
    #sudo pip3 install numpy

    # Update to later CMake release.
    wget https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Linux-x86_64.sh || exit
    sudo mkdir /opt/cmake || exit
    sudo sh cmake-3.15.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license || exit

  # elif [ "${TRAVIS_DIST}" = bionic ]; then # Ubuntu 18.04 Bionic Beaver
  #  sudo apt-get -qq update
  #  sudo apt-get -qq install -y gfortran wget m4
  #  which gfortran
  #  gcc --version
  #  echo "Hello From BIONIC"

  else
    echo "INVALID LINUX DISTRO: ${TRAVIS_DIST}"
    exit 1
  fi

elif [ "${TRAVIS_OS_NAME}" = osx ]; then # macOS 10.14, XCode 10.2
  # Currently we don't do anything here since we get all dependencies for macOS
  # from the binary download with "ci/get_libmesh.sh" in the next stage.
  echo "Hello From OSX"

else
  echo "INVALID OS: ${TRAVIS_OS_NAME}"
  exit 1
fi