FROM centos:7

ENV SAMTOOLS_VER=1.7
ENV PREFIX_INSTALLATION=/opt/craig
ENV CRAIG_HOME="${PREFIX_INSTALLATION}"
ENV SAMTOOLS_HOME=/opt/samtools
ENV REGTOOLS_HOME=/opt/regtools
ENV PATH="${CRAIG_HOME}/bin:${CRAIG_HOME}/perl/bin:${CRAIG_HOME}/python/bin:${REGTOOLS_HOME}/build:${SAMTOOLS_HOME}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${CRAIG_HOME}/lib"


RUN yum install -y --disableplugin=fastestmirror epel-release && \
  yum clean all && \
  yum update -y --disableplugin=fastestmirror && \
  yum install -y --disableplugin=fastestmirror \
    autoconf \
    automake \
    boost-regex \
    bzip2 \
    bzip2-devel \
    cmake \
    gcc \
    gcc-c++ \
    git \
    libtool \
    make \
    ncurses-devel \
    python2-pip \
    xz-devel \
    zlib-devel && \
  yum clean all && rm -rf /var/cache/yum

WORKDIR /tmp

RUN git clone https://github.com/sparsehash/sparsehash.git
RUN cd sparsehash && \
  ./configure && \
  make install

RUN git clone https://github.com/griffithlab/regtools "${REGTOOLS_HOME}"
RUN cd "${REGTOOLS_HOME}" && \
  mkdir build && \
  cd build/ && \
  cmake .. && \
  make

RUN git clone https://github.com/axl-bernal/CraiG.git
RUN cd CraiG  && pwd && ls -l && \
  ./autogen.sh  && \
  ./configure --prefix="${PREFIX_INSTALLATION}" CXXFLAGS="${CXXFLAGS} -std=c++11" --enable-opt=yes --enable-mpi=no && \
  make && make install && make installcheck && \
  if [[ -f python/requirements.txt ]]; then pip install -r python/requirements.txt; fi

RUN curl -LO "https://gigenet.dl.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2" && \
  rm -rf "samtools-${SAMTOOLS_VER}" && \
  tar xf "samtools-${SAMTOOLS_VER}.tar.bz2" && \
  pushd "samtools-${SAMTOOLS_VER}" && \
  ./configure --prefix=/opt/samtools && \
  make && \
  make install && \
  popd

RUN pip install numpy
