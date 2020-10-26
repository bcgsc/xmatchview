FROM  ubuntu:20.04 as builder

USER  root

ARG VER_MINIMAP2="2.17"

ENV OPT=/opt/xmatchview
# have to split as OPT needs to exist to use in the next ENV
ENV PATH=$OPT/bin:$PATH

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends curl ca-certificates

WORKDIR /install_tmp

RUN mkdir -p $OPT/bin

RUN curl --retry 10 -o distro.tar.bz2 -sSL https://github.com/lh3/minimap2/releases/download/v${VER_MINIMAP2}/minimap2-${VER_MINIMAP2}_x64-linux.tar.bz2 \
&& mkdir -p distro \
&& tar -C distro --strip-components 1 -xjf distro.tar.bz2 \
&& cp distro/minimap2 $OPT/bin/.

COPY tarballs/fonts $OPT/fonts
RUN chmod +r -R $OPT/fonts

# add tests to OPT bin
RUN mkdir -p $OPT/tests
COPY test $OPT/tests/test
COPY test-hive $OPT/tests/test-hive
RUN chmod +r -R $OPT/tests
RUN chmod +x $OPT/tests/test/*.sh

# add the scripts to the OPT bin
COPY *.py $OPT/bin/
COPY *.sh $OPT/bin/
RUN chmod +rx $OPT/bin/*.sh $OPT/bin/*.py

FROM  ubuntu:20.04

 LABEL maintainer="Rene Warren (rwarren@bcgsc.ca)"\
       description="xmatchview"

ENV OPT=/opt/xmatchview
# have to split as OPT needs to exist to use in the next ENV
ENV PATH=$OPT/bin:/opt/cross_match/bin:$PATH \
    XM_FONTS=$OPT/fonts

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
python3 python3-pip \
unattended-upgrades && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN pip3 install pillow

RUN mkdir -p $OPT
COPY --from=builder $OPT $OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
