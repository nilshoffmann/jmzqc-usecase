FROM jupyter/minimal-notebook

LABEL maintainer="Nils Hoffmann <n.hoffmann@fz-juelich.de>"

USER root

# Install dependencies
RUN apt-get update && apt-get install -y \
  software-properties-common \
  curl

# Install Zulu OpenJdk 21 (LTS)
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 0xB1998361219BD9C9
# download and install the package that adds 
# the Azul APT repository to the list of sources 
RUN curl -O https://cdn.azul.com/zulu/bin/zulu-repo_1.0.0-3_all.deb
# install the package
RUN apt install ./zulu-repo_1.0.0-3_all.deb
# update the package sources
RUN apt update
RUN apt install -y zulu21-jdk
RUN rm zulu-repo_1.0.0-3_all.deb

USER $NB_USER

# Download the kernel release
RUN curl -Ls https://sh.jbang.dev | bash -s - app setup --force --fresh --verbose

RUN $HOME/.jbang/bin/jbang version
RUN $HOME/.jbang/bin/jbang trust add https://github.com/jupyter-java/
RUN $HOME/.jbang/bin/jbang install-kernel@jupyter-java rapaio

# Cleanup
# RUN rm ijava-kernel.zip

# Add README.md and sample notebooks to the home directory.
ADD "README.md" $HOME
ADD "*.ipynb" $HOME/

WORKDIR $HOME

RUN jupyter trust *.ipynb

# Set user back to non-privileged user.
USER $NB_USER
