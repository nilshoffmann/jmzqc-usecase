# jmzqc-usecase


## Building the Docker Container

  ./mvnw jib:dockerBuild

## Runing the Docker Container

  docker run -it --volume=$(pwd):/tmp --rm jmzqc-usecase -f /tmp/20181113_010_autoQC01.mzML -o /tmp/proteomics-usecase.mzQC
