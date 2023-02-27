# jmzqc-usecase


## Building the Docker Container

 ./mvnw -B jib:dockerBuild --file pom.xml

## Runing the Docker Container

  docker run -it --volume=$(pwd):/tmp --rm jmzqc-usecase -f /tmp/20181113_010_autoQC01.mzML -o /tmp/proteomics-usecase.mzQC
