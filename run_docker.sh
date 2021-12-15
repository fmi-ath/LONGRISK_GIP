# Start GRASS GIS docker and bind this folder, then enter the container

name="longrisk-gip"
srcbindpoint=/usr/src/app/
databindpoint=/data/
datafolder="$(realpath ../data/)"
outputfolder="$(realpath ../out/)"
# Limit container's resource usage
cpus='20'
memory='32g'
memory_swap='33g' # True swap amount is memory_swap - memory.

docker run -it --rm \
    --user=$(id -u):$(id -g) \
    --name=$name \
    --cpus=$cpus \
    --memory=$memory \
    --memory-swap=$memory_swap \
    --volume ${PWD}:$srcbindpoint \
    --volume $datafolder:$databindpoint:ro \
    --volume $outputfolder:/tmp/gip/ \
    --env HOME=$srcbindpoint \
    --env LOGNAME=$LOGNAME \
    --env GISRC=$srcbindpoint/.grass7/rc \
    $(id -n -u)/grass-itzi \
    /bin/bash

unset name srcbindpoint databindpoint datafolder outputfolder cpus memory_swap memory
