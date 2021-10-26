# Start GRASS GIS docker and bind this folder, then enter the container

name="longrisk-gip"
srcbindpoint=/usr/src/app/
databindpoint=/data/
datafolder=$HOME/data/longrisk/
# Limit container's resource usage
cpus='12.5'
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
    --volume $HOME/data/storage/tmp/longrisk/:/tmp/gip/ \
    --env HOME=$srcbindpoint \
    --env LOGNAME=$LOGNAME \
    --env GISRC=$srcbindpoint/.grass7/rc \
    $(id -n -u)/grass-itzi \
    /bin/bash

# docker run -it --rm --user=$(id -u):$(id -g) \
#     --volume /your/test/grassdata/:/data --env HOME=/data/ grassgis80 \
#         grass /data/nc_basic_spm_grass7/PERMANENT --exec g.region -p

unset name srcbindpoint databindpoint datafolder cpus memory_swap memory
