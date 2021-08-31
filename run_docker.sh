# Start GRASS GIS docker and bind this folder, then enter the container

name="longrisk-gip"
srcbindpoint=/usr/src/app/
databindpoint=/data/
datafolder=${PWD%/*}/data

docker run -it --rm \
    --user=$(id -u):$(id -g) \
    --name=$name \
    --volume ${PWD}:$srcbindpoint \
    --volume $datafolder:$databindpoint:ro \
    --env HOME=$srcbindpoint \
    --env LOGNAME=$LOGNAME \
    --env GISRC=$srcbindpoint/.grass7/rc \
    $(id -n -u)/grass-itzi \
    /bin/bash

# docker run -it --rm --user=$(id -u):$(id -g) \
#     --volume /your/test/grassdata/:/data --env HOME=/data/ grassgis80 \
#         grass /data/nc_basic_spm_grass7/PERMANENT --exec g.region -p

unset name srcbindpoint databindpoint datafolder
