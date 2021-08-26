# Start GRASS GIS docker and bind this folder, then enter the container

name="longrisk-gip"
bindpoint=/data/

docker run -it --rm \
    --user=$(id -u):$(id -g) \
    --name=$name \
    --volume $(pwd):$bindpoint \
    --env HOME=/home/anon/ \
    $(id -n -u)/grass-itzi \
    /bin/bash

# docker run -it --rm --user=$(id -u):$(id -g) \
#     --volume /your/test/grassdata/:/data --env HOME=/data/ grassgis80 \
#         grass /data/nc_basic_spm_grass7/PERMANENT --exec g.region -p

unset name bindpoint
