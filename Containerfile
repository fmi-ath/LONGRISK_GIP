FROM mundialis/grass-py3-pdal:stable-ubuntu

# TODO: Create and switch to non-root user

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN python3 -m pip install -r requirements.txt
