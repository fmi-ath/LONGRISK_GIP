---
title: Using with Docker
author: Petteri Karsisto
year: 2021
---

# Running LONGRISK_GIP with Docker

This document includes information on how to run LONGRISK_GIP with Docker.

## Setting up the model parameters

You'll need to modify the model configuration to reflect the file system *inside* the container.
Meaning that to paths should reflect the volume binds (mounts) instead of the host system paths.
The volume bind paths are defined in run_docker.sh (see [Run](#Run) below)

## Building the container image

Here are quick instructions and explanations of the Docker-related files. See Docker's documentation
for more details about the Docker commands and their options.

### build_docker.sh

This script creates a container image based on the commands given in `Containerfile`. You'll need
to run this script when you have to update the image we're running. Basically when you've updated
the `Containerfile` or whatever files you are mentioning in the commands (e.g. `requirements.txt`)

The `-t` option names resulting image as `your-user-name/grass-itzi` - which can be used later on.
~~~bash
echo "Building container"
docker build -t $(id -n -u)/grass-itzi --rm -f ./Containerfile .
~~~

### Containerfile

Container files contain instructions on how the Docker engine should build the container image.
These files are known also as `Dockerfile`s, which is also their default name. However, there are
other container engines alongside Docker, so `Containerfile` would be more generic name.

Anyways, this short file simply states that "use the `mundialis/grass-py-pdal:stable-ubuntu` image
as the base image, then install specified Python libraries with `pip` on top of it".
~~~dockerfile
FROM mundialis/grass-py3-pdal:stable-ubuntu

WORKDIR /usr/src/app

COPY requirements.txt ./

RUN python3 -m pip install -r requirements.txt
~~~

### .dockerignore

Use .dockerignore file to ignore certain files or folders during the image building. It's similar
to .gitignore files, if you're familiar with them. In our case, lacking this file isn't a very big
deal as we are including only one new file in the container image - the `requirements.txt` - but
it's still a good practice to have this file regardless (it also helps in reducing the data amount
sent to the docker server during the `docker build` command).
~~~ignore
# Ignore .git
.git

# Ignore secrets and other private information
.gitignore
.dockerignore
.grass7

# Ignore files that are irrelevant for running the container
*.py[cdo]
*_docker.sh

# Ignore volatile files and runtime-generated files
GRASS_itzi
.*_history
~~~

## Run

The file `run_docker.sh` has the configuration and command to run the container. The configuration
parameters are explained below.

Once you run this script, you should enter the container automatically. The bash command line will
complain about missing group name and user name (`I have no name!`), but these can be safely
ignored. If this interactivity is not desired, then replace the `-it` options with `-d` (detached).
~~~bash
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

unset name srcbindpoint databindpoint datafolder cpus memory_swap memory
~~~

### Why `--user`?

We want to write the output files to somewhere we can access them later on without having to rely
on the container. Easiest way of doing this is to write the files back to the host system via
mounts. However, default user in Docker containers is `root`, which causes us issues. While we can
(typically) read the files written by `root`, we cannot move or delete them (without `sudo` rights).
To prevent this, we should execute the programs inside the container as a non-root user. Setting the
"inside user"'s user id and group id the same as ours we can then use the output files like they
were owned by us (which is technically true).

For storing the data that is intended to be used within or shared between containers, please look
up on how to use docker volumes.

### Data paths

- `srcbindpoint`: Mount the current folder to this path inside the container. This path is also set
  as the home folder.
- `databindpoint`: Mount the input data folder to this path inside the container. The folder is
  mounted as read-only.
- `datafolder`: Location of the input data in the host system
- `outputfolder`: A path to a folder where the LONGRISK_GIP generated files should be stored in the
  host system.

### Environment variables

We set certain additional environment variables in the `docker run` command. These variables are
passed to the container and are available for the processes within the container.
- `HOME`: We want to have a home folder when we enter the container with a non-existing user
- `LOGNAME`: Python function [getpass.getuser](https://docs.python.org/3/library/getpass.html#getpass.getuser)
  checks for this variable to deduce current user name. As we enter with a non-existing user,
  we don't have a username entry in `passwd` file so we need to provide it when entering.
- `GISRC`: Required for GRASS to know where its configuration file resides.

### Limiting Docker / container resource usage

Limit the container's resource usage by setting `--cpus`, `--memory`, and `--memory-swap`.
See Docker's documentation for detailed information.

**Too long; didn't read:**
- `--cpus=10` allows container to use up to as many threads as 10 CPUs have. Values such as `1.5` are allowed to fine-tune the thread amount.
- `--memory='32g'` allows the container to use up to 32GB of RAM.
- `--memory-swap='33g'` allows the container to use up to 33GB of *RAM and SWAP combined*. This setup gives 32GB of RAM and (33GB - 32GB =) 1GB of SWAP. If `--memory-swap == --memory`, then SWAP is disabled.

## Troubleshooting

### GISRC is not set / does not exist

1. If GISRC is empty, set it to point at the rc file.
    - Default is `$HOME/.grass7/rc`
2. If such file does not exist, then create one. It looks something like this:
    ~~~
    GISDBASE: /usr/src/app/GRASS_itzi/grassdata
    LOCATION_NAME: FINLAND
    MAPSET: Helsinki
    ~~~
   You can also start up a grass session and create the file that way, e.g.
   `grass /usr/src/app/GRASS_itzi/grassdata/FINLAND/Helsinki`
