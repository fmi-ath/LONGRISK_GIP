echo "Building container"
docker build -t $(id -n -u)/grass-itzi --rm -f ./Containerfile .
