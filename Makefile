# Build Docker image
.PHONY: docker_build
docker_build:
	docker build -t rapidqt .

# Run Docker container with robot-related environment (ROS, etc.)
.PHONY: docker_run_sim
docker_run:
	docker run -it --network host   --env="DISPLAY" --env="QT_X11_NO_MITSHM=1" --device=/dev/video0:/dev/video0 --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" --entrypoint /bin/bash rapidqt
