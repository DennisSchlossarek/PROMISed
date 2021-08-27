

virtual machine engine

https://www.virtualbox.org/

then download and install ubuntu 20.04 desktop version

here to install docker on ubuntu

https://docs.docker.com/engine/install/ubuntu/


then

go to folder in with Dockerfile and your application

from here run a command to build a docker ( a dot `.` is not a mistake)

```
  sudo docker build -t promised .
```

after a build please restart a docker engine

```
  sudo systemctl restart docker.service
```

and then run a dockcer contaner in command line attached to terminal

```
  sudo docker run -p 30000:3838 promised
```

30000 is just a port number - it can be any port but better higher number then lower



