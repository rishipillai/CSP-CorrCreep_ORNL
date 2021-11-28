# Corrosion and Creep Simulation
- Written by Rishi, Matt Lee (lees4@ornl.gov)

# How to use

Please install Docker Desktop for Windows from the following link: https://hub.docker.com/search?q=&type=edition&offering=community&sort=updated_at&order=desc When installation is complete, open a terminal window, and go to Docker directory where Dockerfile is located:

Then run the following commands to build docker image
```
docker build . -t corrsim 
```

Then, please run a command similar to:
```
docker run -d -p8888:8888 -p8088:8088 -v [Absolute path to store your data]:/web_interface/output corrsim
```

For instance, if you want to save output files in /Dcouments/Result then, run
```
docker run -d -p8888:8888 -p8088:8088 -v /Documents/Result:/web_interface/output corrsim
```

Open your web browser and go to http://localhost:8088 to start the GUI for the simulation