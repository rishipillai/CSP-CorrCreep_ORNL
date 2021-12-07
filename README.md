# Corrosion and Creep Simulation
- Written by Rishi Pillai (pillairr@ornl.gov), Marie Romedenne (romedennem@ornl.gov) and  Matt Lee (lees4@ornl.gov)


# How to use in console mode if all required modules and python are pre-installed

- python creep.py data/*.txt output/$.csv
- python diffusion.py data/*.txt output/$.csv
- where * is CreepProp_282.txt for creep.py and CorrProp_282.txtfor diffusion.py and when material is 282
- where $ is user defined output file name


# How to use as web interface

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
