FROM ubuntu

EXPOSE 8088 8888
RUN apt-get update
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install tornado
RUN pip install matplotlib
RUN pip install pandas
RUN pip install requests
RUN pip install numpy
RUN pip install scipy

ADD web_interface /web_interface
ADD data /data
ADD run.sh /run.sh

CMD /run.sh

