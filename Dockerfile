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

ADD templates/index.html /web_interface/templates/index.html
ADD api.py /web_interface/api.py
ADD ui_server_test.py /web_interface/ui_server_test.py
ADD creep.py /web_interface/creep.py
ADD diffusion.py /web_interface/diffusion.py
ADD data /web_interface/data
ADD run.sh /run.sh
RUN ["chmod", "+x", "/run.sh"]

CMD /run.sh

