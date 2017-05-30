FROM ubuntu:zesty

RUN echo "deb http://archive.ubuntu.com/ubuntu/ zesty main universe" >> /etc/apt/sources.list
RUN apt-get update -y
RUN apt-get install -y python3 python3-pip python-dev build-essential git postgresql python3-tk

COPY ./requirements.txt /app/requirements.txt
ADD https://s3-us-west-2.amazonaws.com/biof-hackathon-arvos/htsint.sql.gz /app/htsint.sql.gz
RUN gunzip -c /app/htsint.sql.gz > /app/htsint.sql
COPY ./dbconfig.log /root/.htsint/dbconfig.log
WORKDIR /app

RUN pip3 install -r requirements.txt

USER postgres

RUN /etc/init.d/postgresql start &&\
	psql --command "CREATE USER docker WITH SUPERUSER PASSWORD 'docker';" &&\
	psql --command "CREATE DATABASE docker WITH OWNER docker;" &&\
	psql --command "GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO docker;" &&\
	psql docker < htsint.sql
CMD ["/etc/init.d/postgresql", "start"]

USER root
COPY . /app

RUN git clone https://github.com/ajrichards/htsint.git /htsint
WORKDIR /htsint
RUN python3 setup.py install

EXPOSE 80
WORKDIR /app/web_app
CMD ["python3", "app.py"]
